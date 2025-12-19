"""
Enrichment analysis functionality for renalprog.

This module provides proper differential expression analysis using PyDESeq2
followed by Gene Set Enrichment Analysis (GSEA) for pathway analysis.

Key Components:
- PyDESeq2 differential expression analysis on RNA-seq count data
- GSEA (Gene Set Enrichment Analysis) execution and result processing
- Pathway enrichment heatmap generation

Pipeline Overview:
1. Convert log2(RSEM+1) data back to integer RSEM counts
2. Run PyDESeq2 to get statistically valid log2FoldChange values
3. Generate ranked gene lists (.rnk files) for GSEA
4. Execute GSEA in parallel for pathway enrichment
5. Combine and visualize results

This implementation follows the validated workflow from:
- D:\Repos\My_BRCA\src\fun_gsea\funs.py
- D:\Repos\My_BRCA\src_deseq_and_gsea_NCSR\py_deseq.py

IMPORTANT: This module uses PyDESeq2 for proper differential expression analysis.
DO NOT use simple log-fold change calculations on preprocessed data, as this
produces biologically invalid results.

Author: Renalprog Team (migrated from Guillermo Prol-Castelo)
Date: December 19, 2025
License: Apache 2.0
"""

import os
import glob
import re
import subprocess
import logging
import time
from pathlib import Path
from typing import List, Dict, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    import matplotlib.figure
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import warnings

import pandas as pd
import numpy as np

# Import PyDESeq2 for proper differential expression analysis
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from renalprog.config import PATHS

# Suppress warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", category=RuntimeWarning, message="^invalid value encountered in log")

# Set environment variable for NUMEXPR to allow more threads
os.environ['NUMEXPR_MAX_THREADS'] = '112'

logger = logging.getLogger(__name__)


class EnrichmentPipeline:
    """
    Main pipeline for dynamic enrichment analysis using PyDESeq2 and GSEA.

    This class orchestrates the complete enrichment analysis workflow:

    1. **PyDESeq2 Differential Expression Analysis**
       - Converts log2(RSEM+1) data back to integer RSEM counts
       - Runs PyDESeq2 for each trajectory timepoint vs controls
       - Generates statistically valid log2FoldChange values

    2. **GSEA Pathway Enrichment**
       - Creates ranked gene lists (.rnk files) from DESeq2 results
       - Executes GSEA in parallel with ReactomePathways gene sets
       - Collects pathway enrichment scores (NES, p-values, FDR)

    3. **Result Processing and Visualization**
       - Combines GSEA results across all trajectories and timepoints
       - Generates pathway enrichment heatmaps

    IMPORTANT: This pipeline uses PyDESeq2 for proper differential expression.
    DO NOT bypass this with simple fold-change calculations.
    """

    def __init__(
        self,
        trajectory_dir: str,
        output_dir: str,
        cancer_type: str = 'kirc',
        data_dir: Optional[str] = None,
        metadata_dir: Optional[str] = None,
        control_data_dir: Optional[str] = None,
        control_metadata_dir: Optional[str] = None,
        gsea_path: str = './GSEA_4.3.2/gsea-cli.sh',
        pathways_file: str = 'data/external/ReactomePathways.gmt',
        n_threads: int = 4
    ):
        """
        Initialize enrichment pipeline.

        Args:
            trajectory_dir: Directory containing trajectory CSV files
            output_dir: Output directory for enrichment results
            cancer_type: Cancer type ('kirc', 'lobular', 'ductal')
            data_dir: Path to preprocessed RNA-seq data
            metadata_dir: Path to clinical metadata
            control_data_dir: Path to control RNA-seq data
            control_metadata_dir: Path to control metadata
            gsea_path: Path to GSEA CLI tool
            pathways_file: Path to pathways GMT file
            n_threads: Number of parallel threads
        """
        self.trajectory_dir = Path(trajectory_dir)
        self.output_dir = Path(output_dir)
        self.cancer_type = cancer_type
        self.n_threads = n_threads
        self.gsea_path = Path(gsea_path)
        self.pathways_file = Path(pathways_file)

        # Set default data paths if not provided
        if data_dir is None:
            # Find latest preprocessed data
            data_dir = self._find_latest_preprocessed_data()
        if metadata_dir is None:
            metadata_dir = self._find_latest_preprocessed_metadata()
        if control_data_dir is None:
            control_data_dir = PATHS['processed'] / 'controls' / 'KIRC' / 'rnaseq_control.csv'
        if control_metadata_dir is None:
            control_metadata_dir = PATHS['processed'] / 'controls' / 'KIRC' / 'clinical_control.csv'

        self.data_dir = Path(data_dir)
        self.metadata_dir = Path(metadata_dir)
        self.control_data_dir = Path(control_data_dir)
        self.control_metadata_dir = Path(control_metadata_dir)

        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.deseq_dir = self.output_dir / 'deseq'
        self.gsea_dir = self.output_dir / 'gsea'
        self.deseq_dir.mkdir(exist_ok=True)
        self.gsea_dir.mkdir(exist_ok=True)

        logger.info(f"Initialized EnrichmentPipeline for {cancer_type}")
        logger.info(f"  Trajectory dir: {self.trajectory_dir}")
        logger.info(f"  Output dir: {self.output_dir}")
        logger.info(f"  Threads: {self.n_threads}")

    def _find_latest_preprocessed_data(self) -> Path:
        """Find the latest preprocessed RNA-seq data."""
        pattern = str(PATHS['interim'] / '*_preprocessed_KIRC' / 'rnaseq_maha.csv')
        files = sorted(glob.glob(pattern), reverse=True)
        if files:
            return Path(files[0])
        # Fallback to default
        return PATHS['interim'] / '20240930_preprocessed_KIRC' / 'rnaseq_maha.csv'

    def _find_latest_preprocessed_metadata(self) -> Path:
        """Find the latest preprocessed metadata."""
        pattern = str(PATHS['interim'] / '*_preprocessed_KIRC' / 'CuratedClinicalData.csv')
        files = sorted(glob.glob(pattern), reverse=True)
        if files:
            return Path(files[0])
        # Fallback to default
        return PATHS['interim'] / '20240930_preprocessed_KIRC' / 'CuratedClinicalData.csv'

    def run(
        self,
        skip_deseq: bool = False,
        skip_gsea: bool = False,
        cleanup: bool = False
    ):
        """
        Run the complete enrichment pipeline.

        Args:
            skip_deseq: Skip DESeq processing (use if already completed)
            skip_gsea: Skip GSEA analysis (use if already completed)
            cleanup: Remove intermediate files after processing
        """
        logger.info("Starting enrichment analysis pipeline")

        # Step 1: Process trajectories for DESeq
        if not skip_deseq:
            logger.info("Step 1/4: Processing trajectories for DESeq...")
            self._run_deseq_processing()
        else:
            logger.info("Step 1/4: Skipping DESeq processing")

        # Step 2: Run GSEA in parallel
        if not skip_gsea:
            logger.info("Step 2/4: Running GSEA analysis...")
            self._run_gsea_parallel()
        else:
            logger.info("Step 2/4: Skipping GSEA analysis")

        # Step 3: Combine GSEA results
        logger.info("Step 3/4: Combining GSEA results...")
        final_df = self._combine_gsea_results()

        # Save final dataset
        output_path = self.output_dir / 'trajectory_enrichment.csv'
        final_df.to_csv(output_path, index=False)
        logger.info(f"Final enrichment dataset saved to: {output_path}")
        logger.info(f"Dataset shape: {final_df.shape}")

        # Step 4: Generate pathway enrichment heatmap
        logger.info("Step 4/4: Generating pathway enrichment heatmap...")
        try:
            heatmap_data, heatmap_fig = generate_pathway_heatmap(
                enrichment_df=final_df,
                output_dir=self.output_dir
            )
            logger.info(f"Heatmap generated with {heatmap_data.shape[0]} significant pathways")
        except Exception as e:
            logger.error(f"Failed to generate heatmap: {e}", exc_info=True)
            logger.warning("Continuing without heatmap...")

        # Cleanup if requested
        if cleanup:
            logger.info("Cleaning up intermediate files...")
            self._cleanup()

        return final_df

    def _run_deseq_processing(self):
        """Process all trajectory files for DESeq analysis."""
        # Find all trajectory CSV files
        trajectory_files = list(self.trajectory_dir.glob('*.csv'))

        if not trajectory_files:
            raise ValueError(f"No trajectory files found in {self.trajectory_dir}")

        logger.info(f"Found {len(trajectory_files)} trajectory files")

        # Load data once
        rnaseq_data, clinical_data, control_data, control_metadata, gene_list = self._load_data()

        # Process files in parallel
        with ProcessPoolExecutor(max_workers=self.n_threads) as executor:
            futures = []
            for traj_file in trajectory_files:
                future = executor.submit(
                    process_trajectory_file,
                    traj_file=traj_file,
                    rnaseq_data=rnaseq_data,
                    clinical_data=clinical_data,
                    control_data=control_data,
                    control_metadata=control_metadata,
                    gene_list=gene_list,
                    output_dir=self.deseq_dir,
                    cancer_type=self.cancer_type,
                    gsea_path=self.gsea_path,
                    pathways_file=self.pathways_file
                )
                futures.append(future)

            # Track progress
            for future in tqdm(as_completed(futures), total=len(futures), desc="DESeq processing"):
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Error processing file: {e}")

        logger.info("DESeq processing complete")

        # Validate that files were created
        rnk_files = list(self.deseq_dir.rglob('*.rnk'))
        cmd_files = list(self.deseq_dir.glob('*.cmd'))

        logger.info(f"Created {len(rnk_files)} .rnk files")
        logger.info(f"Created {len(cmd_files)} .cmd files")

        if len(rnk_files) == 0:
            logger.error("No .rnk files were created during DESeq processing")
            raise ValueError("DESeq processing failed to create rank files")

        if len(cmd_files) == 0:
            logger.error("No .cmd files were created during DESeq processing")
            raise ValueError("DESeq processing failed to create command files")

    def _load_data(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, np.ndarray]:
        """Load preprocessed data, clinical data, and controls."""
        logger.info("Loading preprocessed data...")

        # Load RNA-seq data
        rnaseq_data = pd.read_csv(self.data_dir, index_col=0)

        # Load clinical data
        clinical_data = pd.read_csv(self.metadata_dir, index_col=0)
        clinical_data = pd.DataFrame(clinical_data['ajcc_pathologic_tumor_stage'])

        # For KIRC, convert stages to early/late
        if self.cancer_type == 'kirc':
            clinical_data.replace({
                'Stage I': 'early',
                'Stage II': 'early',
                'Stage III': 'late',
                'Stage IV': 'late'
            }, inplace=True)

        # Ensure correct shape
        if rnaseq_data.shape[1] != clinical_data.shape[0]:
            rnaseq_data = rnaseq_data.T

        # Load control data
        control_data = pd.read_csv(self.control_data_dir, index_col=0)
        control_metadata = pd.read_csv(self.control_metadata_dir, index_col=0)

        # Get gene list
        gene_list = rnaseq_data.index.values

        # Ensure control data has same genes (genes should be in index/rows)
        try:
            # If control_data has genes as columns, transpose it first
            if control_data.shape[0] != len(gene_list):
                control_data = control_data.T

            # Now select only the genes that are in gene_list, maintaining order
            control_data = control_data.loc[gene_list]
        except KeyError as e:
            raise ValueError(f"Control data missing some genes from preprocessed data: {e}")

        logger.info(f"Data loaded: RNA-seq shape={rnaseq_data.shape}, "
                   f"Clinical shape={clinical_data.shape}, "
                   f"Control shape={control_data.shape}")

        return rnaseq_data, clinical_data, control_data, control_metadata, gene_list

    def _run_gsea_parallel(self):
        """Run GSEA commands in parallel."""
        # Find all GSEA command files
        cmd_files = list(self.deseq_dir.glob('*.cmd'))

        if not cmd_files:
            raise ValueError(f"No GSEA command files found in {self.deseq_dir}")

        logger.info(f"Found {len(cmd_files)} GSEA command files")

        # Read all commands
        all_commands = []
        for cmd_file in cmd_files:
            with open(cmd_file, 'r') as f:
                commands = [line.strip() for line in f if line.strip()]
                all_commands.extend(commands)

        logger.info(f"Total GSEA commands to run: {len(all_commands)}")

        # Run commands in parallel
        failed_commands = []
        with ProcessPoolExecutor(max_workers=self.n_threads) as executor:
            futures = {}
            for cmd in all_commands:
                future = executor.submit(run_gsea_command, cmd)
                futures[future] = cmd

            # Track progress
            failed = 0
            for future in tqdm(as_completed(futures), total=len(futures), desc="GSEA analysis"):
                cmd = futures[future]
                try:
                    success = future.result()
                    if not success:
                        failed += 1
                        failed_commands.append(cmd)
                except Exception as e:
                    logger.error(f"Error running GSEA: {e}")
                    failed += 1
                    failed_commands.append(cmd)

        if failed > 0:
            logger.warning(f"{failed} GSEA commands failed")
            # Save failed commands for debugging
            failed_cmd_file = self.output_dir / 'failed_gsea_commands.txt'
            with open(failed_cmd_file, 'w') as f:
                f.write('\n'.join(failed_commands))
            logger.warning(f"Failed commands saved to: {failed_cmd_file}")

        logger.info("GSEA analysis complete")

        # Cleanup GSEA temp files
        self._cleanup_gsea_files()

    def _cleanup_gsea_files(self):
        """Remove unnecessary GSEA output files."""
        logger.info("Cleaning up GSEA temporary files...")

        # Remove HTML, PNG, CSS, RPT files
        for pattern in ['*.html', '*.png', '*.css', '*.rpt']:
            for file in self.deseq_dir.rglob(pattern):
                try:
                    file.unlink()
                except Exception as e:
                    logger.debug(f"Could not delete {file}: {e}")

        # Remove edb directories
        for edb_dir in self.deseq_dir.rglob('edb'):
            try:
                import shutil
                shutil.rmtree(edb_dir)
            except Exception as e:
                logger.debug(f"Could not delete {edb_dir}: {e}")

    def _combine_gsea_results(self) -> pd.DataFrame:
        """Combine all GSEA results into a single dataset."""
        # Load pathways
        pathways = load_pathways_from_gmt(self.pathways_file)

        # Diagnostic: show directory structure
        logger.info(f"Scanning GSEA results in: {self.deseq_dir}")
        logger.info("Directory structure:")
        for item in self.deseq_dir.iterdir():
            if item.is_dir():
                logger.info(f"  {item.name}/")
                for subitem in item.iterdir():
                    if subitem.is_dir():
                        gsea_dirs = list(subitem.glob('gsea_tp*'))
                        logger.info(f"    {subitem.name}/ ({len(gsea_dirs)} gsea_tp* dirs)")
                        # Show deeper structure for debugging
                        for gsea_dir in gsea_dirs[:2]:  # Show first 2
                            subdirs = [d.name for d in gsea_dir.iterdir() if d.is_dir()]
                            files = [f.name for f in gsea_dir.iterdir() if f.is_file()]
                            logger.debug(f"      {gsea_dir.name}/ subdirs={subdirs[:3]}, files={files[:3]}")
            else:
                logger.info(f"  {item.name}")

        # Find all trajectory directories (transition/patient structure)
        # First, look for transition directories
        transition_dirs = [d for d in self.deseq_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]

        logger.info(f"Found {len(transition_dirs)} transition directories")

        all_results = []
        failed_trajectories = []
        skipped_trajectories = []

        for transition_dir in transition_dirs:
            # Find patient directories within each transition
            patient_dirs = [d for d in transition_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
            logger.info(f"  Found {len(patient_dirs)} patient directories in {transition_dir.name}")

            for patient_dir in tqdm(patient_dirs, desc=f"Processing {transition_dir.name}"):
                try:
                    result = process_trajectory_results(patient_dir, pathways)
                    if result is not None:
                        all_results.append(result)
                    else:
                        skipped_trajectories.append(str(patient_dir))
                        logger.warning(f"No results for {patient_dir}")
                except Exception as e:
                    failed_trajectories.append((str(patient_dir), str(e)))
                    logger.error(f"Error processing {patient_dir}: {e}", exc_info=True)

        # Report statistics
        logger.info(f"Successfully processed: {len(all_results)} trajectories")
        logger.info(f"Skipped (no GSEA dirs): {len(skipped_trajectories)} trajectories")
        logger.info(f"Failed (exceptions): {len(failed_trajectories)} trajectories")

        if failed_trajectories:
            logger.warning("Failed trajectories:")
            for path, error in failed_trajectories[:10]:  # Show first 10
                logger.warning(f"  {path}: {error}")

        if not all_results:
            logger.error(f"No GSEA results found in {self.deseq_dir}")
            logger.error("Directory structure should be: deseq_dir/transition/patient/gsea_tp*/")
            logger.error(f"Transition dirs found: {len(transition_dirs)}")
            logger.error(f"Total patient dirs scanned: {len(skipped_trajectories) + len(failed_trajectories)}")

            # Provide more specific error message
            if failed_trajectories:
                logger.error(f"All trajectories failed with errors. First error: {failed_trajectories[0][1]}")
            elif skipped_trajectories:
                logger.error("All trajectories were skipped (no gsea_tp* directories found)")

            raise ValueError("No GSEA results found to combine")

        # Combine all results
        final_df = pd.concat(all_results, axis=0, ignore_index=True)

        return final_df

    def _cleanup(self):
        """Remove intermediate files."""
        import shutil

        # Remove command files
        for cmd_file in self.deseq_dir.glob('*.cmd'):
            cmd_file.unlink()

        logger.info("Cleanup complete")


def process_trajectory_file(
    traj_file: Path,
    rnaseq_data: pd.DataFrame,
    clinical_data: pd.DataFrame,
    control_data: pd.DataFrame,
    control_metadata: pd.DataFrame,
    gene_list: np.ndarray,
    output_dir: Path,
    cancer_type: str,
    gsea_path: Path,
    pathways_file: Path
) -> None:
    """
    Process a single trajectory file for DESeq2 differential expression analysis.

    This function:
    1. Loads the synthetic trajectory
    2. For each timepoint, runs PyDESeq2 analysis vs controls
    3. Generates .rnk rank files and GSEA commands

    Args:
        traj_file: Path to trajectory CSV file
        rnaseq_data: Preprocessed RNA-seq data (in log2(RSEM+1) format)
        clinical_data: Clinical metadata
        control_data: Control RNA-seq data (in log2(RSEM+1) format)
        control_metadata: Control metadata
        gene_list: List of genes
        output_dir: Output directory
        cancer_type: Cancer type
        gsea_path: Path to GSEA CLI
        pathways_file: Path to pathways GMT file
    """
    # Load trajectory
    trajectory = pd.read_csv(traj_file, index_col=0).T

    # Get patient and transition info from filename
    patient_name = traj_file.stem
    transition = traj_file.parent.name

    # Create output subdirectory
    traj_output_dir = output_dir / transition / patient_name
    traj_output_dir.mkdir(parents=True, exist_ok=True)

    # Create reports directory for GSEA outputs
    reports_dir = traj_output_dir / 'reports'
    reports_dir.mkdir(exist_ok=True)

    # Process each timepoint
    gsea_commands = []

    for interpol_idx in range(trajectory.shape[1]):
        # Get timepoint data
        timepoint_data = pd.DataFrame(trajectory.iloc[:, interpol_idx])

        # Create unique sample name
        sample_name = f'{patient_name}_{interpol_idx}'

        # Run PyDESeq2 analysis to get proper log2FoldChange values
        try:
            log2fc = run_deseq2_analysis(
                sample_data=timepoint_data,
                control_data=control_data,
                control_metadata=control_metadata,
                gene_list=gene_list,
                sample_name=sample_name,
                stage_transition=transition
            )
        except Exception as e:
            logger.error(f"DESeq2 analysis failed for {sample_name}: {e}")
            logger.debug(f"  Timepoint data shape: {timepoint_data.shape}")
            logger.debug(f"  Control data shape: {control_data.shape}")
            logger.debug(f"  Gene list length: {len(gene_list)}")
            continue

        # Save rank file (.rnk format for GSEA)
        # Format: gene_name\tlog2FoldChange (tab-separated, no header)
        rnk_file = traj_output_dir / f'{patient_name}_{interpol_idx}.rnk'

        # Write with proper format for GSEA (gene name, tab, value)
        with open(rnk_file, 'w') as f:
            f.write('#\n')  # GSEA ignores first line starting with #
            for gene, fc_value in log2fc.items():
                f.write(f'{gene}\t{fc_value}\n')

        # Generate GSEA command
        gsea_output = reports_dir

        cmd = generate_gsea_command(
            gsea_path=gsea_path,
            rnk_file=rnk_file,
            gmt_file=pathways_file,
            output_dir=gsea_output,
            label=f'{patient_name}_{interpol_idx}'
        )

        gsea_commands.append(cmd)

    # Save GSEA commands to file
    cmd_file = output_dir / f'gsea_commands_{transition}_{patient_name}.cmd'
    with open(cmd_file, 'w') as f:
        f.write('\n'.join(gsea_commands))


def run_deseq2_analysis(
    sample_data: pd.DataFrame,
    control_data: pd.DataFrame,
    control_metadata: pd.DataFrame,
    gene_list: np.ndarray,
    sample_name: str,
    stage_transition: str
) -> pd.Series:
    """
    Perform DESeq2 differential expression analysis between sample and controls.

    This function properly:
    1. Converts log2(RSEM+1) data back to RSEM integer counts
    2. Runs PyDESeq2 analysis to get statistically valid log2FoldChange values
    3. Returns ranked gene list for GSEA

    Args:
        sample_data: Sample expression data (genes x 1) in log2(RSEM+1) format
        control_data: Control expression data (genes x samples) in log2(RSEM+1) format
        control_metadata: Control clinical metadata with stage information
        gene_list: List of genes
        sample_name: Name/ID of the sample
        stage_transition: Stage transition label (e.g., 'early_to_late', 'I_to_II')

    Returns:
        Series of log2FoldChange values sorted for GSEA input
    """
    # Debug logging
    logger.debug(f"DESeq2 analysis for {sample_name}")
    logger.debug(f"  Sample data shape: {sample_data.shape}")
    logger.debug(f"  Control data shape: {control_data.shape}")
    logger.debug(f"  Gene list length: {len(gene_list)}")

    # 1. Convert from log2(RSEM+1) back to RSEM counts
    # The preprocessed data is in log2(RSEM+1) format, so we reverse this transformation
    # and round to integers as required by DESeq2
    sample_counts = 2 ** sample_data.values.flatten() - 1
    control_counts = 2 ** control_data.values - 1

    # Clip negative values to 0 (can occur due to numerical precision or interpolation)
    sample_counts = np.clip(sample_counts, 0, None)
    control_counts = np.clip(control_counts, 0, None)

    # Round to integers as required by DESeq2
    sample_counts = np.round(sample_counts).astype(int)
    control_counts = np.round(control_counts).astype(int)

    # Ensure control_counts is 2D (genes x samples)
    if control_counts.ndim == 1:
        control_counts = control_counts[:, np.newaxis]

    # Validate counts
    if np.any(sample_counts < 0):
        logger.error(f"Sample has negative counts after conversion: min={sample_counts.min()}")
        raise ValueError("Sample counts contain negative values after conversion")
    if np.any(control_counts < 0):
        logger.error(f"Controls have negative counts after conversion: min={control_counts.min()}")
        raise ValueError("Control counts contain negative values after conversion")

    logger.debug(f"  Sample counts: min={sample_counts.min()}, max={sample_counts.max()}, mean={sample_counts.mean():.1f}")
    logger.debug(f"  Control counts: min={control_counts.min()}, max={control_counts.max()}, mean={control_counts.mean():.1f}")

    # Ensure sample_counts has the same number of genes as control_counts
    if len(sample_counts) != control_counts.shape[0]:
        raise ValueError(
            f"Sample has {len(sample_counts)} genes but controls have {control_counts.shape[0]} genes. "
            f"Sample shape: {sample_data.shape}, Control shape: {control_data.shape}"
        )

    # 2. Combine sample with controls
    # Create count matrix: genes (rows) x samples (columns)
    # sample_counts is 1D (n_genes,), control_counts is 2D (n_genes, n_controls)
    # We need to reshape sample_counts to (n_genes, 1) to stack horizontally
    counts_matrix = np.column_stack([sample_counts.reshape(-1, 1), control_counts])

    counts_df = pd.DataFrame(
        counts_matrix,
        index=gene_list,
        columns=[sample_name] + list(control_data.columns)
    )

    # 3. Create metadata DataFrame with condition labels
    metadata_df = pd.DataFrame({
        'condition': [stage_transition] + list(control_metadata['ajcc_pathologic_tumor_stage'].values)
    }, index=counts_df.columns)

    # 4. Run PyDESeq2 analysis
    # Initialize DESeqDataSet
    dds = DeseqDataSet(
        counts=counts_df.T,  # DESeq2 expects samples as rows, genes as columns
        metadata=metadata_df,
        design_factors='condition',
        refit_cooks=True,
        quiet=True  # Suppress verbose output in parallel processing
    )

    # Fit dispersions and log-fold changes
    dds.deseq2()

    # Statistical analysis
    stat_res = DeseqStats(
        dds,
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
        quiet=True
    )
    stat_res.summary()

    # 5. Extract log2FoldChange values
    results_df = stat_res.results_df

    # Ensure we have log2FoldChange column
    if 'log2FoldChange' not in results_df.columns:
        raise ValueError(f"DESeq2 results missing 'log2FoldChange' column. Available columns: {results_df.columns.tolist()}")

    # Extract log2FoldChange and sort by absolute value (for GSEA ranking)
    log2fc = results_df['log2FoldChange'].copy()

    # Replace NaN values with 0 (genes with no change or insufficient data)
    log2fc = log2fc.fillna(0)

    # Sort by absolute value descending (GSEA expects ranked list)
    log2fc = log2fc.sort_values(ascending=False)

    logger.debug(f"  DESeq2 complete: {len(log2fc)} genes, range [{log2fc.min():.2f}, {log2fc.max():.2f}]")

    return log2fc


def generate_gsea_command(
    gsea_path: Path,
    rnk_file: Path,
    gmt_file: Path,
    output_dir: Path,
    label: Optional[str] = None
) -> str:
    """
    Generate GSEA command line without quotes around paths (standard GSEA CLI format).

    Args:
        gsea_path: Path to GSEA CLI script
        rnk_file: Path to ranked gene list file
        gmt_file: Path to gene sets GMT file
        output_dir: Output directory for GSEA results
        label: Custom report label (if None, uses rnk filename stem)

    Returns:
        GSEA command string
    """
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Convert paths to absolute paths to avoid issues
    gsea_path = gsea_path.resolve()
    rnk_file = rnk_file.resolve()
    gmt_file = gmt_file.resolve()
    output_dir = output_dir.resolve()

    # Generate label from RNK filename if not provided
    if label is None:
        label = rnk_file.stem.replace('_foldchange', '')

    # Build command - no quotes around paths, as per standard GSEA CLI usage
    cmd = (
        f'{gsea_path} GSEAPreranked '
        f'-rnk {rnk_file} '
        f'-gmx {gmt_file} '
        f'-collapse No_Collapse '
        f'-mode Max_probe '
        f'-norm meandiv '
        f'-nperm 1000 '
        f'-rnd_seed timestamp '
        f'-scoring_scheme weighted '
        f'-rpt_label {label} '
        f'-create_svgs false '
        f'-include_only_symbols true '
        f'-make_sets false '
        f'-plot_top_x 20 '
        f'-set_max 500 '
        f'-set_min 15 '
        f'-zip_report false '
        f'-out {output_dir}'
    )

    return cmd


def run_gsea_command(cmd: str) -> bool:
    """
    Run a single GSEA command.

    Args:
        cmd: GSEA command string

    Returns:
        True if successful, False otherwise
    """
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout per command
        )

        if result.returncode != 0:
            logger.error(f"GSEA command failed with return code {result.returncode}")
            logger.error(f"Command: {cmd[:200]}...")  # Truncate long commands
            logger.error(f"STDERR: {result.stderr[:500]}")  # First 500 chars
            if result.stdout:
                logger.error(f"STDOUT: {result.stdout[:500]}")
            return False

        # Check for Java errors in stdout (GSEA sometimes returns 0 even on error)
        if result.stdout and ("error" in result.stdout.lower() or "exception" in result.stdout.lower()):
            logger.warning(f"GSEA command may have failed (found error in output)")
            logger.warning(f"Command: {cmd[:200]}...")
            logger.warning(f"Output: {result.stdout[:500]}")  # First 500 chars
            return False

        # Check for GSEA-specific success indicators
        if result.stdout and "Enrichment score" in result.stdout:
            logger.debug(f"GSEA command completed successfully")
            return True
        elif result.stdout:
            # Log a sample of output for debugging
            logger.debug(f"GSEA output sample: {result.stdout[:200]}")

        return True

    except subprocess.TimeoutExpired:
        logger.error(f"GSEA command timed out after 600s")
        logger.error(f"Command: {cmd[:200]}...")
        return False
    except Exception as e:
        logger.error(f"Error running GSEA command: {e}")
        logger.error(f"Command: {cmd[:200]}...")
        return False


def process_static_patient(
    patient_id: str,
    patient_stage: str,
    rnaseq_data: pd.DataFrame,
    control_data: pd.DataFrame,
    control_metadata: pd.DataFrame,
    gene_list: np.ndarray,
    output_dir: Path,
    gsea_path: Path,
    pathways_file: Path
) -> Optional[str]:
    """
    Process a single real patient for static enrichment analysis.

    This function runs PyDESeq2 on a single patient vs controls and generates
    GSEA command for pathway enrichment.

    Args:
        patient_id: Patient identifier
        patient_stage: Patient cancer stage
        rnaseq_data: Patient RNA-seq data (in log2(RSEM+1) format)
        control_data: Control RNA-seq data (in log2(RSEM+1) format)
        control_metadata: Control metadata
        gene_list: List of genes
        output_dir: Output directory
        gsea_path: Path to GSEA CLI
        pathways_file: Path to pathways GMT file

    Returns:
        GSEA command string if successful, None otherwise
    """
    # Create output directory for this stage
    stage_dir = output_dir / patient_stage
    stage_dir.mkdir(parents=True, exist_ok=True)

    patient_dir = stage_dir / patient_id
    patient_dir.mkdir(exist_ok=True)

    reports_dir = patient_dir / 'reports'
    reports_dir.mkdir(exist_ok=True)

    # Get patient data
    patient_data = pd.DataFrame(rnaseq_data.loc[patient_id])

    # Run PyDESeq2 analysis
    try:
        log2fc = run_deseq2_analysis(
            sample_data=patient_data,
            control_data=control_data,
            control_metadata=control_metadata,
            gene_list=gene_list,
            sample_name=patient_id,
            stage_transition=patient_stage
        )
    except Exception as e:
        logger.error(f"DESeq2 analysis failed for patient {patient_id}: {e}")
        return None

    # Save rank file
    rnk_file = patient_dir / f'{patient_id}.rnk'
    with open(rnk_file, 'w') as f:
        f.write('#\n')
        for gene, fc_value in log2fc.items():
            f.write(f'{gene}\t{fc_value}\n')

    # Generate GSEA command
    cmd = generate_gsea_command(
        gsea_path=gsea_path,
        rnk_file=rnk_file,
        gmt_file=pathways_file,
        output_dir=reports_dir,
        label=patient_id
    )

    return cmd


def load_pathways_from_gmt(gmt_file: Path) -> List[str]:
    """
    Load pathway names from GMT file.

    Args:
        gmt_file: Path to GMT file

    Returns:
        List of pathway names
    """
    pathways = []

    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                pathway_name = parts[0]
                pathways.append(pathway_name)

    return pathways


def process_trajectory_results(
    traj_dir: Path,
    pathways: List[str]
) -> Optional[pd.DataFrame]:
    """
    Process GSEA results for a single trajectory.

    Args:
        traj_dir: Directory containing GSEA results for one trajectory
        pathways: List of all pathway names

    Returns:
        DataFrame with combined results for this trajectory
    """
    # Get patient and transition info
    patient_name = traj_dir.name
    transition = traj_dir.parent.name

    # Find all GSEA output directories
    gsea_dirs = sorted(traj_dir.glob('gsea_tp*'))

    if not gsea_dirs:
        logger.warning(f"No GSEA results found in {traj_dir}")
        return None

    results = []

    for idx, gsea_dir in enumerate(gsea_dirs):
        try:
            # Read GSEA reports
            result = read_gsea_reports(
                gsea_dir=gsea_dir,
                patient_name=patient_name,
                interpol_index=idx,
                transition=transition,
                pathways=pathways
            )
            results.append(result)
        except Exception as e:
            logger.error(f"Error reading GSEA results from {gsea_dir}: {e}")

    if not results:
        return None

    # Combine all timepoints
    combined = pd.concat(results, axis=0, ignore_index=True)

    return combined


def read_gsea_reports(
    gsea_dir: Path,
    patient_name: str,
    interpol_index: int,
    transition: str,
    pathways: List[str]
) -> pd.DataFrame:
    """
    Read GSEA positive and negative reports and combine them.

    Args:
        gsea_dir: Directory containing GSEA output
        patient_name: Patient identifier
        interpol_index: Interpolation index (timepoint)
        transition: Transition label (e.g., 'early_to_late')
        pathways: List of all pathway names

    Returns:
        DataFrame with columns: Patient, Idx, Transition, NAME, ES, NES, FDR q-val
    """
    # Find the actual GSEA output directory (it has a timestamp suffix)
    gsea_subdirs = [d for d in gsea_dir.iterdir() if d.is_dir()]

    if not gsea_subdirs:
        logger.warning(f"No GSEA output subdirectory found in {gsea_dir}")
        # List what's actually there
        contents = list(gsea_dir.iterdir())
        logger.warning(f"Contents of {gsea_dir}: {[c.name for c in contents]}")
        raise FileNotFoundError(f"No GSEA output subdirectory in {gsea_dir}")

    # Use the first (should be only) subdirectory
    gsea_output = gsea_subdirs[0]
    logger.debug(f"Reading GSEA results from {gsea_output}")

    # Find report files
    pos_files = list(gsea_output.glob('gsea_report_for_na_pos*.tsv'))
    neg_files = list(gsea_output.glob('gsea_report_for_na_neg*.tsv'))

    if not pos_files:
        logger.error(f"No positive report file found in {gsea_output}")
        # List all TSV files to help diagnose
        all_tsvs = list(gsea_output.glob('*.tsv'))
        all_files = list(gsea_output.glob('*'))
        logger.error(f"TSV files found: {[f.name for f in all_tsvs]}")
        logger.error(f"All files: {[f.name for f in all_files[:20]]}")  # First 20
        raise FileNotFoundError(f"GSEA positive report not found in {gsea_output}")

    if not neg_files:
        logger.error(f"No negative report file found in {gsea_output}")
        # List all TSV files to help diagnose
        all_tsvs = list(gsea_output.glob('*.tsv'))
        logger.error(f"TSV files found: {[f.name for f in all_tsvs]}")
        raise FileNotFoundError(f"GSEA negative report not found in {gsea_output}")

    # Read reports
    try:
        pos_report = pd.read_csv(pos_files[0], sep='\t')
        neg_report = pd.read_csv(neg_files[0], sep='\t')
    except Exception as e:
        logger.error(f"Error reading report files: {e}")
        logger.error(f"Positive file: {pos_files[0]}")
        logger.error(f"Negative file: {neg_files[0]}")
        raise

    # Check if reports have expected columns
    required_cols = ['NAME', 'ES', 'NES', 'FDR q-val']
    for col in required_cols:
        if col not in pos_report.columns:
            logger.error(f"Missing column '{col}' in positive report")
            logger.error(f"Available columns: {pos_report.columns.tolist()}")
            raise ValueError(f"Missing required column '{col}' in GSEA positive report")
        if col not in neg_report.columns:
            logger.error(f"Missing column '{col}' in negative report")
            logger.error(f"Available columns: {neg_report.columns.tolist()}")
            raise ValueError(f"Missing required column '{col}' in GSEA negative report")

    # Keep only relevant columns
    pos_report = pos_report[required_cols]
    neg_report = neg_report[required_cols]

    # Combine
    combined = pd.concat([pos_report, neg_report], axis=0, ignore_index=True)

    # Add metadata columns
    combined.insert(0, 'Patient', patient_name)
    combined.insert(1, 'Idx', interpol_index)
    combined.insert(2, 'Transition', transition)

    # Add missing pathways
    combined = add_missing_pathways(
        patient=patient_name,
        transition=transition,
        idx=interpol_index,
        report=combined,
        pathways=pathways
    )

    return combined


def add_missing_pathways(
    patient: str,
    transition: str,
    idx: int,
    report: pd.DataFrame,
    pathways: List[str]
) -> pd.DataFrame:
    """
    Add missing pathways to report with NaN values.

    Args:
        patient: Patient identifier
        transition: Transition label
        idx: Interpolation index
        report: GSEA report DataFrame
        pathways: List of all pathway names

    Returns:
        DataFrame including missing pathways
    """
    # Find missing pathways
    missing = [p for p in pathways if p not in report['NAME'].values]

    if not missing:
        return report

    # Create DataFrame for missing pathways
    missing_df = pd.DataFrame({
        'Patient': patient,
        'Idx': idx,
        'Transition': transition,
        'NAME': missing,
        'ES': np.nan,
        'NES': np.nan,
        'FDR q-val': np.nan
    })

    # Combine
    combined = pd.concat([report, missing_df], axis=0, ignore_index=True)

    return combined


# Convenience functions for use in other scripts
def process_trajectories_for_deseq(
    trajectory_dir: str,
    output_dir: str,
    cancer_type: str = 'kirc',
    n_threads: int = 4
) -> None:
    """
    Convenience function to process trajectories for DESeq.

    Args:
        trajectory_dir: Directory containing trajectory CSV files
        output_dir: Output directory
        cancer_type: Cancer type
        n_threads: Number of threads
    """
    pipeline = EnrichmentPipeline(
        trajectory_dir=trajectory_dir,
        output_dir=output_dir,
        cancer_type=cancer_type,
        n_threads=n_threads
    )
    pipeline._run_deseq_processing()


def run_gsea_parallel(
    deseq_dir: str,
    n_threads: int = 4
) -> None:
    """
    Convenience function to run GSEA in parallel.

    Args:
        deseq_dir: Directory containing DESeq results and GSEA commands
        n_threads: Number of threads
    """
    pipeline = EnrichmentPipeline(
        trajectory_dir='.',  # Not used
        output_dir=Path(deseq_dir).parent,
        n_threads=n_threads
    )
    pipeline.deseq_dir = Path(deseq_dir)
    pipeline._run_gsea_parallel()


def combine_gsea_results(
    deseq_dir: str,
    pathways_file: str = 'data/external/ReactomePathways.gmt'
) -> pd.DataFrame:
    """
    Convenience function to combine GSEA results.

    Args:
        deseq_dir: Directory containing GSEA results
        pathways_file: Path to pathways GMT file

    Returns:
        Combined enrichment DataFrame
    """
    pipeline = EnrichmentPipeline(
        trajectory_dir='.',  # Not used
        output_dir=Path(deseq_dir).parent,
        pathways_file=pathways_file
    )
    pipeline.deseq_dir = Path(deseq_dir)
    return pipeline._combine_gsea_results()


def generate_pathway_heatmap(
    enrichment_df: pd.DataFrame,
    output_dir: str,
    fdr_threshold: float = 0.05,
    colorbar: bool = True,
    legend: bool = False,
    yticks_fontsize: int = 12,
    show: bool = False
) -> Tuple[pd.DataFrame, Dict[str, 'matplotlib.figure.Figure']]:
    """
    Generate multiple pathway enrichment heatmaps from GSEA results.

    This function creates several heatmaps showing the sum of NES (Normalized Enrichment Score)
    across all trajectories for each pathway at each timepoint:

    1. Top 50 most changing pathways (first vs last timepoint)
    2. Top 50 most upregulated pathways (average NES > 0)
    3. Top 50 most downregulated pathways (average NES < 0)
    4. Selected pathways (high-level Reactome + literature pathways)

    The heatmaps have:
    - Rows: Pathway names
    - Columns: Timepoints (pseudo-time from early to late)
    - Values: Sum of NES across all trajectories at each timepoint

    Args:
        enrichment_df: DataFrame with columns [Patient, Idx, Transition, NAME, ES, NES, FDR q-val]
        output_dir: Output directory for heatmap files
        fdr_threshold: FDR q-value threshold for significance (default: 0.05)
        colorbar: Whether to show colorbar (default: True)
        legend: Whether to show legend (default: False)
        yticks_fontsize: Font size for y-axis tick labels (default: 12)
        show: Whether to display the plot (default: False)

    Returns:
        Tuple of (heatmap_data, figures_dict):
            - heatmap_data: DataFrame with summed NES values (pathways × timepoints)
            - figures_dict: Dictionary mapping figure names to Matplotlib Figure objects

    Example:
        >>> enrichment_df = pd.read_csv('trajectory_enrichment.csv')
        >>> heatmap_data, figs = generate_pathway_heatmap(
        ...     enrichment_df=enrichment_df,
        ...     output_dir='results/',
        ...     fdr_threshold=0.05
        ... )
        >>> print(f"Generated {len(figs)} heatmaps")
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.colors as mcolors

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Generating pathway enrichment heatmaps (FDR < {fdr_threshold})...")

    # Define pathway lists
    highest_pathways = [
        "Autophagy",
        "Cell Cycle",
        "Cell-Cell communication",
        "Cellular responses to stimuli",
        "Chromatin organization",
        "Circadian Clock",
        "DNA Repair",
        "DNA Replication",
        "Developmental Biology",
        "Digestion and absorption",
        "Disease",
        "Drug ADME",
        "Extracellular matrix organization",
        "Gene expression (Transcription)",
        "Hemostasis",
        "Immune System",
        "Metabolism",
        "Metabolism of RNA",
        "Metabolism of proteins",
        "Muscle contraction",
        "Neuronal System",
        "Organelle biogenesis and maintenance",
        "Programmed Cell Death",
        "Protein localization",
        "Reproduction",
        "Sensory Perception",
        "Signal Transduction",
        "Transport of small molecules",
        "Vesicle-mediated transport",
    ]

    pathways_literature = [
        # VHL/HIF pathway
        "CELLULAR RESPONSE TO HYPOXIA",
        "OXYGEN-DEPENDENT PROLINE HYDROXYLATION OF HYPOXIA-INDUCIBLE FACTOR ALPHA",
        "REGULATION OF GENE EXPRESSION BY HYPOXIA-INDUCIBLE FACTOR",
        # PI3K/AKT/MTOR Pathway
        "PI3K/AKT ACTIVATION",
        "PI3K/AKT SIGNALING IN CANCER",
        "MTOR SIGNALLING",
        # Warburg effect
        "TP53 REGULATES METABOLIC GENES",
        "GLYCOLYSIS",
        "GLUCOSE METABOLISM",
        # TCA/Krebs cycle
        "CITRIC ACID CYCLE (TCA CYCLE)",
        "THE CITRIC ACID (TCA) CYCLE AND RESPIRATORY ELECTRON TRANSPORT",
        # Pentose phosphate pathway
        "NFE2L2 REGULATES PENTOSE PHOSPHATE PATHWAY GENES",
        "PENTOSE PHOSPHATE PATHWAY",
        "PENTOSE PHOSPHATE PATHWAY DISEASE",
        # Fatty Acid Metabolism
        "FATTY ACID METABOLISM",
        # Glutamine metabolism
        "GLUTAMATE AND GLUTAMINE METABOLISM",
        # EGFR
        "SIGNALING BY EGFR",
        "SIGNALING BY EGFR IN CANCER",
        "EGFR DOWNREGULATION",
        # TGF-β signaling
        "SIGNALING BY TGF-BETA RECEPTOR COMPLEX",
        "TGF-BETA RECEPTOR SIGNALING IN EMT (EPITHELIAL TO MESENCHYMAL TRANSITION)",
        "SIGNALING BY TGF-BETA RECEPTOR COMPLEX IN CANCER",
        "SIGNALING BY TGFB FAMILY MEMBERS",
        "TGF-BETA RECEPTOR SIGNALING ACTIVATES SMADS",
        # Wnt/β-catenin pathway
        "BETA-CATENIN INDEPENDENT WNT SIGNALING",
        "SIGNALING BY WNT",
        # SLIT-2-ROBO1 pathways
        "REGULATION OF EXPRESSION OF SLITS AND ROBOS",
        # DNA repair
        "DNA REPAIR",
        # Energy homeostasis
        "ION HOMEOSTASIS",
        # Apoptosis
        "APOPTOSIS",
        # Angiogenesis
        "SIGNALING BY VEGF"
    ]

    # Step 1: Ensure numeric types for NES and FDR q-val
    enrichment_df = enrichment_df.copy()
    enrichment_df['NES'] = pd.to_numeric(enrichment_df['NES'], errors='coerce')
    enrichment_df['FDR q-val'] = pd.to_numeric(enrichment_df['FDR q-val'], errors='coerce')
    enrichment_df['Idx'] = pd.to_numeric(enrichment_df['Idx'], errors='coerce')

    # Log any rows that had non-numeric values
    invalid_nes = enrichment_df['NES'].isna().sum()
    invalid_fdr = enrichment_df['FDR q-val'].isna().sum()
    if invalid_nes > 0:
        logger.warning(f"Found {invalid_nes} rows with non-numeric NES values (converted to NaN)")
    if invalid_fdr > 0:
        logger.warning(f"Found {invalid_fdr} rows with non-numeric FDR q-val values (converted to NaN)")

    # Step 2: Filter by FDR threshold
    significant = enrichment_df[enrichment_df['FDR q-val'] < fdr_threshold].copy()

    logger.info(f"Found {significant.shape[0]} significant pathway enrichments (FDR < {fdr_threshold})")

    if significant.empty:
        logger.warning("No significant pathways found. Cannot generate heatmap.")
        # Return empty results
        empty_df = pd.DataFrame()
        empty_dict = {}
        return empty_df, empty_dict

    # Step 3: Group by Timepoint (Idx) and Pathway (NAME), sum NES across all trajectories
    pathway_summary = significant.groupby(['Idx', 'NAME'])['NES'].sum().reset_index()

    logger.info(f"Aggregated results for {pathway_summary['Idx'].nunique()} timepoints "
                f"and {pathway_summary['NAME'].nunique()} pathways")

    # Step 4: Pivot to create matrix (pathways × timepoints)
    heatmap_data = pathway_summary.pivot(
        index='NAME',
        columns='Idx',
        values='NES'
    ).fillna(0)  # Fill missing with 0

    logger.info(f"Full heatmap dimensions: {heatmap_data.shape[0]} pathways × {heatmap_data.shape[1]} timepoints")

    # Save full summary data
    summary_file = output_dir / 'pathway_nes_summary.csv'
    heatmap_data.to_csv(summary_file)
    logger.info(f"Saved pathway NES summary to: {summary_file}")

    # Dictionary to store all figures
    figures = {}

    # Helper function to create and save a heatmap
    def plot_heatmap_regulation(df_plot, unique_pathways, cmap_here='viridis',
                               save_name=None, colorbar_title='Sum of NES'):
        """Plot heatmap following paper_figures.ipynb style"""
        # Generate a range of locations for the ticks
        tick_locations = range(len(unique_pathways))
        z_min, z_max = df_plot.min().min(), df_plot.max().max()

        # Make the range symmetric around 0
        if z_min < 0 and z_max > 0:
            abs_max = max(abs(z_min), abs(z_max))
            z_min, z_max = -abs_max, abs_max
            norm = mcolors.TwoSlopeNorm(vmin=z_min, vcenter=0, vmax=z_max)
        else:
            # If all values are positive or all negative, use regular normalization
            norm = mcolors.Normalize(vmin=z_min, vmax=z_max)
            logger.warning(f"Cannot center colormap at zero: range [{z_min:.3f}, {z_max:.3f}] does not cross zero")

        fig, ax = plt.subplots(figsize=(30, 10))

        # Make heatmap
        cax = ax.imshow(df_plot.values, cmap=cmap_here, norm=norm, aspect='auto')

        # Set the y-ticks
        plt.yticks(tick_locations, unique_pathways, fontsize=yticks_fontsize)

        # Set the x-ticks at specific positions
        num_timepoints = df_plot.shape[1]
        ax.set_xticks([0, num_timepoints - 1])
        ax.set_xticklabels(['early', 'late'], fontsize=yticks_fontsize*1.33, rotation=45)
        ax.set_xlabel('Pseudo-Time', fontsize=yticks_fontsize*1.33)

        # Get x and y axis range
        ymin, ymax = ax.get_ylim()
        xmin, xmax = ax.get_xlim()

        # Custom x and y ticks
        x_custom = np.arange(xmin, xmax, step=1)
        y_custom = np.arange(ymax, ymin, step=1)

        # set minor ticks at custom locations:
        ax.set_xticks(x_custom, minor=True)
        ax.set_yticks(y_custom, minor=True)

        # Add grid lines at both major and minor ticks
        plt.grid(False, which='major')
        plt.grid(True, which='minor', color='black', linestyle='-', linewidth=1)

        # Remove ticks
        ax.tick_params(axis='both', which='minor', length=0)

        # Add colorbar
        if colorbar:
            cbar = plt.colorbar(cax, shrink=0.7)
            cbar.set_label(colorbar_title, rotation=270, labelpad=20, fontsize=16)

        # Add legend
        if legend:
            colors = np.append(plt.get_cmap(cmap_here)([0, 0.5, 1]), np.array([[1, 1, 1, 1]]), axis=0)
            labels = ['Downregulated', 'No change', 'Upregulated', 'No data']
            patches = [
                mpatches.Patch(facecolor=colors[i], label=labels[i], edgecolor='black')
                for i in range(len(labels))
            ]
            ax.legend(handles=patches, bbox_to_anchor=(1.05, 1),
                      loc=2, borderaxespad=0., title='Regulation',
                      fontsize=yticks_fontsize*2, title_fontsize=24)

        # Save figures
        if save_name:
            plt.savefig(output_dir / f'{save_name}.pdf', bbox_inches='tight')
            plt.savefig(output_dir / f'{save_name}.png', bbox_inches='tight', dpi=600)
            plt.savefig(output_dir / f'{save_name}.svg', bbox_inches='tight',
                       format='svg', transparent=True)
            logger.info(f"Saved heatmap to: {output_dir / save_name}.{{pdf,png,svg}}")

        if show:
            plt.show()
        else:
            plt.close()

        return fig

    # 1. Top 50 most changing pathways (first vs last timepoint)
    logger.info("Creating heatmap 1/5: Top 50 most changing pathways...")
    first_col = heatmap_data.columns[0]
    last_col = heatmap_data.columns[-1]
    change = (heatmap_data[last_col] - heatmap_data[first_col]).abs()
    top_changing = change.nlargest(50).index.tolist()

    df_top_changing = heatmap_data.loc[top_changing]
    fig1 = plot_heatmap_regulation(
        df_top_changing,
        top_changing,
        cmap_here='RdBu_r',
        save_name='top50_most_changing_pathways',
        colorbar_title='Sum of NES'
    )
    figures['top50_changing'] = fig1

    # 2. Top 50 most upregulated pathways (average NES > 0)
    logger.info("Creating heatmap 2/5: Top 50 most upregulated pathways...")
    avg_nes = heatmap_data.mean(axis=1)
    upregulated = avg_nes[avg_nes > 0].nlargest(50).index.tolist()

    df_upregulated = heatmap_data.loc[upregulated]
    fig2 = plot_heatmap_regulation(
        df_upregulated,
        upregulated,
        cmap_here='YlGn',
        save_name='top50_most_upregulated_pathways',
        colorbar_title='Sum of NES'
    )
    figures['top50_upregulated'] = fig2

    # 3. Top 50 most downregulated pathways (average NES < 0)
    logger.info("Creating heatmap 3/5: Top 50 most downregulated pathways...")
    downregulated = avg_nes[avg_nes < 0].nsmallest(50).index.tolist()

    df_downregulated = heatmap_data.loc[downregulated]
    fig3 = plot_heatmap_regulation(
        df_downregulated,
        downregulated,
        cmap_here='YlOrBr',
        save_name='top50_most_downregulated_pathways',
        colorbar_title='Sum of NES'
    )
    figures['top50_downregulated'] = fig3

    # 4. High-level pathways (29 pathways from Reactome highest level)
    logger.info("Creating heatmap 4/5: High-level pathways...")
    available_highest = [p for p in highest_pathways if p in heatmap_data.index]

    if available_highest:
        df_highest = heatmap_data.loc[available_highest]
        fig4 = plot_heatmap_regulation(
            df_highest,
            available_highest,
            cmap_here='RdBu_r',
            save_name='selected_pathways_highest_level',
            colorbar_title='Sum of NES'
        )
        figures['selected_highest_level'] = fig4
        logger.info(f"Found {len(available_highest)}/{len(highest_pathways)} high-level pathways in data")
    else:
        logger.warning("No high-level pathways found in the data")

    # 5. Literature pathways (33 pathways from literature review)
    logger.info("Creating heatmap 5/5: Literature pathways...")
    available_literature = [p for p in pathways_literature if p in heatmap_data.index]

    if available_literature:
        df_literature = heatmap_data.loc[available_literature]
        fig5 = plot_heatmap_regulation(
            df_literature,
            available_literature,
            cmap_here='RdBu_r',
            save_name='selected_pathways_literature',
            colorbar_title='Sum of NES'
        )
        figures['selected_literature'] = fig5
        logger.info(f"Found {len(available_literature)}/{len(pathways_literature)} literature pathways in data")
    else:
        logger.warning("No literature pathways found in the data")

    logger.info(f"Pathway heatmap generation complete. Created {len(figures)} heatmaps.")

    return heatmap_data, figures

