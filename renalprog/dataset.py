"""
Dataset loading and train/test splitting functionality for renalprog.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from typing import Tuple, Optional
import logging
import requests
import gzip
from renalprog.config import PreprocessingConfig
from renalprog.utils import set_seed

logger = logging.getLogger(__name__)

def download_data(
        destination: Path = "data/raw",
        remove_gz: bool = True,
        timeout: int = 300) -> Tuple[Path,Path,Path]:
    """
    Download the KIRC datasets to the specified destination.

    Args:
        destination: Path to directory where dataset should be saved
        remove_gz: Whether to remove .gz files after extraction
        timeout: Request timeout in seconds
    """
    # Ensure save directory exists
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)

    datasets = [
        ("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
         "EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",
         "Gene expression RNAseq - Batch effects normalized mRNA data"),
        ("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
         "TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
         "TCGA phenotype data"),
        ("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
         "Survival_SupplementalTable_S1_20171025_xena_sp",
         "Clinical survival data")
    ]

    for url, description in datasets:
        try:
            filename = url.split('/')[-1]
            file_path = destination / filename
            is_gzipped = filename.endswith('.gz')

            logger.info(f"Downloading dataset: {description}...")
            response = requests.get(url, timeout=timeout, stream=True)
            response.raise_for_status()

            total_size = int(response.headers.get('content-length', 0))
            with open(file_path, 'wb') as file:
                downloaded = 0
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        file.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            print(f"\rDownloading {filename}: {percent:.1f}%", end='', flush=True)

            print(f"\nFile downloaded successfully: {file_path}")

            if is_gzipped:
                extracted_path = destination / filename.replace('.gz', '')
                logger.info("Extracting compressed file...")
                with gzip.open(file_path, 'rb') as f_in, open(extracted_path, 'wb') as f_out:
                    f_out.write(f_in.read())
                logger.info(f"Successfully extracted file to: {extracted_path}")

                if remove_gz:
                    file_path.unlink()
                    logger.info("Removed compressed .gz file")

        except requests.RequestException as e:
            logger.error(f"Failed to download {description}: {e}")
            raise
        except IOError as e:
            logger.error(f"File I/O error for {description}: {e}")
            raise

    # Return paths to downloaded files
    rnaseq_path = destination / "EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
    clinical_path = destination / "Survival_SupplementalTable_S1_20171025_xena_sp"
    phenotype_path = destination / "TCGA_phenotype_denseDataOnlyDownload.tsv"

    return rnaseq_path, clinical_path, phenotype_path


def download_preprocessed_from_zenodo(
        rnaseq_url: str,
        clinical_url: str,
        output_dir: Path,
        rnaseq_filename: str = "preprocessed_rnaseq.csv",
        clinical_filename: str = "clinical_data.csv",
        timeout: int = 300) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Download preprocessed RNAseq and clinical data from Zenodo.

    This function downloads already preprocessed datasets from Zenodo repositories,
    which bypasses the need for raw data processing. The downloaded data includes
    normalized RNAseq expression matrices and curated clinical metadata.

    Args:
        rnaseq_url: Direct download URL for the preprocessed RNAseq CSV file from Zenodo.
            Example: 'https://zenodo.org/records/17987300/files/Static_KIRC.csv?download=1'
        clinical_url: Direct download URL for the clinical metadata CSV file from Zenodo.
            Example: 'https://zenodo.org/records/17987300/files/nodes_metadata.csv?download=1'
        output_dir: Directory where downloaded files will be saved.
        rnaseq_filename: Filename to save the RNAseq data as (default: "preprocessed_rnaseq.csv").
        clinical_filename: Filename to save the clinical data as (default: "clinical_data.csv").
        timeout: Request timeout in seconds (default: 300).

    Returns:
        Tuple of (rnaseq_data, clinical_data) as pandas DataFrames:
            - rnaseq_data: Gene expression matrix (genes × samples)
            - clinical_data: Clinical metadata (samples × features)

    Raises:
        requests.RequestException: If download fails
        IOError: If file I/O operations fail
        ValueError: If downloaded data cannot be parsed

    Notes:
        - Downloaded files are automatically saved to the specified output directory
        - Progress is logged during download
        - Data is validated after download to ensure proper format

    Examples:
        >>> # Download preprocessed KIRC data
        >>> rnaseq, clinical = download_preprocessed_from_zenodo(
        ...     rnaseq_url='https://zenodo.org/records/17987300/files/Static_KIRC.csv?download=1',
        ...     clinical_url='https://zenodo.org/records/17987300/files/nodes_metadata.csv?download=1',
        ...     output_dir=Path('data/interim/preprocessed_KIRC_data')
        ... )
    """
    # Ensure output directory exists
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("="*80)
    logger.info("Downloading preprocessed data from Zenodo")
    logger.info("="*80)

    # =========================================================================
    # Download RNAseq data
    # =========================================================================
    rnaseq_path = output_dir / rnaseq_filename
    logger.info("\n[DOWNLOADING RNASEQ DATA]")
    logger.info(f"  Source URL: {rnaseq_url}")
    logger.info(f"  Destination: {rnaseq_path}")

    try:
        response = requests.get(rnaseq_url, timeout=timeout, stream=True)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))
        with open(rnaseq_path, 'wb') as file:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    file.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print(f"\rDownloading RNAseq: {percent:.1f}%", end='', flush=True)

        print()  # New line after progress
        logger.info("  ✓ Successfully downloaded RNAseq data")

        # Load and validate the data
        rnaseq_data = pd.read_csv(rnaseq_path, index_col=0)
        logger.info(f"  Shape: {rnaseq_data.shape[0]:,} genes × {rnaseq_data.shape[1]:,} samples")

    except requests.RequestException as e:
        logger.error(f"Failed to download RNAseq data: {e}")
        raise
    except Exception as e:
        logger.error(f"Failed to load RNAseq data: {e}")
        raise

    # =========================================================================
    # Download clinical data
    # =========================================================================
    clinical_path = output_dir / clinical_filename
    logger.info("\n[DOWNLOADING CLINICAL DATA]")
    logger.info(f"  Source URL: {clinical_url}")
    logger.info(f"  Destination: {clinical_path}")

    try:
        response = requests.get(clinical_url, timeout=timeout, stream=True)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))
        with open(clinical_path, 'wb') as file:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    file.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print(f"\rDownloading clinical: {percent:.1f}%", end='', flush=True)

        print()  # New line after progress
        logger.info("  ✓ Successfully downloaded clinical data")

        # Load and validate the data
        clinical_data = pd.read_csv(clinical_path, index_col=0)
        logger.info(f"  Shape: {clinical_data.shape[0]:,} samples × {clinical_data.shape[1]:,} features")

    except requests.RequestException as e:
        logger.error(f"Failed to download clinical data: {e}")
        raise
    except Exception as e:
        logger.error(f"Failed to load clinical data: {e}")
        raise

    # =========================================================================
    # Validate data consistency
    # =========================================================================
    logger.info("\n[DATA VALIDATION]")

    # Check for sample overlap
    rnaseq_samples = set(rnaseq_data.columns)
    clinical_samples = set(clinical_data.index)
    common_samples = rnaseq_samples & clinical_samples

    logger.info(f"  RNAseq samples: {len(rnaseq_samples):,}")
    logger.info(f"  Clinical samples: {len(clinical_samples):,}")
    logger.info(f"  Common samples: {len(common_samples):,}")

    if len(common_samples) == 0:
        logger.warning("  ⚠ WARNING: No common samples found between RNAseq and clinical data!")
    elif len(common_samples) < len(rnaseq_samples):
        logger.warning(f"  ⚠ WARNING: Only {len(common_samples)}/{len(rnaseq_samples)} RNAseq samples have clinical data")

    logger.info("\n" + "="*80)
    logger.info("Preprocessed data download completed successfully")
    logger.info("="*80 + "\n")

    return rnaseq_data, clinical_data


def process_downloaded_data(
        rnaseq_path: Path = "data/raw/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena",
        clinical_path: Path = "data/raw/Survival_SupplementalTable_S1_20171025_xena_sp",
        phenotype_path: Path = "data/raw/TCGA_phenotype_denseDataOnlyDownload.tsv",
        cancer_type: str = "KIRC",
        output_dir: Path = "data/raw",
        extract_controls: bool = True,
        early_late: bool = False) -> Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    """
    Process TCGA Pan-Cancer Atlas data for a specific cancer type.

    This function performs comprehensive data processing and quality control on TCGA
    datasets, including cancer type filtering, sample type selection, stage harmonization,
    and optional binary classification mapping.

    Processing Steps:
        1. Load raw RNA-seq, clinical, and phenotype data from TCGA Xena Hub
        2. Filter samples by cancer type (KIRC or BRCA)
        3. Select primary tumor samples only (exclude metastatic/recurrent)
           and controls (solid tissue normal samples)
        4. Remove ambiguous stage annotations (Stage X, discrepancies, missing)
        5. Harmonize substages (e.g., Stage IA/IB/IC → Stage I)
        6. Optionally map stages to binary early/late classification
        7. Ensure consistency across all three datasets
        8. Save processed data in CSV format

    Args:
        rnaseq_path: Path to RNA-seq expression matrix (genes × samples).
            Expected format: tab-delimited file with gene IDs as rows and
            sample IDs (TCGA barcodes) as columns.
        clinical_path: Path to clinical survival data.
            Expected format: tab-delimited file with survival information,
            stage annotations, and patient metadata.
        phenotype_path: Path to phenotype annotations.
            Expected format: tab-delimited file with sample type information
            (Primary Tumor, Metastatic, etc.).
        cancer_type: Cancer type abbreviation. Supported values:
            - "KIRC": Kidney Renal Clear Cell Carcinoma
            - "BRCA": Breast Invasive Carcinoma (filters to female patients only)
        output_dir: Directory where processed CSV files will be saved.
        extract_controls: If True, extract and save control samples, named as
            Solid Tissue Normal, into a separate "controls" subdirectory within
            output_dir.
        early_late: If True, map AJCC stages to binary classification:
            - "early": Stage I and Stage II
            - "late": Stage III and Stage IV
            If False, retain original stage granularity (I, II, III, IV).

    Returns:
        Tuple of Paths to processed files:
            - rnaseq_path: Processed RNA-seq expression matrix
            - clinical_path: Processed clinical annotations
            - phenotype_path: Processed phenotype data
            - control_path: Processed control sample data (if applicable)

    Raises:
        FileNotFoundError: If input files do not exist
        KeyError: If required columns are missing from input data
        ValueError: If cancer_type is not supported

    Notes:
        - For BRCA, only female patients with ductal or lobular carcinomas are retained
        - AJCC pathologic tumor stage is used as the primary staging system
        - Substages (A, B, C) are collapsed to main stages for statistical power
        - All three output datasets maintain consistent sample identifiers

    Examples:
        >>> # Process KIRC data with 4-stage classification
        >>> paths = process_downloaded_data(
        ...     rnaseq_path="data/raw/expression.xena",
        ...     clinical_path="data/raw/survival.tsv",
        ...     phenotype_path="data/raw/phenotype.tsv",
        ...     cancer_type="KIRC",
        ...     output_dir="data/processed",
        ...     early_late=False
        ... )

        >>> # Process KIRC data with binary early/late classification
        >>> paths = process_downloaded_data(
        ...     cancer_type="KIRC",
        ...     output_dir="data/processed",
        ...     early_late=True
        ... )
    """
    # =========================================================================
    # STEP 1: Load raw data from TCGA sources
    # =========================================================================
    logger.info("="*80)
    logger.info(f"Starting data processing pipeline for {cancer_type}")
    logger.info("="*80)

    rnaseq = pd.read_table(rnaseq_path, index_col=0)
    clinical = pd.read_csv(clinical_path, sep="\t", index_col=0)
    pheno = pd.read_table(phenotype_path, sep="\t", index_col=0)

    # Remove redundant _PATIENT column from clinical data
    clinical.drop("_PATIENT", axis=1, inplace=True, errors='ignore')

    logger.info("\n[INITIAL DATA SHAPES]")
    logger.info(f"  RNA-seq:   {rnaseq.shape[0]:>6,} genes × {rnaseq.shape[1]:>5,} samples")
    logger.info(f"  Clinical:  {clinical.shape[0]:>6,} patients × {clinical.shape[1]:>3,} features")
    logger.info(f"  Phenotype: {pheno.shape[0]:>6,} samples × {pheno.shape[1]:>3,} features")

    # =========================================================================
    # STEP 2: Filter by cancer type
    # =========================================================================
    logger.info(f"\n[FILTERING BY CANCER TYPE: {cancer_type}]")

    if cancer_type == 'BRCA':
        # For breast cancer, only include female patients
        clinical = clinical[
            (clinical["gender"] == "FEMALE") &
            (clinical["cancer type abbreviation"] == cancer_type)
        ]
        logger.info(f"  Filter: Female patients with {cancer_type}")
    elif cancer_type == 'KIRC':
        clinical = clinical[clinical["cancer type abbreviation"] == cancer_type]
        logger.info(f"  Filter: Patients with {cancer_type}")
    else:
        logger.warning(f"  Cancer type '{cancer_type}' may not be fully supported")
        clinical = clinical[clinical["cancer type abbreviation"] == cancer_type]

    # Synchronize datasets to common samples
    pheno = pheno[pheno.index.isin(clinical.index)]
    rnaseq = rnaseq.loc[:, rnaseq.columns.isin(clinical.index)]
    pheno = pheno[pheno.index.isin(rnaseq.columns)]
    clinical = clinical[clinical.index.isin(rnaseq.columns)]

    logger.info("\n  Post-filter shapes:")
    logger.info(f"    RNA-seq:   {rnaseq.shape[0]:>6,} genes × {rnaseq.shape[1]:>5,} samples")
    logger.info(f"    Clinical:  {clinical.shape[0]:>6,} patients")
    logger.info(f"    Phenotype: {pheno.shape[0]:>6,} samples")

    # =========================================================================
    # STEP 3: Select primary tumor samples only
    # =========================================================================
    logger.info("\n[FILTERING BY SAMPLE TYPE: Primary Tumor]")

    # Log sample type distribution before filtering
    sample_type_counts = pheno["sample_type"].value_counts()
    logger.info("  Sample type distribution:")
    for sample_type, count in sample_type_counts.items():
        logger.info(f"    {sample_type:<30} {count:>5,} samples")

    # Extract controls samples if specified
    if extract_controls:
        logger.info("\n[EXTRACT CONTROLS]")
        pheno_controls = pheno[pheno["sample_type"] == "Solid Tissue Normal"]
        clinical_controls = clinical[clinical.index.isin(pheno.index)]
        rnaseq_controls = rnaseq.loc[:, rnaseq.columns.isin(pheno_controls.index)]
        logger.info(f"  Found {pheno_controls.shape[0]:>6,} control samples")
        # Save control data
        control_output_dir = Path(output_dir) / "controls"
        control_output_dir.mkdir(parents=True, exist_ok=True)

        rnaseq_controls.to_csv(control_output_dir / f"{cancer_type}_control_rnaseq.csv")
        clinical_controls.to_csv(control_output_dir / f"{cancer_type}_control_clinical.csv")
        pheno_controls.to_csv(control_output_dir / f"{cancer_type}_control_phenotype.csv")

        logger.info(f"  Saved control data to {control_output_dir}")

    pheno = pheno[pheno["sample_type"] == "Primary Tumor"]
    clinical = clinical[clinical.index.isin(pheno.index)]
    rnaseq = rnaseq.loc[:, rnaseq.columns.isin(pheno.index)]

    logger.info(f"\n  Retained {len(pheno):,} primary tumor samples")
    logger.info(f"    RNA-seq:   {rnaseq.shape[0]:>6,} genes × {rnaseq.shape[1]:>5,} samples")
    logger.info(f"    Clinical:  {clinical.shape[0]:>6,} patients")
    logger.info(f"    Phenotype: {pheno.shape[0]:>6,} samples")

    # =========================================================================
    # STEP 4: Remove ambiguous stage annotations
    # =========================================================================
    logger.info("\n[FILTERING BY STAGE QUALITY]")

    stages_remove = ['Stage X', '[Discrepancy]', np.nan]
    initial_stage_counts = clinical["ajcc_pathologic_tumor_stage"].value_counts(dropna=False)

    logger.info("  Initial stage distribution:")
    _log_stage_table(initial_stage_counts)

    # Filter out problematic stage annotations
    clinical_redux = clinical[~clinical["ajcc_pathologic_tumor_stage"].isin(stages_remove)]
    removed_count = len(clinical) - len(clinical_redux)
    logger.info(f"\n  Removed {removed_count} samples with ambiguous staging")

    # =========================================================================
    # STEP 5: Harmonize substages (collapse A/B/C to main stage)
    # =========================================================================
    logger.info("\n[HARMONIZING SUBSTAGES]")

    stages_available = clinical_redux["ajcc_pathologic_tumor_stage"].unique()
    has_substages = any(
        str(stage).endswith(("A", "B", "C"))
        for stage in stages_available
        if pd.notna(stage)
    )

    if has_substages:
        logger.info("  Substages detected (e.g., Stage IA, Stage IIB)")
        logger.info("  Collapsing substages to main stages for statistical power")

        # Collapse substages by removing trailing A/B/C letters
        stages_clump = [
            str(stage)[:-1] if str(stage).endswith(("A", "B", "C")) else stage
            for stage in clinical_redux["ajcc_pathologic_tumor_stage"]
        ]
        clinical_redux["ajcc_pathologic_tumor_stage"] = stages_clump

        post_clump_counts = clinical_redux["ajcc_pathologic_tumor_stage"].value_counts().sort_index()
        logger.info("\n  Stage distribution after harmonization:")
        _log_stage_table(post_clump_counts)
    else:
        logger.info("  No substages detected; stage labels are already harmonized")

    # Synchronize datasets after stage filtering
    pheno_redux = pheno[pheno.index.isin(clinical_redux.index)]
    rnaseq_redux = rnaseq.loc[:, rnaseq.columns.isin(clinical_redux.index)]

    # =========================================================================
    # STEP 6: Additional cancer-specific filtering (BRCA only)
    # =========================================================================
    if cancer_type == 'BRCA':
        logger.info("\n[BRCA-SPECIFIC FILTERING: Ductal and Lobular Carcinomas]")
        ductal_lobular = ['Infiltrating Ductal Carcinoma', 'Infiltrating Lobular Carcinoma']
        patients_ductal_lobular = clinical_redux.index[
            clinical_redux["histological_type"].isin(ductal_lobular)
        ]

        clinical_redux = clinical_redux[clinical_redux.index.isin(patients_ductal_lobular)]
        pheno_redux = pheno_redux[pheno_redux.index.isin(patients_ductal_lobular)]
        rnaseq_redux = rnaseq_redux.loc[:, rnaseq_redux.columns.isin(patients_ductal_lobular)]

        logger.info(f"  Retained {len(patients_ductal_lobular):,} ductal/lobular samples")

    # =========================================================================
    # STEP 7: Optional binary classification mapping (early vs. late)
    # =========================================================================
    if early_late:
        logger.info("\n[MAPPING TO BINARY CLASSIFICATION: Early vs. Late]")
        logger.info("  Mapping scheme:")
        logger.info("    Early: Stage I, Stage II")
        logger.info("    Late:  Stage III, Stage IV")

        mapped_stages = map_stages_to_early_late(clinical_redux["ajcc_pathologic_tumor_stage"])
        valid_mask = mapped_stages.notna()

        clinical_redux = clinical_redux.loc[valid_mask, :]
        clinical_redux["ajcc_pathologic_tumor_stage"] = mapped_stages[valid_mask]
        pheno_redux = pheno_redux[pheno_redux.index.isin(clinical_redux.index)]
        rnaseq_redux = rnaseq_redux.loc[:, rnaseq_redux.columns.isin(clinical_redux.index)]

        final_stage_counts = clinical_redux["ajcc_pathologic_tumor_stage"].value_counts().sort_index()
        logger.info("\n  Final binary classification distribution:")
        _log_stage_table(final_stage_counts)
    else:
        final_stage_counts = clinical_redux["ajcc_pathologic_tumor_stage"].value_counts().sort_index()
        logger.info("\n  Final stage distribution (4-class):")
        _log_stage_table(final_stage_counts)

    # =========================================================================
    # STEP 8: Final statistics and data saving
    # =========================================================================
    logger.info("\n[FINAL PROCESSED DATA SHAPES]")
    logger.info(f"  RNA-seq:   {rnaseq_redux.shape[0]:>6,} genes × {rnaseq_redux.shape[1]:>5,} samples")
    logger.info(f"  Clinical:  {clinical_redux.shape[0]:>6,} patients × {clinical_redux.shape[1]:>3,} features")
    logger.info(f"  Phenotype: {pheno_redux.shape[0]:>6,} samples × {pheno_redux.shape[1]:>3,} features")

    # Calculate and log filtering statistics
    initial_samples = rnaseq.shape[1]
    final_samples = rnaseq_redux.shape[1]
    retention_rate = (final_samples / initial_samples) * 100 if initial_samples > 0 else 0

    logger.info("\n[FILTERING SUMMARY]")
    logger.info(f"  Initial samples:  {initial_samples:>5,}")
    logger.info(f"  Final samples:    {final_samples:>5,}")
    logger.info(f"  Retention rate:   {retention_rate:>5.1f}%")
    logger.info(f"  Samples removed:  {initial_samples - final_samples:>5,}")
    if extract_controls:
        logger.info(f"  Control samples found:  {rnaseq_controls.shape[1]:>5,}")
    # Save processed data to disk
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    path_rnaseq = output_dir / f"{cancer_type}_rnaseq.csv"
    path_clinical = output_dir / f"{cancer_type}_clinical.csv"
    path_phenotype = output_dir / f"{cancer_type}_phenotype.csv"

    logger.info("\n[SAVING PROCESSED DATA]")
    logger.info(f"  Output directory: {output_dir}")

    rnaseq_redux.to_csv(path_rnaseq)
    logger.info(f"  ✓ Saved RNA-seq:   {path_rnaseq.name}")

    clinical_redux.to_csv(path_clinical)
    logger.info(f"  ✓ Saved clinical:  {path_clinical.name}")

    pheno_redux.to_csv(path_phenotype)
    logger.info(f"  ✓ Saved phenotype: {path_phenotype.name}")

    logger.info("\n" + "="*80)
    logger.info(f"Data processing completed successfully for {cancer_type}")
    logger.info("="*80 + "\n")

    return rnaseq_redux, clinical_redux, pheno_redux


def _log_stage_table(stage_counts: pd.Series) -> None:
    """
    Log a beautified table of stage value counts.

    Args:
        stage_counts: Series with stage labels as index and counts as values
    """
    if len(stage_counts) == 0:
        logger.info("    No stage data available")
        return

    # Calculate proportions
    total = stage_counts.sum()

    # Create formatted table
    logger.info("    " + "-" * 50)
    logger.info(f"    {'Stage':<15} {'Count':>10} {'Proportion':>12} {'Bar':>10}")
    logger.info("    " + "-" * 50)

    for stage, count in stage_counts.items():
        proportion = (count / total) * 100
        # Create a simple bar chart with asterisks
        bar_length = int(proportion / 5)  # Scale to fit
        bar = "█" * bar_length

        stage_str = str(stage) if pd.notna(stage) else "Unknown"
        logger.info(f"    {stage_str:<15} {count:>10,} {proportion:>11.1f}% {bar}")

    logger.info("    " + "-" * 50)
    logger.info(f"    {'TOTAL':<15} {total:>10,} {100.0:>11.1f}%")
    logger.info("    " + "-" * 50)


def load_rnaseq_data(path: Path) -> pd.DataFrame:
    """
    Load RNA-seq data from CSV file.
    
    Args:
        path: Path to RNA-seq CSV file (genes as rows, samples as columns)
        
    Returns:
        DataFrame with RNA-seq data
    """
    logger.info(f"Loading RNA-seq data from {path}")
    data = pd.read_csv(path, index_col=0)
    logger.info(f"Loaded RNA-seq data with shape: {data.shape}")
    return data


def load_clinical_data(path: Path, stage_column: str = "ajcc_pathologic_tumor_stage", early_late = True) -> pd.Series:
    """
    Load clinical metadata.
    
    Args:
        path: Path to clinical data CSV file
        stage_column: Name of column containing stage information
        
    Returns:
        Series with clinical stages indexed by sample ID
    """
    logger.info(f"Loading clinical data from {path}")
    data = pd.read_csv(path, index_col=0)
    
    if stage_column not in data.columns:
        raise ValueError(f"Stage column '{stage_column}' not found in clinical data")
    
    stages = data[stage_column]
    logger.info(f"Loaded clinical data with {len(stages)} samples")

    if early_late:
        # Check if the stages are already in early/late format
        unique_stages = stages.dropna().unique()
        if set(unique_stages).issubset({"early", "late"}):
            logger.info("Stages are already in early/late format; no mapping needed")
        else:
            stages = map_stages_to_early_late(stages)
            logger.info("Mapped stages to binary early/late classification")
    logger.info(f"Stage distribution:\n{stages.value_counts()}")
    
    return stages


def map_stages_to_early_late(stages: pd.Series) -> pd.Series:
    """
    Map detailed stages (I, II, III, IV) to binary early/late classification.
    
    Args:
        stages: Series with stage labels (e.g., "Stage I", "Stage II", etc.)
        
    Returns:
        Series with mapped stages ("early" or "late")
    """
    stage_mapping = PreprocessingConfig.STAGE_MAPPING
    mapped_stages = stages.map(stage_mapping)
    
    # Check for unmapped values
    unmapped = mapped_stages.isna() & stages.notna()
    if unmapped.any():
        logger.warning(f"Found {unmapped.sum()} unmapped stage values:")
        logger.warning(stages[unmapped].unique())
    
    return mapped_stages


def create_train_test_split(
    rnaseq_path: Path,
    clinical_path: Path,
    stage_column: str = "ajcc_pathologic_tumor_stage",
    test_size: float = 0.2,
    seed: int = 2023,
    use_onehot: bool = True,
    output_dir: Optional[Path] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray, pd.DataFrame, pd.Series]:
    """
    Create stratified train/test split of KIRC data.
    
    Args:
        rnaseq_path: Path to RNA-seq CSV file
        clinical_path: Path to clinical CSV file
        stage_column: Name of column containing stage information (default: "ajcc_pathologic_tumor_stage")
        test_size: Fraction of data to use for testing (default: 0.2)
        seed: Random seed for reproducibility (default: 2023)
        use_onehot: Whether to one-hot encode the labels (default: True)
        output_dir: Optional directory to save split data
        
    Returns:
        Tuple of (X_train, X_test, y_train, y_test, full_rnaseq, full_clinical)
    """
    set_seed(seed)
    
    # Load data
    rnaseq = load_rnaseq_data(rnaseq_path)
    clinical = load_clinical_data(clinical_path,stage_column)
    
    # Ensure samples match between RNA-seq and clinical data
    common_samples = rnaseq.columns.intersection(clinical.index)
    if len(common_samples) < len(rnaseq.columns):
        logger.warning(
            f"Only {len(common_samples)} of {len(rnaseq.columns)} samples "
            f"have clinical data. Filtering to common samples."
        )
    
    rnaseq = rnaseq[common_samples]
    clinical = clinical[common_samples]
    
    # Transpose to have samples as rows
    rnaseq_t = rnaseq.T
    
    # Prepare labels for stratification
    if use_onehot:
        # Get unique stages in sorted order for consistent encoding
        categories = sorted(clinical.unique())
        ohe = OneHotEncoder(
            categories=[categories],
            handle_unknown='ignore',
            sparse_output=False,
            dtype=np.int8
        )
        y = ohe.fit_transform(clinical.values.reshape(-1, 1))
        logger.info(f"One-hot encoded labels with {y.shape[1]} classes")
    else:
        y = clinical.values
    
    # Perform stratified split
    X_train, X_test, y_train, y_test = train_test_split(
        rnaseq_t,
        y,
        test_size=test_size,
        stratify=y if use_onehot else clinical,
        random_state=seed
    )
    
    logger.info(f"Train set: {X_train.shape[0]} samples")
    logger.info(f"Test set: {X_test.shape[0]} samples")
    
    # Save split data
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save split data
        X_train.to_csv(output_dir / "X_train.csv")
        X_test.to_csv(output_dir / "X_test.csv")
        
        # Save labels WITH patient IDs as indices (matching X_train and X_test)
        if use_onehot:
            # Create DataFrame with patient IDs as index
            y_train_df = pd.DataFrame(y_train, index=X_train.index)
            y_test_df = pd.DataFrame(y_test, index=X_test.index)
            y_train_df.to_csv(output_dir / "y_train.csv")
            y_test_df.to_csv(output_dir / "y_test.csv")
        else:
            # Create Series with patient IDs as index
            y_train_series = pd.Series(y_train, index=X_train.index)
            y_test_series = pd.Series(y_test, index=X_test.index)
            y_train_series.to_csv(output_dir / "y_train.csv")
            y_test_series.to_csv(output_dir / "y_test.csv")

        # Save full data for reference
        rnaseq.to_csv(output_dir / "data.csv")
        clinical.to_csv(output_dir / "metadata.csv")
        
        # Save split statistics
        _save_split_statistics(clinical, y_train, y_test, output_dir, use_onehot, categories if use_onehot else None)
        
        logger.info(f"Saved train/test split to {output_dir}")
    
    return X_train, X_test, y_train, y_test, rnaseq, clinical


def _save_split_statistics(
    clinical: pd.Series,
    y_train: np.ndarray,
    y_test: np.ndarray,
    output_dir: Path,
    use_onehot: bool,
    categories: Optional[list] = None
):
    """Save statistics about the train/test split."""
    
    # Convert one-hot back to labels if needed
    if use_onehot:
        y_train_labels = pd.Series([categories[i] for i in y_train.argmax(axis=1)])
        y_test_labels = pd.Series([categories[i] for i in y_test.argmax(axis=1)])
    else:
        y_train_labels = pd.Series(y_train)
        y_test_labels = pd.Series(y_test)
    
    # Create statistics dataframe
    stats = pd.DataFrame({
        "original_count": clinical.value_counts(),
        "train_count": y_train_labels.value_counts(),
        "test_count": y_test_labels.value_counts()
    })
    
    # Add proportions
    stats["original_proportion"] = stats["original_count"] / stats["original_count"].sum()
    stats["train_proportion"] = stats["train_count"] / stats["train_count"].sum()
    stats["test_proportion"] = stats["test_count"] / stats["test_count"].sum()
    
    stats.to_csv(output_dir / "split_statistics.csv")
    logger.info(f"Split statistics:\n{stats}")


def load_train_test_split(split_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load previously saved train/test split.
    
    Args:
        split_dir: Directory containing saved split files
        
    Returns:
        Tuple of (X_train, X_test, y_train, y_test)
    """
    split_dir = Path(split_dir)
    
    X_train = pd.read_csv(split_dir / "X_train.csv", index_col=0)
    X_test = pd.read_csv(split_dir / "X_test.csv", index_col=0)
    y_train = pd.read_csv(split_dir / "y_train.csv", index_col=0)
    y_test = pd.read_csv(split_dir / "y_test.csv", index_col=0)
    
    logger.info(f"Loaded train/test split from {split_dir}")
    logger.info(f"Train: {X_train.shape}, Test: {X_test.shape}")
    
    return X_train, X_test, y_train, y_test
