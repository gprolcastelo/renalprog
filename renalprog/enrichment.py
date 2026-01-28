r"""
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

IMPORTANT: This module uses PyDESeq2 for proper differential expression analysis.
DO NOT use simple log-fold change calculations on preprocessed data, as this
produces biologically invalid results.

Author: from Guillermo Prol-Castelo
Date: December 19, 2025
License: Apache 2.0
"""

import os
import logging
from pathlib import Path
from typing import Dict, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    import matplotlib.figure

import pandas as pd
import numpy as np
# Import PyDESeq2 for proper differential expression analysis
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats



logger = logging.getLogger(__name__)


def build_gsea_command(save_path_rnk_fun:str, filename_rnk_fun:str,
                       pathway_file:str='data/external/ReactomePathways.gmt',
                       mode:str='Max_probe',norm:str='meandiv', nperm:int=1000, rnd_seed:str='timestamp',
                       scoring_scheme:str='weighted', set_max:int=500, set_min:int=15):
    """
    Build a command string for running GSEA (Gene Set Enrichment Analysis) with pre-ranked gene lists.

    Parameters:
    - save_path_rnk_fun (str): The directory where the pre-ranked gene list file is stored.
    - filename_rnk_fun (str): The name of the pre-ranked gene list file.
    - pathway_file (str, optional): Path to the gene set file in GMT format. Default is 'data/external/ReactomePathways.gmt'.
    - mode (str, optional): GSEA mode. Default is 'Max_probe'.
    - norm (str, optional): Normalization method. Default is 'meandiv'.
    - nperm (int, optional): Number of permutations. Default is 1000.
    - rnd_seed (str, optional): Random seed for reproducibility. Default is 'timestamp'.
    - scoring_scheme (str, optional): Scoring scheme for GSEA. Default is 'weighted'.
    - set_max (int, optional): Maximum size of gene sets to consider. Default is 500.
    - set_min (int, optional): Minimum size of gene sets to consider. Default is 15.

    Returns:
    str: The GSEA command string.
    """

    # Start of the GSEA command
    start = 'gsea GSEAPreranked '

    # File path of the pre-ranked gene list
    rnk = '-rnk ' + os.path.join(save_path_rnk_fun, filename_rnk_fun) + ' '

    # Path to the gene set file in GMT format
    gmx = f'-gmx {pathway_file} '

    # GSEA options
    options = (
            '-collapse No_Collapse '
            + f'-mode {mode} '
            + f'-norm {norm} -nperm {nperm} '
            + f'-rnd_seed {rnd_seed} '
            + f'-scoring_scheme {scoring_scheme} '
    )

    # Report label option
    label = "-rpt_label " + filename_rnk_fun[:-4] + " "

    # Additional GSEA options
    options2 = (
            '-create_svgs false '
            + '-include_only_symbols true '
            + '-make_sets false '
            + '-plot_top_x 20 '
            + f'-set_max {set_max} '
            + f'-set_min {set_min} '
            + '-zip_report false '
    )

    # Output directory option
    out = '-out ' + os.path.join(save_path_rnk_fun, 'reports/')

    # Concatenate all components to form the complete GSEA command
    gsea_command = start + rnk + gmx + options + label + options2 + out

    return gsea_command


def get_rnk_single_patient(path_above, genes_here, samples_real_l, pat_i, index_pat=None, foldchange=True,
                          pathway_file='data/external/ReactomePathways.gmt',
                          mode='Max_probe', norm='meandiv', nperm=1000, rnd_seed='timestamp',
                          scoring_scheme='weighted', set_max=500, set_min=15):
    """
    Generate GSEA command for a single patient based on gene expression data.

    Parameters:
    - path_above (str): The parent directory for saving patient-specific data.
    - genes_here (list): List of gene names.
    - samples_real_l (pd.DataFrame): DataFrame containing gene expression data for single patient sample.
    - pat_i (str): Patient identifier.
    - index_pat (int, optional): Index of the patient. Default is None.
    - foldchange (bool, optional): If True, log2 fold change values are used; otherwise, raw gene expression values.
    - pathway_file (str, optional): Path to the GMT file containing gene sets. Default is 'data/external/ReactomePathways.gmt'.
    - mode (str, optional): GSEA collapse mode. Default is 'Max_probe'.
    - norm (str, optional): GSEA normalization mode. Default is 'meandiv'.
    - nperm (int, optional): Number of permutations. Default is 1000.
    - rnd_seed (str, optional): Random seed for GSEA. Default is 'timestamp'.
    - scoring_scheme (str, optional): GSEA scoring scheme. Default is 'weighted'.
    - set_max (int, optional): Maximum size of gene sets. Default is 500.
    - set_min (int, optional): Minimum size of gene sets. Default is 15.

    Returns:
    str: The GSEA command string for the given patient and gene expression data.
    """

    #############################################################
    # IF FOLDCHANGE IS TRUE
    #############################################################
    if foldchange:
        # Create directories if they don't exist
        os.makedirs(path_above, exist_ok=True)

        # DataFrame for patient i
        df_i = pd.DataFrame(samples_real_l['log2FoldChange'])
        df_i.index.name = '#'
        # print('path_above:', path_above)
        # print('pat_i:', pat_i)
        path_i = os.path.join(path_above, pat_i)

        os.makedirs(path_i, exist_ok=True)

        os.makedirs(os.path.join(path_i, 'reports'), exist_ok=True)

        # Save DataFrame as .rnk file
        if index_pat is not None:
            df_i.to_csv(os.path.join(path_i, pat_i + '_' + str(index_pat) + '.rnk'), sep='\t')
            # Get GSEA command
            gsea_command_k = build_gsea_command(
                save_path_rnk_fun=path_i,
                filename_rnk_fun=pat_i + '_' + str(index_pat) + '.rnk',
                pathway_file=pathway_file, mode=mode, norm=norm, nperm=nperm, rnd_seed=rnd_seed,
                scoring_scheme=scoring_scheme, set_max=set_max, set_min=set_min
            )
        else:
            df_i.to_csv(os.path.join(path_i, pat_i + '.rnk'), sep='\t')
            # Get GSEA command
            gsea_command_k = build_gsea_command(
                save_path_rnk_fun=path_i,
                filename_rnk_fun=pat_i + '.rnk',
                pathway_file=pathway_file, mode=mode, norm=norm, nperm=nperm, rnd_seed=rnd_seed,
                scoring_scheme=scoring_scheme, set_max=set_max, set_min=set_min
            )



    #############################################################
    # ELSE FOLDCHANGE FALSE
    #############################################################
    else:
        raise "This option has been deprecated. GSEA should use fold change values from DESeq2 analysis."

    return gsea_command_k


def fun_apply_deseq(rnaseq_in, clinical_in):
    '''
    Perform DESeq2 analysis on RNA-seq data.

    Parameters:
    - rnaseq_in (pd.DataFrame): RNA-seq count data.
    - clinical_in (pd.DataFrame): Clinical metadata, including the cancer stage.

    Returns:
    pd.DataFrame: DataFrame containing DESeq2 analysis results.
    '''

    # Re-scale to have RSEM
    # The input values are in log2(RSEM+1) format, so we need to convert them to RSEM, rounding to the nearest integer.
    rsem = pd.DataFrame(np.round(2 ** rnaseq_in - 1), dtype=int)
    # Get maximum value in rnaseq_in
    # Print dtype of rnaseq_in
    # Extract the condition (cancer stage)
    condition = pd.DataFrame(clinical_in['ajcc_pathologic_tumor_stage'])
    condition.columns = ['condition']

    # Construct DESEQDataSet Object
    dds_i = DeseqDataSet(
        counts=rsem,
        metadata=condition,
        design_factors='condition',
        refit_cooks=True,
        quiet=False
    )

    # Once a DeseqDataSet was initialized, run the deseq2() method to fit dispersions and LFCs
    dds_i.deseq2()

    # Statistical analysis
    stat_res = DeseqStats(dds_i, alpha=0.05, cooks_filter=True, independent_filter=True, quiet=False)

    # Summary of results
    stat_res.summary()

    return stat_res.results_df


def fun_single_patient_and_gsea(
        patient_here,
        patient_stage,
        rnaseq_data,
        path_above_in,
        clinical_controls,
        rnaseq_controls,
        gene_list,
        foldchange=True,
        pathway_file='data/external/ReactomePathways.gmt',
        mode='Max_probe',
        norm='meandiv',
        nperm=1000,
        rnd_seed='timestamp',
        scoring_scheme='weighted',
        set_max=500,
        set_min=15
):
    '''
    Perform analysis for a single patient and generate GSEA-related files.

    Parameters:
    - patient_here: Current patient ID.
    - patient_stage: DataFrame containing patient stage information.
    - rnaseq_data: DataFrame containing RNAseq data.
    - path_above_in: Path to the directory containing patient data.
    - clinical_controls: DataFrame containing clinical control data.
    - rnaseq_controls: DataFrame containing RNAseq control data.
    - gene_list: List of genes.
    - foldchange: Boolean flag indicating whether to use fold change data in the analysis. Default is True.
    - pathway_file (str, optional): Path to the GMT file containing gene sets. Default is 'data/external/ReactomePathways.gmt'.
    - mode (str, optional): GSEA collapse mode. Default is 'Max_probe'.
    - norm (str, optional): GSEA normalization mode. Default is 'meandiv'.
    - nperm (int, optional): Number of permutations. Default is 1000.
    - rnd_seed (str, optional): Random seed for GSEA. Default is 'timestamp'.
    - scoring_scheme (str, optional): GSEA scoring scheme. Default is 'weighted'.
    - set_max (int, optional): Maximum size of gene sets. Default is 500.
    - set_min (int, optional): Minimum size of gene sets. Default is 15.

    Returns:
    - bash_file_l: List of bash commands generated during the analysis.
    '''
    # Get stage of the current patient
    # stage_pat_i= str(patient_stage['stage'][patient_stage['name'] == patient_here].values[0])
    stage_pat_i = patient_stage
    # print('Stage ',stage_pat_i)

    # Prepare RNAseq and clinical data for analysis
    rna_kk = pd.DataFrame(rnaseq_data.loc[patient_here].values.reshape(1, -1),
                          index=[patient_here], columns=gene_list)
    rna_kk = pd.concat([rna_kk, rnaseq_controls], axis=0, ignore_index=False)

    cli_kk = pd.concat(
        [
            # pd.DataFrame([stage_pat_i], columns=['ajcc_pathologic_tumor_stage'], index=[patient_here]),
            pd.DataFrame(stage_pat_i, columns=['ajcc_pathologic_tumor_stage'], index=[patient_here]),
            clinical_controls,
        ],
        axis=0,
        ignore_index=False,
    )
    # print('cli_kk head():', cli_kk.head())
    # Use DESeq to get fold change for each interpolated patient
    fold_change_df_i = fun_apply_deseq(rnaseq_in=rna_kk, clinical_in=cli_kk)

    # Generate GSEA-related files and bash commands
    # print('path_above_in:', path_above_in)
    # print('stage_pat_i:', stage_pat_i)
    path_above = os.path.join(path_above_in, stage_pat_i)
    bash_file_l = get_rnk_single_patient(
        path_above=path_above,
        genes_here=gene_list,
        samples_real_l=fold_change_df_i,
        pat_i=patient_here,
        foldchange=foldchange,
        pathway_file=pathway_file,
        mode=mode,
        norm=norm,
        nperm=nperm,
        rnd_seed=rnd_seed,
        scoring_scheme=scoring_scheme,
        set_max=set_max,
        set_min=set_min
    )

    # print('Done with patient ' + stage_pat_i)
    return bash_file_l


def fun_synth_single_patient_and_gsea(
        patient_here,
        stage_trans_i,
        rnaseq_data,
        path_above_in,
        clinical_controls,
        rnaseq_controls,
        gene_list,
        index_pat=None,
        foldchange=True,
        pathway_file='data/external/ReactomePathways.gmt',
        mode='Max_probe',
        norm='meandiv',
        nperm=1000,
        rnd_seed='timestamp',
        scoring_scheme='weighted',
        set_max=500,
        set_min=15
):
    '''
    Perform analysis for a single patient and generate GSEA-related files.

    Parameters:
    - patient_here: str Current synth patient ID (like: pat_to_pat_interpolnum).
    - stage_trans_i: str containing synthetic patient stage transition information. 'I_to_II' or 'II_to_III'.
    - rnaseq_data: DataFrame containing RNAseq data only of synthetic patient to analyze.
    - path_above_in: Path to the directory containing patient data.
    - clinical_controls: DataFrame containing clinical control data.
    - rnaseq_controls: DataFrame containing RNAseq control data.
    - gene_list: List of genes.
    - index_pat: int, optional: Index of the synthetic patient. Default is None.
    - foldchange: Boolean flag indicating whether to use fold change data in the analysis. Default is True.
    - pathway_file (str, optional): Path to the GMT file containing gene sets. Default is 'data/external/ReactomePathways.gmt'.
    - mode (str, optional): GSEA collapse mode. Default is 'Max_probe'.
    - norm (str, optional): GSEA normalization mode. Default is 'meandiv'.
    - nperm (int, optional): Number of permutations. Default is 1000.
    - rnd_seed (str, optional): Random seed for GSEA. Default is 'timestamp'.
    - scoring_scheme (str, optional): GSEA scoring scheme. Default is 'weighted'.
    - set_max (int, optional): Maximum size of gene sets. Default is 500.
    - set_min (int, optional): Minimum size of gene sets. Default is 15.

    Returns:
    - bash_file_l: List of bash commands generated during the analysis.
    '''
    # Get stage of the current patient
    # stage_pat_i = str(patient_stage['stage'][patient_stage['name'] == patient_here].values[0])
    assert stage_trans_i in ['I', 'II', 'III', 'I_to_II', 'II_to_III', '1_to_2', '2_to_3', 'early_to_late',
                             'early_to_early'], 'Stage not recognized'

    # Prepare RNAseq and clinical data for analysis
    # rna_kk = pd.DataFrame(rnaseq_data.loc[patient_here].values.reshape(1, -1),
    #                       index=[patient_here], columns=gene_list)
    rna_kk = pd.DataFrame(rnaseq_data.values.reshape(1, -1))
    rna_kk.columns = gene_list
    rna_kk.index = [patient_here]
    rna_kk = pd.concat([rna_kk, rnaseq_controls], axis=0, ignore_index=False)
    assert ~rna_kk.isna().any().any(), "NAs detected in input data. Make sure all data is numeric"
    # Build clinical data with the patient's stage and the control data
    # The string stage_pat_i is transformed into a DataFrame with the same index as the
    # patient, and is then concatenated with the control data.
    cli_kk = pd.concat(
        [
            pd.DataFrame([stage_trans_i],
                         columns=['ajcc_pathologic_tumor_stage'],
                         index=[patient_here]),
            clinical_controls,
        ],
        axis=0,
        ignore_index=False,
    )

    # Use DESeq to get fold change for each interpolated patient
    print(f"Applying DESeq for synthetic patient {patient_here} with stage transition {stage_trans_i}")
    fold_change_df_i = fun_apply_deseq(rnaseq_in=rna_kk, clinical_in=cli_kk)

    # Generate GSEA-related files and bash commands
    bash_file_l = get_rnk_single_patient(
        path_above=os.path.join(path_above_in, stage_trans_i),
        genes_here=gene_list,
        samples_real_l=fold_change_df_i,
        pat_i=patient_here,
        index_pat=index_pat,
        foldchange=foldchange,
        pathway_file=pathway_file,
        mode=mode,
        norm=norm,
        nperm=nperm,
        rnd_seed=rnd_seed,
        scoring_scheme=scoring_scheme,
        set_max=set_max,
        set_min=set_min
    )

    # print('Done with patient ' + stage_trans_i)
    return bash_file_l


def generate_pathway_heatmap(
    enrichment_df: pd.DataFrame,
    output_dir: str,
    fdr_threshold: float = 0.05,
    colorbar: bool = True,
    legend: bool = False,
    yticks_fontsize: int = 12,
    show: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, "matplotlib.figure.Figure"]]:
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
        "SIGNALING BY VEGF",
    ]

    # Step 1: Ensure numeric types for NES and FDR q-val
    enrichment_df = enrichment_df.copy()
    enrichment_df["NES"] = pd.to_numeric(enrichment_df["NES"], errors="coerce")
    enrichment_df["FDR q-val"] = pd.to_numeric(
        enrichment_df["FDR q-val"], errors="coerce"
    )
    enrichment_df["Idx"] = pd.to_numeric(enrichment_df["Idx"], errors="coerce")

    # Log any rows that had non-numeric values
    invalid_nes = enrichment_df["NES"].isna().sum()
    invalid_fdr = enrichment_df["FDR q-val"].isna().sum()
    if invalid_nes > 0:
        logger.warning(
            f"Found {invalid_nes} rows with non-numeric NES values (converted to NaN)"
        )
    if invalid_fdr > 0:
        logger.warning(
            f"Found {invalid_fdr} rows with non-numeric FDR q-val values (converted to NaN)"
        )

    # Step 2: Filter by FDR threshold
    significant = enrichment_df[enrichment_df["FDR q-val"] < fdr_threshold].copy()

    logger.info(
        f"Found {significant.shape[0]} significant pathway enrichments (FDR < {fdr_threshold})"
    )

    if significant.empty:
        logger.warning("No significant pathways found. Cannot generate heatmap.")
        # Return empty results
        empty_df = pd.DataFrame()
        empty_dict = {}
        return empty_df, empty_dict

    # Step 3: Group by Timepoint (Idx) and Pathway (NAME), sum NES across all trajectories
    pathway_summary = significant.groupby(["Idx", "NAME"])["NES"].sum().reset_index()

    logger.info(
        f"Aggregated results for {pathway_summary['Idx'].nunique()} timepoints "
        f"and {pathway_summary['NAME'].nunique()} pathways"
    )

    # Step 4: Pivot to create matrix (pathways × timepoints)
    heatmap_data = pathway_summary.pivot(
        index="NAME", columns="Idx", values="NES"
    ).fillna(0)  # Fill missing with 0

    logger.info(
        f"Full heatmap dimensions: {heatmap_data.shape[0]} pathways × {heatmap_data.shape[1]} timepoints"
    )

    # Save full summary data
    summary_file = output_dir / "pathway_nes_summary.csv"
    heatmap_data.to_csv(summary_file)
    logger.info(f"Saved pathway NES summary to: {summary_file}")

    # Dictionary to store all figures
    figures = {}

    # Helper function to create and save a heatmap
    def plot_heatmap_regulation(
        df_plot,
        unique_pathways,
        cmap_here="viridis",
        save_name=None,
        colorbar_title="Sum of NES",
    ):
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
            logger.warning(
                f"Cannot center colormap at zero: range [{z_min:.3f}, {z_max:.3f}] does not cross zero"
            )

        fig, ax = plt.subplots(figsize=(30, 10))

        # Make heatmap
        cax = ax.imshow(df_plot.values, cmap=cmap_here, norm=norm, aspect="auto")

        # Set the y-ticks
        plt.yticks(tick_locations, unique_pathways, fontsize=yticks_fontsize)

        # Set the x-ticks at specific positions
        num_timepoints = df_plot.shape[1]
        ax.set_xticks([0, num_timepoints - 1])
        ax.set_xticklabels(
            ["early", "late"], fontsize=yticks_fontsize * 1.33, rotation=45
        )
        ax.set_xlabel("Pseudo-Time", fontsize=yticks_fontsize * 1.33)

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
        plt.grid(False, which="major")
        plt.grid(True, which="minor", color="black", linestyle="-", linewidth=1)

        # Remove ticks
        ax.tick_params(axis="both", which="minor", length=0)

        # Add colorbar
        if colorbar:
            cbar = plt.colorbar(cax, shrink=0.7)
            cbar.set_label(colorbar_title, rotation=270, labelpad=20, fontsize=16)

        # Add legend
        if legend:
            colors = np.append(
                plt.get_cmap(cmap_here)([0, 0.5, 1]), np.array([[1, 1, 1, 1]]), axis=0
            )
            labels = ["Downregulated", "No change", "Upregulated", "No data"]
            patches = [
                mpatches.Patch(facecolor=colors[i], label=labels[i], edgecolor="black")
                for i in range(len(labels))
            ]
            ax.legend(
                handles=patches,
                bbox_to_anchor=(1.05, 1),
                loc=2,
                borderaxespad=0.0,
                title="Regulation",
                fontsize=yticks_fontsize * 2,
                title_fontsize=24,
            )

        # Save figures
        if save_name:
            plt.savefig(output_dir / f"{save_name}.pdf", bbox_inches="tight")
            plt.savefig(output_dir / f"{save_name}.png", bbox_inches="tight", dpi=600)
            plt.savefig(
                output_dir / f"{save_name}.svg",
                bbox_inches="tight",
                format="svg",
                transparent=True,
            )
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
        cmap_here="RdBu_r",
        save_name="top50_most_changing_pathways",
        colorbar_title="Sum of NES",
    )
    figures["top50_changing"] = fig1

    # 2. Top 50 most upregulated pathways (average NES > 0)
    logger.info("Creating heatmap 2/5: Top 50 most upregulated pathways...")
    avg_nes = heatmap_data.mean(axis=1)
    upregulated = avg_nes[avg_nes > 0].nlargest(50).index.tolist()

    df_upregulated = heatmap_data.loc[upregulated]
    fig2 = plot_heatmap_regulation(
        df_upregulated,
        upregulated,
        cmap_here="YlGn",
        save_name="top50_most_upregulated_pathways",
        colorbar_title="Sum of NES",
    )
    figures["top50_upregulated"] = fig2

    # 3. Top 50 most downregulated pathways (average NES < 0)
    logger.info("Creating heatmap 3/5: Top 50 most downregulated pathways...")
    downregulated = avg_nes[avg_nes < 0].nsmallest(50).index.tolist()

    df_downregulated = heatmap_data.loc[downregulated]
    fig3 = plot_heatmap_regulation(
        df_downregulated,
        downregulated,
        cmap_here="YlOrBr",
        save_name="top50_most_downregulated_pathways",
        colorbar_title="Sum of NES",
    )
    figures["top50_downregulated"] = fig3

    # 4. High-level pathways (29 pathways from Reactome highest level)
    logger.info("Creating heatmap 4/5: High-level pathways...")
    available_highest = [p for p in highest_pathways if p in heatmap_data.index]

    if available_highest:
        df_highest = heatmap_data.loc[available_highest]
        fig4 = plot_heatmap_regulation(
            df_highest,
            available_highest,
            cmap_here="RdBu_r",
            save_name="selected_pathways_highest_level",
            colorbar_title="Sum of NES",
        )
        figures["selected_highest_level"] = fig4
        logger.info(
            f"Found {len(available_highest)}/{len(highest_pathways)} high-level pathways in data"
        )
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
            cmap_here="RdBu_r",
            save_name="selected_pathways_literature",
            colorbar_title="Sum of NES",
        )
        figures["selected_literature"] = fig5
        logger.info(
            f"Found {len(available_literature)}/{len(pathways_literature)} literature pathways in data"
        )
    else:
        logger.warning("No literature pathways found in the data")

    logger.info(
        f"Pathway heatmap generation complete. Created {len(figures)} heatmaps."
    )

    return heatmap_data, figures
