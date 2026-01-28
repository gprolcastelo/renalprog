"""
Model prediction and trajectory generation functionality for renalprog.

This module will contain:
- VAE inference functions
- Synthetic trajectory generation
- Patient connection identification
- Control trajectory generation
- Trajectory classification
- SDMetrics evaluation

NOTE: This is a stub file. Full implementation requires migrating code from:
- src_deseq_and_gsea_NCSR/synthetic_data_generation.py (trajectory generation)
- src/data/fun_interpol.py (interpolation functions)
- notebooks/4_1_trajectories.ipynb (patient connections)
- src_deseq_and_gsea_NCSR/sdanalyses.py (reconstruction evaluation)
"""

import torch
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
import logging
from sklearn.preprocessing import MinMaxScaler


import os
from tqdm import tqdm
import plotly.express as px
from sdv.io.local import CSVHandler
from sdv.metadata import Metadata
from sdmetrics.reports.single_table import QualityReport
from sdmetrics.single_column import KSComplement
from sdmetrics.single_column import BoundaryAdherence

logger = logging.getLogger(__name__)


def apply_vae(
    model: torch.nn.Module, data: pd.DataFrame, device: str = "cpu"
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Apply VAE to data to get reconstruction and latent representation.

    Args:
        model: Trained VAE model
        data: Input data (samples x genes)
        device: Device to run inference on

    Returns:
        Tuple of (reconstruction, mu, logvar, z)
    """
    model.eval()
    model = model.to(device)

    # Convert to tensor
    if isinstance(data, pd.DataFrame):
        data_tensor = torch.tensor(data.values, dtype=torch.float32).to(device)
    else:
        data_tensor = torch.tensor(data, dtype=torch.float32).to(device)

    with torch.no_grad():
        reconstruction, mu, logvar, z = model(data_tensor)

    # Convert back to numpy
    reconstruction = reconstruction.cpu().numpy()
    mu = mu.cpu().numpy()
    logvar = logvar.cpu().numpy() if logvar is not None else None
    z = z.cpu().numpy()

    return reconstruction, mu, logvar, z




def get_metadata(test_path: Path) -> Dict[str, Any]:
    """
    Extract metadata from test dataset in SDMetrics format.

    This function loads a CSV file and extracts its column metadata using SDMetrics'
    automatic detection. The metadata describes the structure and data types of
    the dataset, which is required for SDMetrics quality evaluation.

    Args:
        test_path: Path to the CSV file containing the test dataset.
                   Can be a string or Path object.

    Returns:
        Dictionary containing metadata with column names and data types.
        Format: {'columns': {col_name: {'sdtype': type}}}

    Note:
        - The CSV is loaded with index_col=0 to avoid treating the index as a feature
        - Both real and synthetic data must share the same metadata structure
        - This ensures SDMetrics can properly validate and compare the datasets

    Example:
        >>> metadata = get_metadata("data/X_test.csv")
        >>> print(metadata['columns'].keys())
        dict_keys(['gene1', 'gene2', ...])
    """
    from pathlib import Path as pathlib_Path

    logger.info("Extracting metadata for SDMetrics evaluation")
    logger.info(f"Loading data from: {test_path}")

    # Convert to Path object for consistent handling across platforms
    test_path = pathlib_Path(test_path)

    # Load data using SDMetrics CSV handler
    # CRITICAL: index_col=0 prevents the index from being treated as a feature column
    # This would cause metadata mismatch errors if the index is included
    connector = CSVHandler()
    real_data = connector.read(
        folder_name=str(test_path.parent),
        file_names=[test_path.name],
        read_csv_parameters={
            "index_col": 0,  # Use first column as index, not as feature
            "parse_dates": False,  # Don't parse dates (all numeric gene expression)
            "encoding": "latin-1",  # Standard encoding for CSV files
        },
    )

    # Auto-detect metadata from the loaded data
    metadata = Metadata.detect_from_dataframes(data=real_data)

    # Extract table-specific metadata (removes wrapper structure)
    # The key 'X_test' matches the filename without extension
    metadata_use = metadata.to_dict()["tables"]["X_test"]

    logger.info(f"Extracted metadata for {len(metadata_use['columns'])} genes")

    return metadata_use


def diagnostic_metrics(
    real_data: pd.DataFrame,
    synthetic_data: pd.DataFrame,
    save_path_data: Path,
    save_path_figures: Optional[Path] = None,
) -> pd.Series:
    """
    Calculate diagnostic metrics to assess synthetic data quality.

    This function computes the Boundary Adherence metric for each gene, which
    measures whether synthetic values respect the min/max boundaries of real data.
    This is a critical diagnostic to detect mode collapse or distribution shift.

    Args:
        real_data: Real gene expression data (samples × genes)
        synthetic_data: Synthetic/reconstructed gene expression data (samples × genes)
        save_path_data: Directory path to save metric CSV results
        save_path_figures: Optional directory path to save visualization plots

    Returns:
        Series with boundary adherence scores per gene (index=gene, values=scores).
        Scores range from 0.0 (worst) to 1.0 (best).

    Metric Details:
        Boundary Adherence per Gene:
        - 1.0 (best): All synthetic values are within [min, max] of real data
        - 0.0 (worst): No synthetic values fall within real data boundaries
        - Values between 0-1 indicate partial adherence

    Saves:
        - {save_path_data}/boundary_adherence_per_gene.csv: Per-gene scores
        - {save_path_figures}/boundary_adherence_per_gene.html: Interactive histogram
        - {save_path_figures}/boundary_adherence_per_gene.{png,pdf,svg}: Static plots

    Example:
        >>> ba_scores = diagnostic_metrics(X_real, X_synthetic, "results/", "figures/")
        >>> print(f"Mean adherence: {ba_scores.mean():.4f}")
        Mean adherence: 0.9823
    """

    logger.info("=" * 60)
    logger.info("DIAGNOSTIC METRICS: Boundary Adherence")
    logger.info("=" * 60)
    logger.info(
        f"Evaluating {real_data.shape[1]} genes across {real_data.shape[0]} samples"
    )

    # Calculate boundary adherence for each gene
    # This measures what percentage of synthetic values fall within the
    # [min, max] range observed in the real data
    ba_dict = {}

    for gene_i in tqdm(real_data.columns, desc="Computing Boundary Adherence"):
        # Compute metric: % of synthetic values within [min, max] of real values
        ba_i = BoundaryAdherence.compute(
            real_data=real_data[gene_i], synthetic_data=synthetic_data[gene_i]
        )
        ba_dict[gene_i] = ba_i

    # Convert to Series for easy analysis
    df_ba = pd.Series(ba_dict, name="boundary_adherence")

    # Save results
    output_csv = os.path.join(save_path_data, "boundary_adherence_per_gene.csv")
    df_ba.to_csv(output_csv)
    logger.info(f"Saved results to: {output_csv}")

    # Log summary statistics
    logger.info(f"Mean Boundary Adherence: {df_ba.mean():.4f}")
    logger.info(f"Median Boundary Adherence: {df_ba.median():.4f}")
    logger.info(f"Min Boundary Adherence: {df_ba.min():.4f}")
    logger.info(f"Max Boundary Adherence: {df_ba.max():.4f}")
    logger.info(
        f"Genes with perfect adherence (1.0): {(df_ba == 1.0).sum()}/{len(df_ba)} ({100 * (df_ba == 1.0).sum() / len(df_ba):.1f}%)"
    )

    # Generate visualizations if output directory provided
    if save_path_figures is not None:
        logger.info("Generating visualizations...")

        # Create interactive histogram
        fig = px.histogram(
            df_ba,
            x="boundary_adherence",
            nbins=50,
            title="Distribution of Boundary Adherence Scores per Gene",
            labels={"boundary_adherence": "Boundary Adherence Score"},
            template="plotly_white",
        )
        fig.update_layout(
            xaxis_title="Boundary Adherence Score",
            yaxis_title="Number of Genes",
            showlegend=False,
        )

        # Save in multiple formats
        html_path = os.path.join(save_path_figures, "boundary_adherence_per_gene.html")
        fig.write_html(html_path)
        logger.info(f"  Saved interactive plot: {html_path}")

        for format_ext in ["png", "pdf", "svg"]:
            img_path = os.path.join(
                save_path_figures, f"boundary_adherence_per_gene.{format_ext}"
            )
            fig.write_image(img_path, scale=2)
            logger.info(f"  Saved {format_ext.upper()} plot: {img_path}")

    logger.info("Diagnostic metrics calculation complete")
    logger.info("=" * 60)

    return df_ba


def quality_metrics(
    real_data: pd.DataFrame,
    synthetic_data: pd.DataFrame,
    metadata: Dict[str, Any],
    save_path_data: Path,
    save_path_figures: Optional[Path] = None,
) -> pd.Series:
    """
    Calculate quality metrics to assess synthetic data fidelity.

    This function computes two key metrics:
    1. Quality Report: Overall assessment of column shapes and pair-wise trends
    2. KS Complement: Per-gene similarity of marginal distributions

    These metrics evaluate how well the synthetic data captures the statistical
    properties of the real data, beyond just staying within boundaries.

    Args:
        real_data: Real gene expression data (samples × genes)
        synthetic_data: Synthetic/reconstructed gene expression data (samples × genes)
        metadata: Metadata dictionary from get_metadata() for SDMetrics
        save_path_data: Directory path to save metric results
        save_path_figures: Optional directory path to save visualization plots

    Returns:
        Series with KS Complement scores per gene (index=gene, values=scores).
        Scores range from 0.0 (worst) to 1.0 (best).

    Metrics Details:
        Column Shapes (in Quality Report):
        - Measures overall distribution similarity per column
        - Higher scores indicate better shape matching

        Column Pair Trends (in Quality Report):
        - Measures correlation and relationship preservation
        - Higher scores indicate better trend matching

        KS Complement (per gene):
        - 1.0 (best): Real and synthetic distributions are identical
        - 0.0 (worst): Distributions are maximally different
        - Based on Kolmogorov-Smirnov test

    Saves:
        - {save_path_data}/quality_report.pkl: Full SDMetrics quality report
        - {save_path_data}/ks_complement_per_gene.csv: Per-gene KS scores
        - {save_path_figures}/ks_complement_per_gene.html: Interactive histogram
        - {save_path_figures}/ks_complement_per_gene.{png,pdf,svg}: Static plots

    Example:
        >>> ks_scores = quality_metrics(X_real, X_synth, metadata, "results/", "figs/")
        >>> print(f"Mean KS Complement: {ks_scores.mean():.4f}")
        Mean KS Complement: 0.8756
    """
    import os
    from tqdm import tqdm
    import plotly.express as px

    logger.info("=" * 60)
    logger.info("QUALITY METRICS: Distribution Similarity")
    logger.info("=" * 60)

    # Generate comprehensive quality report
    # This evaluates:
    # 1. Column Shapes: How well distributions match per gene
    # 2. Column Pair Trends: How well correlations are preserved
    logger.info("Generating SDMetrics Quality Report...")
    q_report = QualityReport()
    q_report.generate(real_data, synthetic_data, metadata)

    # Save quality report object for later analysis
    report_path = os.path.join(save_path_data, "quality_report.pkl")
    q_report.save(report_path)
    logger.info(f"Saved quality report to: {report_path}")

    # Get overall quality score from the report
    overall_score = q_report.get_score()
    logger.info(f"Overall Quality Score: {overall_score:.4f}")

    # Calculate KS Complement for each gene
    # This measures similarity of marginal distributions (1D histograms)
    # KS Complement = 1 - KS statistic, where KS statistic measures max difference between CDFs
    logger.info(f"Computing KS Complement for {real_data.shape[1]} genes...")

    ks_dict = {}
    for gene_i in tqdm(real_data.columns, desc="Computing KS Complement"):
        # KS Complement measures how similar the empirical cumulative distribution functions are
        # Higher values mean the distributions are more similar
        ks_i = KSComplement.compute(
            real_data=real_data[gene_i], synthetic_data=synthetic_data[gene_i]
        )
        ks_dict[gene_i] = ks_i

    # Convert to Series for analysis
    df_ks = pd.Series(ks_dict, name="ks_complement")

    # Save results
    output_csv = os.path.join(save_path_data, "ks_complement_per_gene.csv")
    df_ks.to_csv(output_csv)
    logger.info(f"Saved results to: {output_csv}")

    # Log summary statistics
    logger.info(f"Mean KS Complement: {df_ks.mean():.4f}")
    logger.info(f"Median KS Complement: {df_ks.median():.4f}")
    logger.info(f"Min KS Complement: {df_ks.min():.4f}")
    logger.info(f"Max KS Complement: {df_ks.max():.4f}")
    logger.info(
        f"Genes with KS > 0.9: {(df_ks > 0.9).sum()}/{len(df_ks)} ({100 * (df_ks > 0.9).sum() / len(df_ks):.1f}%)"
    )

    # Generate visualizations if output directory provided
    if save_path_figures is not None:
        logger.info("Generating visualizations...")

        # Create interactive histogram
        fig = px.histogram(
            df_ks,
            x="ks_complement",
            nbins=50,
            title="Distribution of KS Complement Scores per Gene",
            labels={"ks_complement": "KS Complement Score"},
            template="plotly_white",
        )
        fig.update_layout(
            xaxis_title="KS Complement Score (Distribution Similarity)",
            yaxis_title="Number of Genes",
            showlegend=False,
        )

        # Add reference line at 0.9 (high quality threshold)
        fig.add_vline(
            x=0.9,
            line_dash="dash",
            line_color="red",
            annotation_text="High Quality (0.9)",
        )

        # Save in multiple formats
        html_path = os.path.join(save_path_figures, "ks_complement_per_gene.html")
        fig.write_html(html_path)
        logger.info(f"  Saved interactive plot: {html_path}")

        for format_ext in ["png", "pdf", "svg"]:
            img_path = os.path.join(
                save_path_figures, f"ks_complement_per_gene.{format_ext}"
            )
            fig.write_image(img_path, scale=2)
            logger.info(f"  Saved {format_ext.upper()} plot: {img_path}")

    logger.info("Quality metrics calculation complete")
    logger.info("=" * 60)

    return df_ks


def evaluate_reconstruction(
    real_data: pd.DataFrame,
    synthetic_data: pd.DataFrame,
    save_path_data: Path,
    save_path_figures: Optional[Path] = None,
    metadata_path: Path = None,
) -> Tuple[pd.Series, pd.Series]:
    """
    Comprehensive evaluation of reconstruction quality using SDMetrics.

    This function orchestrates a complete quality assessment of synthetic/reconstructed
    data by computing both diagnostic and quality metrics. It's the main entry point
    for evaluating VAE reconstructions or VAE+RecNet outputs.

    The evaluation includes:
    1. Boundary Adherence: Do synthetic values stay within real data bounds?
    2. Distribution Similarity: Do synthetic distributions match real distributions?
    3. Quality Report: Overall assessment of column shapes and correlations

    Args:
        real_data: Real gene expression data (samples × genes)
        synthetic_data: Synthetic/reconstructed gene expression data (samples × genes)
        save_path_data: Directory path to save all metric results (CSV, PKL)
        save_path_figures: Optional directory path to save visualization plots
        metadata_path: Path to CSV file used to extract metadata structure
                       (typically the test set CSV file)

    Returns:
        Tuple of (boundary_adherence_series, ks_complement_series):
        - boundary_adherence_series: Series with boundary adherence scores per gene
        - ks_complement_series: Series with KS Complement scores per gene

    Workflow:
        1. Extract metadata from test CSV file
        2. Compute diagnostic metrics (boundary adherence)
        3. Compute quality metrics (KS complement + quality report)
        4. Save all results and visualizations

    Output Files:
        In save_path_data/:
        - boundary_adherence_per_gene.csv: Per-gene boundary scores
        - ks_complement_per_gene.csv: Per-gene distribution similarity
        - quality_report.pkl: Full SDMetrics quality report object

        In save_path_figures/ (if provided):
        - boundary_adherence_per_gene.{html,png,pdf,svg}
        - ks_complement_per_gene.{html,png,pdf,svg}

    Interpretation:
        - Higher scores are better for both metrics (range: 0.0 to 1.0)
        - Boundary Adherence: Checks if synthetic data stays in valid ranges
        - KS Complement: Checks if distributions match (more stringent)
        - Good reconstruction: BA > 0.95, KS > 0.85
        - Excellent reconstruction: BA > 0.99, KS > 0.90

    Example:
        >>> ba_scores, ks_scores = evaluate_reconstruction(
        ...     real_data=X_test,
        ...     synthetic_data=vae_reconstruction,
        ...     save_path_data="results/vae_eval/",
        ...     save_path_figures="figures/vae_eval/",
        ...     metadata_path="data/X_test.csv"
        ... )
        >>> print(f"Mean BA: {ba_scores.mean():.4f}, Mean KS: {ks_scores.mean():.4f}")
        Mean BA: 0.9823, Mean KS: 0.8756

    Note:
        - Ensure real_data and synthetic_data have identical column names and order
        - The metadata_path CSV should have the same structure as real_data
        - This function is used by scripts/pipeline_steps/3_check_reconstruction.py
        - Both DataFrames should have samples as rows and genes as columns

    Raises:
        ValueError: If data shapes don't match or columns don't align
        FileNotFoundError: If metadata_path doesn't exist
    """
    logger.info("=" * 80)
    logger.info("RECONSTRUCTION QUALITY EVALUATION")
    logger.info("=" * 80)
    logger.info(f"Real data shape: {real_data.shape}")
    logger.info(f"Synthetic data shape: {synthetic_data.shape}")
    logger.info(f"Saving results to: {save_path_data}")
    if save_path_figures:
        logger.info(f"Saving figures to: {save_path_figures}")

    # Validate inputs
    if real_data.shape != synthetic_data.shape:
        raise ValueError(
            f"Data shape mismatch: real {real_data.shape} vs synthetic {synthetic_data.shape}"
        )

    if not all(real_data.columns == synthetic_data.columns):
        raise ValueError("Column names must match between real and synthetic data")

    # Step 1: Extract metadata in SDMetrics format
    logger.info("\n[1/3] Extracting metadata...")
    metadata_sd = get_metadata(metadata_path)
    logger.info(f"Metadata extracted for {len(metadata_sd['columns'])} genes")

    # Step 2: Compute diagnostic metrics
    logger.info("\n[2/3] Computing diagnostic metrics...")
    df_ba = diagnostic_metrics(
        real_data=real_data,
        synthetic_data=synthetic_data,
        save_path_data=save_path_data,
        save_path_figures=save_path_figures,
    )

    # Step 3: Compute quality metrics
    logger.info("\n[3/3] Computing quality metrics...")
    df_ks = quality_metrics(
        real_data=real_data,
        synthetic_data=synthetic_data,
        metadata=metadata_sd,
        save_path_data=save_path_data,
        save_path_figures=save_path_figures,
    )

    # Final summary
    logger.info("\n" + "=" * 80)
    logger.info("EVALUATION COMPLETE - SUMMARY")
    logger.info("=" * 80)
    logger.info(
        f"Boundary Adherence - Mean: {df_ba.mean():.4f}, Median: {df_ba.median():.4f}"
    )
    logger.info(
        f"KS Complement      - Mean: {df_ks.mean():.4f}, Median: {df_ks.median():.4f}"
    )

    # Quality assessment
    ba_quality = (
        "Excellent"
        if df_ba.mean() > 0.99
        else "Good"
        if df_ba.mean() > 0.95
        else "Fair"
        if df_ba.mean() > 0.90
        else "Poor"
    )
    ks_quality = (
        "Excellent"
        if df_ks.mean() > 0.90
        else "Good"
        if df_ks.mean() > 0.85
        else "Fair"
        if df_ks.mean() > 0.75
        else "Poor"
    )

    logger.info(
        f"Overall Assessment - Boundary: {ba_quality}, Distribution: {ks_quality}"
    )
    logger.info(f"Results saved to: {save_path_data}")
    logger.info("=" * 80)

    return df_ba, df_ks


def classify_trajectories(
    classifier,
    trajectory_data: Dict[str, pd.DataFrame],
    gene_subset: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Apply stage classifier to synthetic trajectories.

    Args:
        classifier: Trained classifier model
        trajectory_data: Dictionary of patient pair to trajectory DataFrames
        gene_subset: Optional subset of genes to use for classification

    Returns:
        DataFrame with classification results for each trajectory point
    """
    logger.info("Classifying trajectory points")

    # TODO: Implement trajectory classification
    # Migrate from notebooks/kirc_classification_trajectory.ipynb

    raise NotImplementedError(
        "classify_trajectories() needs implementation from "
        "notebooks/kirc_classification_trajectory.ipynb"
    )


def interpolate_latent_linear(
    z_source: np.ndarray, z_target: np.ndarray, n_steps: int = 50
) -> np.ndarray:
    """
    Linear interpolation in latent space.

    Args:
        z_source: Source latent vector
        z_target: Target latent vector
        n_steps: Number of interpolation steps

    Returns:
        Array of interpolated latent vectors (n_steps x latent_dim)
    """
    alphas = np.linspace(0, 1, n_steps)
    interpolated = np.array(
        [(1 - alpha) * z_source + alpha * z_target for alpha in alphas]
    )
    return interpolated


def interpolate_latent_spherical(
    z_source: np.ndarray, z_target: np.ndarray, n_steps: int = 50
) -> np.ndarray:
    """
    Spherical (SLERP) interpolation in latent space.

    Args:
        z_source: Source latent vector
        z_target: Target latent vector
        n_steps: Number of interpolation steps

    Returns:
        Array of interpolated latent vectors (n_steps x latent_dim)
    """
    # Normalize vectors
    z_source_norm = z_source / np.linalg.norm(z_source)
    z_target_norm = z_target / np.linalg.norm(z_target)

    # Calculate angle between vectors
    omega = np.arccos(np.clip(np.dot(z_source_norm, z_target_norm), -1.0, 1.0))

    if omega < 1e-8:
        # Vectors are nearly identical, use linear interpolation
        return interpolate_latent_linear(z_source, z_target, n_steps)

    # SLERP formula
    alphas = np.linspace(0, 1, n_steps)
    interpolated = np.array(
        [
            (np.sin((1 - alpha) * omega) / np.sin(omega)) * z_source
            + (np.sin(alpha * omega) / np.sin(omega)) * z_target
            for alpha in alphas
        ]
    )

    return interpolated


def dynamic_enrichment_analysis(
    trajectory_dir: Path,
    pathways_file: Path,
    output_dir: Path,
    cancer_type: str = "kirc",
) -> pd.DataFrame:
    """
    Perform dynamic enrichment analysis on synthetic trajectories.

    This orchestrates:
    1. DESeq2 analysis on each trajectory point
    2. GSEA on differential expression results
    3. Aggregation of enrichment across trajectories

    Args:
        trajectory_dir: Directory containing trajectory CSV files
        pathways_file: Path to pathway GMT file
        output_dir: Directory to save results
        cancer_type: Cancer type identifier

    Returns:
        DataFrame with aggregated enrichment results
    """
    logger.info(f"Running dynamic enrichment analysis for {cancer_type}")

    # TODO: Implement orchestration
    # Migrate from src_deseq_and_gsea_NCSR/pipeline.sh and related scripts

    raise NotImplementedError(
        "dynamic_enrichment_analysis() needs implementation. "
        "Migrate orchestration from src_deseq_and_gsea_NCSR/pipeline.sh, "
        "py_deseq.py, and trajectory_formatting.py"
    )


# =============================================================================
# Patient Connection and Trajectory Generation Functions
# =============================================================================


def calculate_all_possible_transitions(
    data: pd.DataFrame,
    metadata_selection: pd.DataFrame,
    distance: str = "wasserstein",
    early_late: bool = False,
    negative_trajectory: bool = False,
) -> pd.DataFrame:
    """
    Calculate all possible patient-to-patient transitions for KIRC cancer.

    This function computes pairwise distances between all patients at consecutive
    (or same) cancer stages, considering metadata constraints. Only patients with
    matching gender and race are considered as potential trajectory pairs.

    Parameters
    ----------
    data : pd.DataFrame
        Gene expression data with patients as columns.
    metadata_selection : pd.DataFrame
        Clinical metadata with columns: histological_type, race, gender, stage.
    distance : {'wasserstein', 'euclidean'}, default='wasserstein'
        Distance metric to use for calculating patient-to-patient distances.
    early_late : bool, default=False
        If True, uses early/late stage groupings. If False, uses I-IV stages.
    negative_trajectory : bool, default=False
        If True, generates same-stage transitions (negative controls).
        If False, generates progression transitions (positive trajectories).

    Returns
    -------
    pd.DataFrame
        DataFrame containing all possible transitions with columns:
        - source, target: Patient IDs
        - source_gender, target_gender: Gender
        - source_race, target_race: Race
        - transition: Stage transition label (e.g., '1_to_2', 'early_to_late')
        - distance: Calculated distance between patients

        Sorted by gender, race, transition, and distance.

    Raises
    ------
    ValueError
        If distance metric is not 'wasserstein' or 'euclidean'.

    Notes
    -----
    - For positive trajectories: links I→II, II→III, III→IV or early→late
    - For negative trajectories: links I→I, II→II, III→III, IV→IV or early→early, late→late
    - Only patients with identical gender and race are paired
    """
    # Select distance function
    if distance == "wasserstein":
        from scipy.stats import wasserstein_distance

        distance_fun = wasserstein_distance
    elif distance == "euclidean":
        from scipy.spatial.distance import euclidean

        distance_fun = euclidean
    else:
        raise ValueError(
            'Distance function not implemented. Use either "wasserstein" or "euclidean".'
        )

    # Define stage transitions based on parameters
    if early_late and not negative_trajectory:
        possible_transitions = ["early_to_late"]
        stage_pairs = [["early", "late"]]
    elif early_late and negative_trajectory:
        possible_transitions = ["early_to_early", "late_to_late"]
        stage_pairs = [["early", "early"], ["late", "late"]]
    elif not early_late and not negative_trajectory:
        possible_transitions = ["1_to_2", "2_to_3", "3_to_4"]
        stage_pairs = [["I", "II"], ["II", "III"], ["III", "IV"]]
    elif not early_late and negative_trajectory:
        possible_transitions = ["1_to_1", "2_to_2", "3_to_3", "4_to_4"]
        stage_pairs = [["I", "I"], ["II", "II"], ["III", "III"], ["IV", "IV"]]

    # Calculate all possible transitions
    results = []
    for i_tr, transition in enumerate(possible_transitions):
        source_target_stage = stage_pairs[i_tr]

        # Iterate through all patient pairs at specified stages
        for pat_i in metadata_selection.index[
            metadata_selection["stage"] == source_target_stage[0]
        ]:
            for pat_ii in metadata_selection.index[
                metadata_selection["stage"] == source_target_stage[1]
            ]:
                # Extract metadata for both patients
                source_gender = metadata_selection.at[pat_i, "gender"]
                target_gender = metadata_selection.at[pat_ii, "gender"]
                source_race = metadata_selection.at[pat_i, "race"]
                target_race = metadata_selection.at[pat_ii, "race"]

                # Skip if metadata doesn't match (gender and race must match)
                if not (source_race == target_race and source_gender == target_gender):
                    continue

                # Store transition information
                results_i = {
                    "source": pat_i,
                    "target": pat_ii,
                    "source_gender": source_gender,
                    "target_gender": target_gender,
                    "source_race": source_race,
                    "target_race": target_race,
                    "transition": transition,
                    "distance": distance_fun(data[pat_i], data[pat_ii]),
                }
                results.append(results_i)

    # Convert to DataFrame and sort
    results_df = pd.DataFrame(results)
    results_df.sort_values(
        [
            "source_gender",
            "target_gender",
            "source_race",
            "target_race",
            "transition",
            "distance",
        ],
        inplace=True,
        ignore_index=True,
    )
    return results_df


def link_patients_closest(
    transitions_df: pd.DataFrame,
    start_with_first_stage: bool = True,
    early_late: bool = False,
    closest: bool = True,
) -> pd.DataFrame:
    """
    Link patients by selecting closest (or farthest) matches across stages.

    For each patient at a source stage, this function identifies the closest
    (or farthest) patient at the target stage, considering metadata constraints
    (gender, race). This creates one-to-one patient linkages that form the basis
    for trajectory construction.

    Args:
        transitions_df: DataFrame from calculate_all_possible_transitions()
                        containing all possible patient pairs with distances
        start_with_first_stage: If True, build forward trajectories (early→late)
                                If False, build backward trajectories (late→early)
        early_late: If True, uses early/late groupings. If False, uses I-IV stages
        closest: If True, connect closest patients. If False, connect farthest patients

    Returns:
        DataFrame with selected patient links, containing one row per source patient
        with their optimal target patient match. Includes all columns from transitions_df.

    Selection Strategy:
        - Forward (start_with_first_stage=True): For each source, find optimal target
        - Backward (start_with_first_stage=False): For each target, find optimal source
        - Closest (closest=True): Minimum distance match
        - Farthest (closest=False): Maximum distance match

    Metadata Stratification:
        Links are selected independently within each combination of:
        - Gender (MALE, FEMALE)
        - Race (ASIAN, BLACK OR AFRICAN AMERICAN, WHITE)
        This ensures demographic consistency in trajectories.

    Example:
        >>> links = link_patients_closest(
        ...     transitions_df=all_transitions,
        ...     start_with_first_stage=True,
        ...     closest=True
        ... )
        >>> print(f"Created {len(links)} patient links")
        Created 234 patient links

    Note:
        - Processes transitions in order for forward: I→II→III→IV
        - Processes in reverse for backward: IV→III→II→I
        - Each patient appears at most once as a source in the result
    """
    logger.info("Linking patients by closest/farthest matches")
    logger.info(f"Direction: {'Forward' if start_with_first_stage else 'Backward'}")
    logger.info(f"Strategy: {'Closest' if closest else 'Farthest'}")

    # Define transition order based on direction
    if start_with_first_stage and not early_late:
        transitions_possible = ["1_to_2", "2_to_3", "3_to_4"]
    elif not start_with_first_stage and not early_late:
        transitions_possible = ["3_to_4", "2_to_3", "1_to_2"]
    elif early_late:
        transitions_possible = ["early_to_late"]

    # 0 for closest (smallest distance), -1 for farthest (largest distance)
    idx = 0 if closest else -1

    # Find closest/farthest patient for each source patient
    closest_list = []
    for transition_i in transitions_possible:
        transition_df_i = transitions_df[transitions_df["transition"] == transition_i]

        logger.info(
            f"Processing transition {transition_i}: {len(transition_df_i)} pairs"
        )

        # Iterate through all metadata combinations
        for gender_i in ["FEMALE", "MALE"]:
            df_gender_i = transition_df_i.query(f"source_gender == '{gender_i}'")

            for race_i in ["ASIAN", "BLACK OR AFRICAN AMERICAN", "WHITE"]:
                df_race_i = df_gender_i.query(f"source_race == '{race_i}'")

                if df_race_i.empty:
                    continue

                # Get unique patients to link
                unique_sources = df_race_i["source"].unique()
                unique_targets = df_race_i["target"].unique()
                use_uniques = (
                    unique_sources if start_with_first_stage else unique_targets
                )
                use_column = "source" if start_with_first_stage else "target"

                # Find closest/farthest match for each patient
                for pat_i in use_uniques:
                    pat_matches = df_race_i[df_race_i[use_column] == pat_i]
                    if len(pat_matches) > 0:
                        # Sort by distance and select first (closest) or last (farthest)
                        best_match = pat_matches.sort_values("distance").iloc[idx]
                        closest_list.append(best_match)

    # Convert to DataFrame
    closest_df = pd.DataFrame(closest_list)
    closest_df.reset_index(drop=True, inplace=True)

    logger.info(f"Created {len(closest_df)} patient links")

    return closest_df


def link_patients_random(
    results_df: pd.DataFrame,
    start_with_first_stage: bool = True,
    link_next: int = 5,
    transitions_possible: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Link patients to multiple random targets at the next stage.

    Instead of linking each patient to only their closest match, this function randomly
    samples multiple patients at the next stage to link to each source patient. This
    creates a one-to-many mapping useful for generating multiple trajectory samples.

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame with possible sources and targets, their metadata, and distance.
    start_with_first_stage : bool, default=True
        If True, initiate trajectories with first stage as sources.
        If False, initiate trajectories with last stage as sources.
    link_next : int, default=5
        Number of patients at next stage to randomly link to each patient of current stage.
    transitions_possible : list, optional
        List of transitions to process (e.g., ['1_to_2', '2_to_3']).
        If None, defaults to ['early_to_late'].

    Returns
    -------
    pd.DataFrame
        DataFrame with randomly sampled patient links for each transition.
        Contains multiple rows per source patient (up to link_next).

    Notes
    -----
    - Random sampling is primarily performed for WHITE race patients due to sample size
    - If fewer than link_next targets are available, all available targets are selected
    - Patients from other races are included with all their possible connections
    - Empty DataFrame is returned if no WHITE patients are found
    """
    # Set default transitions if not provided
    if transitions_possible is None:
        transitions_possible = ["early_to_late"]

    # Get unique genders and races
    unique_genders = results_df["source_gender"].unique().tolist()
    # Get unique races
    unique_races = results_df["source_race"].unique().tolist()
    if "WHITE" in unique_races:
        unique_races.remove("WHITE")
    # transition:
    samples = []
    for transition_i in transitions_possible:
        transition_df_i = results_df[results_df["transition"] == transition_i]
        for gender_i in unique_genders:
            df_samples_i = transition_df_i.query(
                f"source_gender == '{gender_i}' & source_race == 'WHITE'"
            )  # we can only do this for the whites since these are the only ones with enough samples
            if df_samples_i.empty:
                print(
                    f"Warning: No WHITE patients found for gender {gender_i} in transition {transition_i}"
                )
                continue
            unique_sources_i = np.unique(df_samples_i["source"]).tolist()
            unique_targets_i = np.unique(df_samples_i["target"]).tolist()
            use_uniques = (
                unique_sources_i if start_with_first_stage else unique_targets_i
            )
            use_source_target = "source" if start_with_first_stage else "target"
            for pat_i in use_uniques:
                sample_i = df_samples_i.loc[df_samples_i[use_source_target] == pat_i]
                if len(sample_i) >= link_next:
                    sample_i = sample_i.sample(
                        link_next
                    )  # Sample a number of patients at next stage to link to each patient of current stage
                else:
                    sample_i = sample_i.sample(
                        len(sample_i)
                    )  # Sample all available patients if less than link_next
                samples.append(sample_i)

    # Check if samples list is empty
    if not samples:
        print("Warning: No samples found for WHITE race. Returning empty DataFrame.")
        return pd.DataFrame(columns=results_df.columns)

    # Turn samples into dataframe:
    samples_df = pd.concat(samples)
    # Add the rest of the races
    if unique_races:
        samples_df = pd.concat(
            [samples_df, results_df[results_df["source_race"].isin(unique_races)]]
        )
    samples_df.reset_index(drop=True, inplace=True)
    return samples_df


def build_trajectory_network(
    patient_links: pd.DataFrame,
) -> Tuple[Dict[str, List[str]], List[List[str]]]:
    """
    Build trajectory network and find all complete disease progression paths.

    Constructs a directed graph from patient links and identifies all possible
    complete trajectories from root nodes (earliest stage patients not appearing
    as targets) to leaf nodes (latest stage patients not appearing as sources).

    Args:
        patient_links: DataFrame with 'source' and 'target' columns from linking functions

    Returns:
        Tuple of:
        - network: Dict mapping each source patient to list of target patients
        - trajectories: List of complete trajectories, where each trajectory is a
                        list of patient IDs ordered from earliest to latest stage

    Network Structure:
        - Adjacency list representation: {source: [target1, target2, ...]}
        - Directed edges from earlier to later stages
        - Allows multiple outgoing edges (one patient → multiple next-stage patients)

    Trajectory Discovery:
        - Uses depth-first search from root nodes
        - Root nodes: Patients in 'source' but not in 'target' (stage I or early)
        - Leaf nodes: Patients in 'target' but not in 'source' (stage IV or late)
        - Each trajectory represents a complete disease progression path

    Example:
        >>> network, trajectories = build_trajectory_network(patient_links)
        >>> print(f"Network has {len(network)} nodes")
        >>> print(f"Found {len(trajectories)} complete trajectories")
        >>> print(f"Example trajectory: {trajectories[0]}")
        Network has 500 nodes
        Found 234 complete trajectories
        Example trajectory: ['PAT001', 'PAT045', 'PAT123', 'PAT289']

    Trajectory Characteristics:
        - Length varies based on how many stages the path spans
        - Typical lengths: 2-4 patients for I→II→III→IV progressions
        - Length 2 for early→late progressions
        - Patients can appear in multiple trajectories

    Note:
        - Cycles are prevented during trajectory search
        - All paths from root to leaf are enumerated
        - Trajectories respect chronological disease progression
    """
    logger.info("Building trajectory network from patient links")

    sources = patient_links["source"]
    targets = patient_links["target"]

    # Build network adjacency list
    network = {}
    for source, target in zip(sources, targets):
        if source not in network:
            network[source] = []
        network[source].append(target)

    logger.info(f"Network built: {len(network)} source nodes")

    # Find root nodes (patients who are sources but never targets)
    unique_sources = set(sources) - set(targets)
    logger.info(f"Found {len(unique_sources)} root nodes (earliest stage patients)")

    # Recursively find all trajectories from each root
    def find_trajectories(
        start_node: str, visited: Optional[List[str]] = None
    ) -> List[List[str]]:
        """Depth-first search to find all paths from start_node to leaf nodes."""
        if visited is None:
            visited = []

        visited.append(start_node)

        # If node has no outgoing edges, this is a leaf node - return path
        if start_node not in network:
            return [visited]

        # Recursively explore all targets
        trajectories = []
        for target in network[start_node]:
            if target not in visited:  # Avoid cycles
                new_visited = visited.copy()
                trajectories.extend(find_trajectories(target, new_visited))

        return trajectories

    # Find all trajectories starting from each root
    all_trajectories = []

    if len(unique_sources) == 0:
        # No clear root nodes - this happens with early→late transitions where
        # patients can be both sources and targets. In this case, each source→target
        # pair is already a complete 2-patient trajectory.
        logger.info("No root nodes found (typical for early→late transitions).")
        logger.info("Using each source→target pair as a complete trajectory.")
        for source, target in zip(sources, targets):
            all_trajectories.append([source, target])
    else:
        # Standard case: multi-stage progressions (I→II→III→IV)
        for source in unique_sources:
            all_trajectories.extend(find_trajectories(source))

    logger.info(
        f"Discovered {len(all_trajectories)} complete disease progression trajectories"
    )

    # Log trajectory length statistics only if we have trajectories
    if len(all_trajectories) > 0:
        traj_lengths = [len(t) for t in all_trajectories]
        logger.info(
            f"Trajectory lengths - Min: {min(traj_lengths)}, Max: {max(traj_lengths)}, "
            f"Mean: {np.mean(traj_lengths):.1f}"
        )
    else:
        logger.warning("No trajectories found!")

    return network, all_trajectories


def generate_trajectory_data(
    vae_model: torch.nn.Module,
    recnet_model: Optional[torch.nn.Module],
    trajectory: List[str],
    gene_data: pd.DataFrame,
    n_timepoints: int = 50,
    interpolation_method: str = "linear",
    device: str = "cpu",
    save_path: Optional[Path] = None,
    scaler: Optional[MinMaxScaler] = None,
) -> pd.DataFrame:
    """
    Generate synthetic gene expression data along a patient trajectory.

    Creates N interpolated time points between consecutive patients in a trajectory
    by performing interpolation in the VAE latent space, then decoding back to
    gene expression space. Optionally applies reconstruction network for refinement.

    Args:
        vae_model: Trained VAE model for encoding/decoding
        recnet_model: Optional reconstruction network for refining VAE output
        trajectory: List of patient IDs in chronological progression order
        gene_data: Gene expression DataFrame (genes × patients)
        n_timepoints: Number of interpolation points between each patient pair
        interpolation_method: 'linear' or 'spherical' interpolation in latent space
        device: Torch device for computation
        save_path: Optional path to save trajectory CSV file
        scaler: Pre-fitted MinMaxScaler from VAE training. If None, will fit on gene_data.

    Returns:
        DataFrame with synthetic gene expression profiles for all time points.
        Shape: (n_timepoints * (len(trajectory)-1), n_genes)
        Index contains time point identifiers

    Workflow:
        1. Extract gene expression for each patient in trajectory
        2. Normalize using the SAME scaler used during VAE training
        3. Encode each patient to VAE latent space
        4. For each consecutive pair:
           a. Interpolate in latent space (linear or spherical)
           b. Decode interpolated points back to gene space
           c. Optionally apply reconstruction network
        5. Concatenate all segments into complete trajectory

    Interpolation Methods:
        linear: Straight-line interpolation in latent space
                z(t) = (1-t)*z_source + t*z_target

        spherical: Spherical linear interpolation (SLERP)
                   Preserves magnitude, interpolates on hypersphere
                   Recommended for normalized latent spaces

    Note:
        CRITICAL: The scaler must be the same one used during VAE training.
        Using a different scaler will produce incorrect latent representations.
        If scaler=None, will fit on all gene_data (all patients), which approximates
        the training distribution.
    """

    logger.info(f"Generating trajectory data for {len(trajectory)} patients")
    logger.info(
        f"Interpolation: {n_timepoints} points × {len(trajectory) - 1} segments"
    )
    logger.info(f"Method: {interpolation_method}")

    # Set models to evaluation mode
    vae_model.eval()
    if recnet_model is not None:
        recnet_model.eval()

    vae_model = vae_model.to(device)
    if recnet_model is not None:
        recnet_model = recnet_model.to(device)

    # Use provided scaler or fit new one on all gene data
    if scaler is None:
        logger.warning("No scaler provided - fitting new scaler on all gene data")
        logger.warning("This may not match VAE training normalization!")
        scaler = MinMaxScaler()
        # gene_data is (genes × patients), need (patients × genes) for scaler
        scaler.fit(gene_data.T.values)
        logger.info(f"Fitted scaler on {gene_data.shape[1]} patients")
    else:
        logger.info("Using provided scaler from VAE training")

    # Select interpolation function
    if interpolation_method == "linear":
        interp_func = interpolate_latent_linear
    elif interpolation_method == "spherical":
        interp_func = interpolate_latent_spherical
    else:
        raise ValueError(f"Unknown interpolation method: {interpolation_method}")

    # Generate synthetic data for each segment of the trajectory
    all_segments = []

    with torch.no_grad():
        for i in range(len(trajectory) - 1):
            source_patient = trajectory[i]
            target_patient = trajectory[i + 1]

            logger.info(
                f"Segment {i + 1}/{len(trajectory) - 1}: {source_patient} → {target_patient}"
            )

            # Get gene expression for source and target
            # gene_data is (genes × patients), so gene_data[patient] is a Series of gene values
            source_expr = gene_data[source_patient].values.reshape(1, -1)  # (1, genes)
            target_expr = gene_data[target_patient].values.reshape(1, -1)  # (1, genes)

            # Normalize data using the provided scaler
            # Scaler expects (n_samples, n_features) = (1, genes)
            source_norm = scaler.transform(source_expr)  # (1, genes)
            target_norm = scaler.transform(target_expr)  # (1, genes)

            # Encode to latent space
            source_tensor = torch.tensor(source_norm, dtype=torch.float32).to(device)
            target_tensor = torch.tensor(target_norm, dtype=torch.float32).to(device)

            _, _, _, z_source = vae_model(source_tensor)
            _, _, _, z_target = vae_model(target_tensor)

            # Interpolate in latent space
            z_source_np = z_source.cpu().numpy().flatten()
            z_target_np = z_target.cpu().numpy().flatten()

            interpolated_z = interp_func(z_source_np, z_target_np, n_timepoints)

            # Decode interpolated latent vectors
            interpolated_z_tensor = torch.tensor(
                interpolated_z, dtype=torch.float32
            ).to(device)
            decoded = vae_model.decoder(interpolated_z_tensor)

            # Denormalize using the same scaler
            # decoded is (n_timepoints, genes), scaler expects (n_samples, n_features)
            decoded_np = decoded.cpu().numpy()  # (n_timepoints, genes)
            segment_data = scaler.inverse_transform(
                decoded_np
            )  # (n_timepoints, genes) - REAL SPACE

            # Apply reconstruction network if provided
            # CRITICAL: RecNet works on REAL SPACE data, not normalized!
            if recnet_model is not None:
                # Convert to tensor and apply RecNet directly to real space data
                segment_tensor = torch.tensor(segment_data, dtype=torch.float32).to(
                    device
                )
                refined = recnet_model(segment_tensor)
                segment_data = (
                    refined.cpu().numpy()
                )  # (n_timepoints, genes) - REAL SPACE

            all_segments.append(segment_data)

    # Concatenate all segments
    trajectory_data = np.vstack(all_segments)

    # Create DataFrame
    trajectory_df = pd.DataFrame(trajectory_data, columns=gene_data.index)

    # Create informative index
    time_indices = []
    for i in range(len(trajectory) - 1):
        for t in range(n_timepoints):
            time_indices.append(f"{trajectory[i]}_to_{trajectory[i + 1]}_t{t:03d}")
    trajectory_df.index = time_indices

    logger.info(f"Generated trajectory data: {trajectory_df.shape}")

    # Save if path provided
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        trajectory_df.to_csv(save_path)
        logger.info(f"Saved trajectory to: {save_path}")

    return trajectory_df
