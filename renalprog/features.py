"""
Feature engineering and preprocessing functionality for renalprog.

Includes functions for:
- Low expression filtering
- Mahalanobis outlier detection
- Gene clustering with tsfresh
- Feature extraction
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import chi2
from sklearn.covariance import MinCovDet
from typing import Tuple, Optional, List
import logging

from renalprog.config import PreprocessingConfig

logger = logging.getLogger(__name__)


def filter_low_expression(
    data: pd.DataFrame,
    mean_threshold: float = 0.5,
    var_threshold: float = 0.5,
    min_sample_fraction: float = 0.2
) -> pd.DataFrame:
    """
    Filter out genes with low expression across samples.
    
    Args:
        data: DataFrame with genes as rows and samples as columns
        mean_threshold: Minimum expression value to consider gene as expressed
        var_threshold: Minimum variance value to consider gene as expressed
        min_sample_fraction: Minimum fraction of non-expressed samples to filter gene
        
    Returns:
        Filtered DataFrame with only genes meeting expression criteria
    """
    logger.info(f"Starting with {data.shape[0]} genes")
    logger.info(f"Initial data shape: {data.shape}")
    # Remove lowly expressed genes
    filtered_data = data[(data == 0).sum(axis=1) / data.shape[1] <= min_sample_fraction]

    filtered_data= filtered_data.iloc[
        (np.mean(filtered_data, axis=1).values >= mean_threshold) &
        (np.var(filtered_data, axis=1).values >= var_threshold)]

    n_removed = data.shape[0] - filtered_data.shape[0]
    logger.info(f"Removed {n_removed} lowly expressed genes")
    logger.info(f"Retained {filtered_data.shape[0]} genes")
    
    return filtered_data


def detect_outliers_mahalanobis(
    data: pd.DataFrame,
    alpha: float = 0.05,
    support_fraction: float = None,
    transpose: bool = False,
    seed: int = 2023
) -> Tuple[pd.DataFrame, List[str], np.ndarray]:
    """
    Detect and remove outlier samples using Mahalanobis distance.
    
    Uses Minimum Covariance Determinant (MCD) for robust covariance estimation,
    then identifies outliers based on Mahalanobis distance and chi-square distribution.
    
    Args:
        data: DataFrame with samples as columns and genes as rows (will be transposed)
        alpha: Significance level for chi-square test (default: 0.05)
        support_fraction: Fraction of samples to use in MCD estimation (default: 0.75)
        transpose: Whether to transpose data before processing (default: True)
        seed: Random seed for reproducibility (default: 2023)
        
    Returns:
        Tuple of:
        - cleaned_data: DataFrame with outliers removed
        - outlier_ids: List of outlier sample IDs
        - mahalanobis_distances: Array of Mahalanobis distances for all samples
    """
    logger.info(f"Detecting outliers with Mahalanobis distance (alpha={alpha})")
    
    # Transpose if needed (we need samples as rows)
    if transpose:
        data_for_mcd = data.T
    else:
        data_for_mcd = data.copy()
    logger.info(f"Data shape for MCD: {data_for_mcd.shape} (features x samples)")
    # Define chi-square cutoff
    n_features = data_for_mcd.shape[0] # this is number of features (genes)
    cutoff = chi2.ppf(1 - alpha, n_features - 1)
    logger.info(f"Chi-square cutoff (df={n_features}, alpha={alpha}): {cutoff:.2f}")
    
    # Minimum Covariance Determinant for robust covariance estimation
    logger.info("Fitting MinCovDet estimator...")
    mcd = MinCovDet(support_fraction=support_fraction,random_state = seed)
    mcd.fit(data_for_mcd) # fit with (n_samples, n_features)
    
    # Calculate Mahalanobis distances
    mahalanobis_distances = mcd.dist_
    
    # Identify outliers
    outlier_mask = mahalanobis_distances > cutoff
    outlier_indices = np.where(outlier_mask)[0]
    outlier_ids = data_for_mcd.index[outlier_indices].tolist()
    
    logger.info(f"Detected {len(outlier_ids)} outlier samples")
    logger.info(f"Outlier IDs: {outlier_ids[:10]}{'...' if len(outlier_ids) > 10 else ''}")
    
    # Remove outliers
    if transpose:
        # Remove columns (samples) from original data
        cleaned_data = data.drop(columns=outlier_ids)
    else:
        # Remove rows from original data
        cleaned_data = data.drop(index=outlier_ids)
    
    logger.info(f"Cleaned data shape: {cleaned_data.shape}")
    
    return cleaned_data, outlier_ids, mahalanobis_distances


def preprocess_rnaseq(
    data: pd.DataFrame,
    filter_expression: bool = True,
    detect_outliers: bool = True,
    log_transform: bool = False,
    **kwargs
) -> Tuple[pd.DataFrame, dict]:
    """
    Complete preprocessing pipeline for RNA-seq data.
    
    Args:
        data: DataFrame with genes as rows and samples as columns
        filter_expression: Whether to filter lowly expressed genes
        detect_outliers: Whether to detect and remove outliers
        log_transform: Whether to apply log transformation
        **kwargs: Additional arguments for filtering and outlier detection
        
    Returns:
        Tuple of:
        - preprocessed_data: Cleaned DataFrame
        - info: Dictionary with preprocessing information
    """
    info = {
        "original_shape": str(data.shape),
        "original_genes": data.shape[0],
        "original_samples": data.shape[1],
        "steps_applied": []
    }
    
    processed_data = data.copy()
    
    # Filter low expression genes
    if filter_expression:
        mean_threshold = kwargs.get("mean_threshold", PreprocessingConfig.MEAN_EXPRESSION_THRESHOLD)
        var_threshold = kwargs.get("var_threshold", PreprocessingConfig.VAR_EXPRESSION_THRESHOLD)
        min_sample_fraction = kwargs.get("min_sample_fraction", PreprocessingConfig.MIN_SAMPLE_FRACTION)
        
        processed_data = filter_low_expression(
            data=processed_data,
            mean_threshold=mean_threshold,
            var_threshold=var_threshold,
            min_sample_fraction=min_sample_fraction
        )
        info["steps_applied"].append("filter_low_expression")
        info["genes_after_filtering"] = processed_data.shape[0]
        info["genes_removed"] = data.shape[0] - processed_data.shape[0]

    # Log transformation
    if log_transform:
        processed_data = np.log1p(processed_data)
        info["steps_applied"].append("log_transform")
    
    # Detect and remove outliers
    if detect_outliers:
        alpha = kwargs.get("alpha", 0.05)
        support_fraction = kwargs.get("support_fraction", None)
        transpose = kwargs.get("transpose", False)
        seed = kwargs.get("seed", 2023)
        processed_data, outlier_ids, maha_distances = detect_outliers_mahalanobis(
            processed_data,
            alpha=alpha,
            support_fraction=support_fraction,
            transpose=transpose,
            seed=seed
        )
        info["steps_applied"].append("detect_outliers_mahalanobis")
        info["outliers_removed"] = len(outlier_ids)
        info["outlier_ids"] = outlier_ids
        info["samples_after_outlier_removal"] = processed_data.shape[1]
    
    info["final_shape"] = str(processed_data.shape)
    info["final_genes"] = processed_data.shape[0]
    info["final_samples"] = processed_data.shape[1]
    info["steps_applied"] = ", ".join(info["steps_applied"])

    logger.info(f"Preprocessing complete: {info['original_shape']} -> {info['final_shape']}")
    
    return processed_data, info


def cluster_genes_tsfresh(
    trajectory_data: pd.DataFrame,
    n_clusters: int = 10,
    feature_extraction_settings: Optional[dict] = None
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Cluster genes based on their temporal patterns using tsfresh features.
    
    This function extracts time series features from gene expression trajectories
    and clusters genes with similar temporal patterns.
    
    Args:
        trajectory_data: DataFrame with trajectory data (time x genes)
        n_clusters: Number of clusters to create
        feature_extraction_settings: Optional tsfresh feature extraction settings
        
    Returns:
        Tuple of:
        - features_df: DataFrame with extracted tsfresh features
        - cluster_labels: Array of cluster assignments for each gene
    """
    try:
        from tsfresh import extract_features
        from tsfresh.utilities.dataframe_functions import impute
        from sklearn.cluster import KMeans
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        logger.error("tsfresh is required for gene clustering. Install with: pip install tsfresh")
        raise
    
    logger.info(f"Clustering genes with tsfresh (n_clusters={n_clusters})")
    
    # Prepare data for tsfresh
    # tsfresh expects: id, time, value columns for each time series
    formatted_data = []
    
    for gene in trajectory_data.columns:
        gene_data = pd.DataFrame({
            'id': gene,
            'time': range(len(trajectory_data)),
            'value': trajectory_data[gene].values
        })
        formatted_data.append(gene_data)
    
    tsfresh_data = pd.concat(formatted_data, ignore_index=True)
    
    # Extract features
    logger.info("Extracting tsfresh features...")
    features = extract_features(
        tsfresh_data,
        column_id='id',
        column_sort='time',
        column_value='value',
        default_fc_parameters=feature_extraction_settings
    )
    
    # Impute missing values
    impute(features)
    
    # Standardize features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)
    
    # Cluster genes
    logger.info("Clustering genes...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=2023, n_init=10)
    cluster_labels = kmeans.fit_predict(features_scaled)
    
    logger.info(f"Clustered {len(features)} genes into {n_clusters} clusters")
    logger.info(f"Cluster sizes: {pd.Series(cluster_labels).value_counts().sort_index().to_dict()}")
    
    return features, cluster_labels


def save_preprocessing_results(
    data: pd.DataFrame,
    metadata: dict,
    output_dir: Path,
    prefix: str = "rnaseq"
):
    """
    Save preprocessed data and metadata to files.
    
    Args:
        data: Preprocessed DataFrame
        metadata: Dictionary with preprocessing metadata
        output_dir: Directory to save files
        prefix: Prefix for output filenames
    """
    import json
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save data
    data_path = output_dir / f"{prefix}_preprocessed.csv"
    data.to_csv(data_path)
    logger.info(f"Saved preprocessed data to {data_path}")
    
    # Save metadata
    metadata_path = output_dir / f"{prefix}_preprocessing_metadata.json"
    
    # Convert non-serializable objects to serializable
    metadata_serializable = metadata.copy()
    if "outlier_ids" in metadata_serializable:
        metadata_serializable["outlier_ids"] = list(metadata_serializable["outlier_ids"])
    
    with open(metadata_path, 'w') as f:
        json.dump(metadata_serializable, f, indent=2)
    logger.info(f"Saved metadata to {metadata_path}")
