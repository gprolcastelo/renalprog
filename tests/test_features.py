"""
Test features module functionality.
"""

import pytest
import pandas as pd
import numpy as np

from renalprog.features import (
    filter_low_expression,
    detect_outliers_mahalanobis,
    preprocess_rnaseq,
)


@pytest.fixture
def sample_expression_data():
    """Create sample gene expression data with some low-expression genes."""
    np.random.seed(2023)

    # Create data with mix of expressed and non-expressed genes
    n_genes = 50
    n_samples = 50  # Increased from 30 to avoid MCD issues

    genes = [f"GENE_{i}" for i in range(n_genes)]
    samples = [f"SAMPLE_{i}" for i in range(n_samples)]

    # Half genes are well-expressed, half are lowly expressed
    data = np.random.rand(n_genes, n_samples)
    data[:25, :] = data[:25, :] * 10 + 5  # Well expressed genes
    data[25:, :] = data[25:, :] * 0.1  # Lowly expressed genes

    df = pd.DataFrame(data, index=genes, columns=samples)
    return df


def test_filter_low_expression(sample_expression_data):
    """Test filtering of lowly expressed genes."""
    filtered = filter_low_expression(
        sample_expression_data,
        mean_threshold=0.5,
        var_threshold=0.5,
        min_sample_fraction=0.5,  # Gene must be expressed in 50% of samples
    )

    # Should filter out roughly half the genes (the lowly expressed ones)
    assert filtered.shape[0] < sample_expression_data.shape[0]
    assert filtered.shape[1] == sample_expression_data.shape[1]
    assert filtered.shape[0] >= 20  # Should keep most of the well-expressed genes


def test_detect_outliers_mahalanobis(sample_expression_data):
    """Test Mahalanobis outlier detection."""
    # Add some outliers
    data_with_outliers = sample_expression_data.copy()
    data_with_outliers.iloc[:, -2:] = (
        data_with_outliers.iloc[:, -2:] * 100
    )  # Make last 2 samples outliers

    cleaned, outlier_ids, distances = detect_outliers_mahalanobis(
        data_with_outliers, alpha=0.05, support_fraction=0.75, transpose=True
    )

    # Should remove some outliers
    assert cleaned.shape[1] < data_with_outliers.shape[1]
    assert len(outlier_ids) > 0
    assert len(distances) == data_with_outliers.shape[1]

    # Check that outliers were actually removed
    for outlier_id in outlier_ids:
        assert outlier_id not in cleaned.columns


def test_preprocess_rnaseq(sample_expression_data):
    """Test complete preprocessing pipeline."""
    # Add outliers
    data_with_outliers = sample_expression_data.copy()
    data_with_outliers.iloc[:, -1] = data_with_outliers.iloc[:, -1] * 100

    processed, metadata = preprocess_rnaseq(
        data_with_outliers,
        filter_expression=True,
        detect_outliers=True,
        mean_threshold=0.5,
        var_threshold=0.5,
        min_sample_fraction=0.5,
        alpha=0.05,
        support_fraction=0.5,  # Use smaller support_fraction to avoid MCD issues
        transpose=True,  # Transpose for outlier detection on samples
    )

    # Check that data was processed
    assert processed.shape[0] <= data_with_outliers.shape[0]  # Genes may be filtered
    assert (
        processed.shape[1] <= data_with_outliers.shape[1]
    )  # Samples possibly filtered

    # Check metadata
    assert "original_shape" in metadata
    assert "final_shape" in metadata
    assert "steps_applied" in metadata
    assert len(metadata["steps_applied"]) > 0


def test_preprocess_no_outlier_detection(sample_expression_data):
    """Test preprocessing without outlier detection."""
    processed, metadata = preprocess_rnaseq(
        sample_expression_data,
        filter_expression=True,
        detect_outliers=False,
        mean_threshold=0.5,
        var_threshold=0.5,
    )

    # Should keep all samples
    assert processed.shape[1] == sample_expression_data.shape[1]
    assert "detect_outliers_mahalanobis" not in metadata["steps_applied"]


def test_filter_low_expression_no_filtering():
    """Test that highly expressed genes are not filtered."""
    # All genes well expressed
    data = pd.DataFrame(
        np.random.rand(10, 20) * 10 + 5,
        index=[f"GENE_{i}" for i in range(10)],
        columns=[f"SAMPLE_{i}" for i in range(20)],
    )

    filtered = filter_low_expression(
        data, mean_threshold=0.5, var_threshold=0.5, min_sample_fraction=0.1
    )

    # Should keep all genes
    assert filtered.shape == data.shape
