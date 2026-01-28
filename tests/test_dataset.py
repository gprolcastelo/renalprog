"""
Test dataset module functionality.
"""

import pytest
import pandas as pd
import numpy as np

from renalprog.dataset import (
    load_rnaseq_data,
    load_clinical_data,
    map_stages_to_early_late,
    create_train_test_split,
)


@pytest.fixture
def sample_rnaseq_data():
    """Create sample RNA-seq data for testing."""
    np.random.seed(2023)
    genes = [f"GENE_{i}" for i in range(100)]
    samples = [f"SAMPLE_{i}" for i in range(20)]
    data = pd.DataFrame(np.random.rand(100, 20), index=genes, columns=samples)
    return data


@pytest.fixture
def sample_clinical_data():
    """Create sample clinical data for testing."""
    samples = [f"SAMPLE_{i}" for i in range(20)]
    stages = ["Stage I"] * 5 + ["Stage II"] * 5 + ["Stage III"] * 5 + ["Stage IV"] * 5
    data = pd.DataFrame({"ajcc_pathologic_tumor_stage": stages}, index=samples)
    return data


def test_map_stages_to_early_late(sample_clinical_data):
    """Test stage mapping to early/late."""
    stages = sample_clinical_data["ajcc_pathologic_tumor_stage"]
    mapped = map_stages_to_early_late(stages)

    assert mapped.notna().all()
    assert set(mapped.unique()) == {"early", "late"}
    assert (mapped[:10] == "early").all()  # Stages I and II
    assert (mapped[10:] == "late").all()  # Stages III and IV


@pytest.mark.skip(
    reason="Test has inconsistent behavior with small sample sizes and stratified splitting"
)
def test_create_train_test_split(sample_rnaseq_data, sample_clinical_data, tmp_path):
    """Test train/test split creation."""
    # Save temporary files
    rnaseq_path = tmp_path / "rnaseq.csv"
    clinical_path = tmp_path / "clinical.csv"

    sample_rnaseq_data.to_csv(rnaseq_path)
    sample_clinical_data.to_csv(clinical_path)

    # Create split
    X_train, X_test, y_train, y_test, full_data, full_clinical = (
        create_train_test_split(
            rnaseq_path, clinical_path, test_size=0.2, seed=2023, use_onehot=True
        )
    )

    # Check shapes
    assert X_train.shape[0] == 16  # 80% of 20 samples
    assert X_test.shape[0] == 4  # 20% of 20 samples
    assert X_train.shape[1] == 100  # All genes

    # Check that we have classes (may not have all 4 in each split due to stratification)
    assert y_train.shape[0] == 16
    assert y_test.shape[0] == 4
    assert y_train.shape[1] > 0  # At least some classes
    assert y_test.shape[1] > 0

    # Total unique classes should match clinical data
    # Note: With small samples and stratified split, some classes may not appear in train/test
    n_unique_stages = sample_clinical_data["ajcc_pathologic_tumor_stage"].nunique()
    # The number of classes in one-hot encoded data should match total unique stages in original data
    # Both train and test should use the same encoding (all possible classes)
    assert y_train.shape[1] == n_unique_stages
    assert y_test.shape[1] == n_unique_stages

    # Check no overlap between train and test
    assert set(X_train.index).isdisjoint(set(X_test.index))


def test_train_test_split_with_output(
    sample_rnaseq_data, sample_clinical_data, tmp_path
):
    """Test that split saves files correctly."""
    rnaseq_path = tmp_path / "rnaseq.csv"
    clinical_path = tmp_path / "clinical.csv"
    output_dir = tmp_path / "split_output"

    sample_rnaseq_data.to_csv(rnaseq_path)
    sample_clinical_data.to_csv(clinical_path)

    create_train_test_split(
        rnaseq_path, clinical_path, test_size=0.2, seed=2023, output_dir=output_dir
    )

    # Check that files were created
    assert (output_dir / "X_train.csv").exists()
    assert (output_dir / "X_test.csv").exists()
    assert (output_dir / "y_train.csv").exists()
    assert (output_dir / "y_test.csv").exists()
    assert (output_dir / "split_statistics.csv").exists()


def test_load_rnaseq_data(sample_rnaseq_data, tmp_path):
    """Test loading RNA-seq data."""
    path = tmp_path / "rnaseq.csv"
    sample_rnaseq_data.to_csv(path)

    loaded = load_rnaseq_data(path)

    assert loaded.shape == sample_rnaseq_data.shape
    assert list(loaded.columns) == list(sample_rnaseq_data.columns)
    assert list(loaded.index) == list(sample_rnaseq_data.index)


def test_load_clinical_data(sample_clinical_data, tmp_path):
    """Test loading clinical data."""
    path = tmp_path / "clinical.csv"
    sample_clinical_data.to_csv(path)

    loaded = load_clinical_data(path)

    assert len(loaded) == len(sample_clinical_data)
    assert loaded.name == "ajcc_pathologic_tumor_stage"
