"""
Shared test fixtures and utilities.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path


@pytest.fixture
def tiny_dataset():
    """
    Create a tiny dataset for quick testing (10 samples, 50 genes).
    """
    np.random.seed(2023)

    genes = [f"GENE_{i:04d}" for i in range(50)]
    samples = [f"TCGA-TEST-{i:02d}" for i in range(10)]

    # Create expression data with some structure
    data = pd.DataFrame(
        np.random.randn(50, 10) * 2 + 5,
        index=genes,
        columns=samples
    )

    # Create clinical data with stages
    clinical = pd.DataFrame({
        "ajcc_pathologic_tumor_stage": ["Stage I", "Stage I", "Stage II", "Stage II", "Stage III", "Stage III", "Stage IV", "Stage IV", "Stage I", "Stage IV"]
    }, index=samples)

    return data, clinical


@pytest.fixture
def small_dataset():
    """
    Create a small dataset for testing (50 samples, 200 genes).
    """
    np.random.seed(2023)

    genes = [f"GENE_{i:04d}" for i in range(200)]
    samples = [f"TCGA-TEST-{i:03d}" for i in range(50)]

    # Create expression data
    data = pd.DataFrame(
        np.random.randn(200, 50) * 2 + 5,
        index=genes,
        columns=samples
    )

    # Create clinical data with balanced stages
    stages = ["Stage I"] * 12 + ["Stage II"] * 13 + ["Stage III"] * 12 + ["Stage IV"] * 13
    clinical = pd.DataFrame({
        "ajcc_pathologic_tumor_stage": stages
    }, index=samples)

    return data, clinical


@pytest.fixture
def mock_vae_model():
    """
    Create a mock VAE model for testing.
    """
    import torch
    import torch.nn as nn

    class MockVAE(nn.Module):
        def __init__(self, input_dim=200, latent_dim=16):
            super().__init__()
            self.input_dim = input_dim
            self.latent_dim = latent_dim
            self.encoder = nn.Linear(input_dim, latent_dim * 2)
            self.decoder = nn.Linear(latent_dim, input_dim)

        def forward(self, x):
            encoded = self.encoder(x)
            mu = encoded[:, :self.latent_dim]
            logvar = encoded[:, self.latent_dim:]
            z = mu  # Simplified, no reparameterization for testing
            reconstruction = self.decoder(z)
            return reconstruction, mu, logvar, z

    return MockVAE()


@pytest.fixture
def temp_data_dir(tmp_path):
    """
    Create a temporary directory structure for data.
    """
    data_dir = tmp_path / "data"
    (data_dir / "raw").mkdir(parents=True)
    (data_dir / "interim").mkdir(parents=True)
    (data_dir / "processed").mkdir(parents=True)
    (data_dir / "external").mkdir(parents=True)

    return data_dir

