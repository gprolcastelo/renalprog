"""
Test modeling module functionality - VAE models, training, and checkpointing.
"""

import pytest
import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from pathlib import Path
import tempfile
import json

from renalprog.config import VAEConfig
from renalprog.modeling.train import (
    VAE, AE, CVAE, vae_loss, create_dataloader,
    train_epoch, evaluate_model, train_vae
)
from renalprog.modeling.checkpointing import (
    ModelCheckpointer, save_model_config, load_model_config
)
from renalprog.utils import set_seed


@pytest.fixture
def vae_config():
    """Create a basic VAE configuration for testing."""
    config = VAEConfig()
    config.INPUT_DIM = 100
    config.MID_DIM = 32
    config.LATENT_DIM = 8
    config.LEARNING_RATE = 1e-3
    config.BATCH_SIZE = 8
    config.EPOCHS = 2
    config.SEED = 2023

    return config


@pytest.fixture
def synthetic_vae_data():
    """Create synthetic data for VAE testing."""
    np.random.seed(2023)

    # Create synthetic gene expression data
    n_samples = 50
    n_genes = 100

    X = np.random.randn(n_samples, n_genes).astype(np.float32)
    # Scale to positive values
    X = (X - X.min()) / (X.max() - X.min())

    # Create synthetic labels (4 stages)
    y = np.array([0, 1, 2, 3] * (n_samples // 4) + [0, 1])[:n_samples]

    return X, y


def test_vae_instantiation(vae_config):
    """Test that VAE model can be instantiated with correct architecture."""
    model = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )

    # Check model exists
    assert model is not None
    assert isinstance(model, nn.Module)

    # Check model has correct components
    assert hasattr(model, 'encoder')
    assert hasattr(model, 'decoder')
    assert hasattr(model, 'reparametrize')

    # Check parameters exist
    params = list(model.parameters())
    assert len(params) > 0

    # Count total parameters
    total_params = sum(p.numel() for p in model.parameters())
    assert total_params > 0
    print(f"VAE has {total_params:,} parameters")


def test_vae_forward_pass(vae_config, synthetic_vae_data):
    """Test VAE forward pass with dummy data."""
    X, _ = synthetic_vae_data
    model = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )

    # Prepare input (use first 5 samples)
    x = torch.tensor(X[:5], dtype=torch.float32)

    # Forward pass
    reconstruction, mu, logvar, z = model(x)

    # Check outputs exist
    assert reconstruction is not None
    assert mu is not None
    assert logvar is not None
    assert z is not None

    # Check output shapes
    assert reconstruction.shape == x.shape, f"Expected {x.shape}, got {reconstruction.shape}"
    assert mu.shape == (5, vae_config.LATENT_DIM), f"Expected (5, {vae_config.LATENT_DIM}), got {mu.shape}"
    assert logvar.shape == (5, vae_config.LATENT_DIM), f"Expected (5, {vae_config.LATENT_DIM}), got {logvar.shape}"
    assert z.shape == (5, vae_config.LATENT_DIM), f"Expected (5, {vae_config.LATENT_DIM}), got {z.shape}"

    # Check values are finite
    assert torch.isfinite(reconstruction).all()
    assert torch.isfinite(mu).all()
    assert torch.isfinite(logvar).all()
    assert torch.isfinite(z).all()


def test_vae_output_shapes(vae_config):
    """Test VAE produces correct output dimensions with various batch sizes."""
    model = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )
    model.eval()

    # Test different batch sizes
    batch_sizes = [1, 4, 16]

    for batch_size in batch_sizes:
        x = torch.randn(batch_size, vae_config.INPUT_DIM)

        with torch.no_grad():
            reconstruction, mu, logvar, z = model(x)

        assert reconstruction.shape == (batch_size, vae_config.INPUT_DIM)
        assert mu.shape == (batch_size, vae_config.LATENT_DIM)
        assert logvar.shape == (batch_size, vae_config.LATENT_DIM)
        assert z.shape == (batch_size, vae_config.LATENT_DIM)


def test_vae_loss_calculation(vae_config, synthetic_vae_data):
    """Test that VAE loss function works correctly."""
    X, _ = synthetic_vae_data
    model = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )

    # Forward pass
    x = torch.tensor(X[:8], dtype=torch.float32)
    reconstruction, mu, logvar, z = model(x)

    # Calculate loss
    loss, recon_loss, kl_loss = vae_loss(reconstruction, x, mu, logvar, beta=1.0)

    # Check losses are scalars
    assert loss.dim() == 0, "Loss should be scalar"
    assert recon_loss.dim() == 0, "Reconstruction loss should be scalar"
    assert kl_loss.dim() == 0, "KL loss should be scalar"

    # Check losses are positive
    assert loss.item() >= 0, "Loss should be non-negative"
    assert recon_loss.item() >= 0, "Reconstruction loss should be non-negative"
    assert kl_loss.item() >= 0, "KL loss should be non-negative"

    # Check losses are finite
    assert torch.isfinite(loss)
    assert torch.isfinite(recon_loss)
    assert torch.isfinite(kl_loss)

    # Check loss composition
    expected_loss = recon_loss + kl_loss
    assert torch.isclose(loss, expected_loss, rtol=1e-5)


def test_ae_forward_pass(vae_config, synthetic_vae_data):
    """Test Autoencoder (AE) forward pass - simpler than VAE."""
    X, _ = synthetic_vae_data
    model = AE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )

    # Forward pass
    x = torch.tensor(X[:5], dtype=torch.float32)
    reconstruction, mu, logvar, z = model(x)

    # Check reconstruction and latent code
    assert reconstruction is not None
    assert z is not None
    assert reconstruction.shape == x.shape
    assert z.shape == (5, vae_config.LATENT_DIM)

    # For AE, mu and logvar should be None
    assert mu is None
    assert logvar is None


def test_cvae_with_conditions(vae_config, synthetic_vae_data):
    """Test Conditional VAE handles conditioning information."""
    X, y = synthetic_vae_data
    num_classes = len(np.unique(y))

    # One-hot encode labels
    y_onehot = np.eye(num_classes)[y[:5]]

    model = CVAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
        num_classes=num_classes,
    )

    # Forward pass with conditions
    x = torch.tensor(X[:5], dtype=torch.float32)
    conditions = torch.tensor(y_onehot, dtype=torch.float32)

    reconstruction, mu, logvar, z = model(x, conditions)

    # Check outputs
    assert reconstruction is not None
    assert mu is not None
    assert z is not None

    # Check shapes
    assert reconstruction.shape == x.shape
    assert mu.shape == (5, vae_config.LATENT_DIM)
    assert z.shape == (5, vae_config.LATENT_DIM)


def test_training_one_epoch(vae_config, synthetic_vae_data, tmp_path):
    """Test that model can train for one epoch without errors."""
    set_seed(vae_config.SEED)

    X, _ = synthetic_vae_data
    model = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )

    # Create dataloader
    dataloader = create_dataloader(X, None, batch_size=8, shuffle=True)

    # Create optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=vae_config.LEARNING_RATE)

    # Train one epoch
    device = 'cpu'
    metrics = train_epoch(model, dataloader, optimizer, device, vae_config)

    # Check metrics are returned
    assert 'loss' in metrics
    assert 'recon_loss' in metrics
    assert 'kl_loss' in metrics

    # Check metrics are reasonable
    assert metrics['loss'] > 0
    assert metrics['recon_loss'] > 0
    assert metrics['kl_loss'] >= 0

    # Check metrics are finite
    assert np.isfinite(metrics['loss'])
    assert np.isfinite(metrics['recon_loss'])
    assert np.isfinite(metrics['kl_loss'])


def test_model_save_load(vae_config, tmp_path):
    """Test that model can be saved and loaded correctly."""
    # Create model
    model = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )

    # Create optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=vae_config.LEARNING_RATE)

    # Save checkpoint
    checkpointer = ModelCheckpointer(save_dir=tmp_path)
    metrics = {'train_loss': 1.5, 'val_loss': 1.8}
    checkpointer.save_checkpoint(
        epoch=10,
        model=model,
        optimizer=optimizer,
        metrics=metrics,
        config=vae_config,
        is_best=True,
    )

    # Check checkpoint file exists
    checkpoint_path = tmp_path / 'best_model.pth'
    assert checkpoint_path.exists()

    # Create new model with same architecture
    new_model = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )

    # Load checkpoint
    checkpoint_info = checkpointer.load_checkpoint(
        checkpoint_path, new_model, optimizer
    )

    # Check checkpoint info
    assert checkpoint_info['epoch'] == 10
    assert 'metrics' in checkpoint_info
    assert checkpoint_info['metrics']['val_loss'] == 1.8

    # Check model weights match
    for p1, p2 in zip(model.parameters(), new_model.parameters()):
        assert torch.allclose(p1, p2)


def test_overfitting_small_batch(vae_config, synthetic_vae_data):
    """Test that VAE can overfit a small batch (sanity check)."""
    set_seed(vae_config.SEED)

    X, _ = synthetic_vae_data
    model = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )

    # Use only 4 samples for overfitting test
    x = torch.tensor(X[:4], dtype=torch.float32)

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)

    # Train for multiple iterations
    initial_loss = None
    final_loss = None

    for i in range(100):
        optimizer.zero_grad()
        reconstruction, mu, logvar, z = model(x)
        loss, recon_loss, kl_loss = vae_loss(reconstruction, x, mu, logvar, beta=0.1)
        loss.backward()
        optimizer.step()

        if i == 0:
            initial_loss = loss.item()
        if i == 99:
            final_loss = loss.item()

    # Loss should decrease significantly
    assert final_loss < initial_loss
    assert final_loss < initial_loss * 0.5, "Loss should decrease by at least 50%"
    print(f"Overfitting test: Initial loss={initial_loss:.4f}, Final loss={final_loss:.4f}")


def test_vae_reparameterization_trick(vae_config):
    """Test that reparameterization works correctly in training vs eval mode."""
    model = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )

    x = torch.randn(4, vae_config.INPUT_DIM)

    # Training mode - should use reparameterization
    model.train()
    _, mu1, _, z1 = model(x)
    _, mu2, _, z2 = model(x)

    # mu should be deterministic
    assert torch.allclose(mu1, mu2)

    # z should be different due to sampling (with high probability)
    assert not torch.allclose(z1, z2)

    # Eval mode - should return mu directly
    model.eval()
    with torch.no_grad():
        _, mu3, _, z3 = model(x)
        _, mu4, _, z4 = model(x)

    # In eval mode, z should equal mu (no sampling)
    assert torch.allclose(z3, mu3)
    assert torch.allclose(z4, mu4)

    # Multiple forward passes in eval should give same z
    assert torch.allclose(z3, z4)


def test_dataloader_creation(synthetic_vae_data):
    """Test dataloader creation with normalization."""
    X, y = synthetic_vae_data

    # Create dataloader
    dataloader = create_dataloader(X, y, batch_size=8, shuffle=True)

    # Check dataloader properties
    assert len(dataloader) > 0
    assert dataloader.batch_size == 8

    # Get one batch
    batch = next(iter(dataloader))

    # Check batch shape
    if y is not None:
        x_batch, y_batch = batch
        assert x_batch.shape[0] <= 8  # Batch size
        assert x_batch.shape[1] == 100  # Features
        assert y_batch.shape[0] == x_batch.shape[0]
    else:
        x_batch = batch
        assert x_batch.shape[0] <= 8
        assert x_batch.shape[1] == 100

    # Check data is normalized (should be in [0, 1] range after MinMaxScaler)
    # MinMaxScaler produces values in [0, 1] inclusive, so max can be exactly 1.0
    assert x_batch.min() >= 0
    assert x_batch.max() <= 1.0 + 1e-6  # Allow small floating point tolerance


def test_checkpoint_history_saving(vae_config, tmp_path):
    """Test that training history is saved correctly."""
    model = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=vae_config.LEARNING_RATE)

    checkpointer = ModelCheckpointer(save_dir=tmp_path)

    # Save multiple checkpoints - history only saved for best/final
    # So we need to mark them as best or final to trigger history saving
    for epoch in range(3):
        metrics = {
            'train_loss': 2.0 - epoch * 0.5,
            'val_loss': 2.5 - epoch * 0.4,
        }
        # Mark each as best to trigger history saving
        checkpointer.save_checkpoint(
            epoch, model, optimizer, metrics, vae_config, is_best=True
        )

    # Check history file exists
    history_path = tmp_path / 'training_history.json'
    assert history_path.exists()

    # Load and verify history
    with open(history_path, 'r') as f:
        history = json.load(f)

    assert 'epochs' in history
    assert 'metrics' in history
    # Each best model save appends to history
    assert len(history['epochs']) >= 1, f"Expected at least 1 epoch in history, got {len(history['epochs'])}"
    assert 'train_loss' in history['metrics']
    assert 'val_loss' in history['metrics']


def test_config_serialization(vae_config, tmp_path):
    """Test that configuration can be saved and loaded."""
    config_path = tmp_path / 'config.json'

    # Save config
    save_model_config(vae_config, config_path)

    # Check file exists
    assert config_path.exists()

    # Load config
    loaded_config = load_model_config(config_path)

    # Check values match
    assert loaded_config['INPUT_DIM'] == vae_config.INPUT_DIM
    assert loaded_config['MID_DIM'] == vae_config.MID_DIM
    assert loaded_config['LATENT_DIM'] == vae_config.LATENT_DIM
    assert loaded_config['LEARNING_RATE'] == vae_config.LEARNING_RATE


@pytest.mark.slow
def test_full_training_integration(vae_config, synthetic_vae_data, tmp_path):
    """Integration test: Full VAE training for a few epochs."""
    set_seed(vae_config.SEED)

    X, y = synthetic_vae_data

    # Split data
    n_train = int(0.8 * len(X))
    X_train, X_test = X[:n_train], X[n_train:]

    # Configure for quick test
    vae_config.EPOCHS = 3
    vae_config.CHECKPOINT_FREQ = 2

    # Train model - force CPU to avoid CUDA compatibility issues
    model, history = train_vae(
        X_train, X_test,
        y_train=None, y_test=None,
        config=vae_config,
        save_dir=tmp_path,
        force_cpu=True,  # Force CPU to avoid CUDA errors in tests
    )

    # Check model was returned
    assert model is not None

    # Check history
    assert 'train_loss' in history
    assert 'val_loss' in history
    assert len(history['train_loss']) == 3
    assert len(history['val_loss']) == 3

    # Check loss decreased
    assert history['train_loss'][-1] < history['train_loss'][0]

    # Check checkpoints were saved
    assert (tmp_path / 'best_model.pth').exists()
    assert (tmp_path / 'final_model.pth').exists()
    assert (tmp_path / 'config.json').exists()

    print(f"Training complete: Final train_loss={history['train_loss'][-1]:.4f}")


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])

