# Models API

The `modeling` module provides neural network architectures and training functions for variational autoencoders (VAEs).

## Overview

This module includes:

- VAE architectures (standard, conditional, simple)
- Training and evaluation functions
- Loss functions (reconstruction, KL divergence)
- Checkpoint management
- Post-processing networks

## Model Architectures

### VAE

Standard Variational Autoencoder with encoder-decoder architecture.

::: renalprog.modeling.train.VAE

**Example Usage:**

```python
import torch
from renalprog.modeling.train import VAE

# Create VAE model
model = VAE(
    input_dim=20000,  # Number of genes
    mid_dim=1024,     # Hidden layer size
    features=128,     # Latent dimension
    dropout=0.1
)

# Forward pass
x = torch.randn(32, 20000)  # Batch of gene expression
reconstruction, mu, log_var, z = model(x)
```

### CVAE

Conditional VAE that incorporates clinical covariates.

::: renalprog.modeling.train.CVAE

**Example Usage:**

```python
from renalprog.modeling.train import CVAE

# Create conditional VAE
model = ConditionalVAE(
    input_dim=20000,
    mid_dim=1024,
    features=128,
    condition_dim=2,  # e.g., one-hot encoded stage
    dropout=0.1
)

# Forward pass with condition
x = torch.randn(32, 20000)
condition = torch.randn(32, 2)  # Clinical covariates
reconstruction, mu, log_var, z = model(x, condition)
```

### AE

Simplified autoencoder without variational component.

::: renalprog.modeling.train.AE

## Loss Functions

### vae_loss

Complete VAE loss combining reconstruction and KL divergence.

::: renalprog.modeling.train.vae_loss

### reconstruction_loss

MSE-based reconstruction loss.

::: renalprog.modeling.train.reconstruction_loss

### kl_divergence

KL divergence between latent distribution and prior.

::: renalprog.modeling.train.kl_divergence

## Training Functions

### train_vae

Main training function for VAE models.

::: renalprog.modeling.train.train_vae

**Example Usage:**

```python
from renalprog.modeling.train import train_vae
from pathlib import Path
import pandas as pd

# Load training data
train_expr = pd.read_csv("data/interim/split/train_expression.tsv", sep="\t", index_col=0)
test_expr = pd.read_csv("data/interim/split/test_expression.tsv", sep="\t", index_col=0)

# Train VAE
history, best_model, checkpoints = train_vae(
    train_data=train_expr.values,
    val_data=test_expr.values,
    input_dim=train_expr.shape[1],
    mid_dim=1024,
    features=128,
    output_dir=Path("models/my_vae"),
    n_epochs=100,
    batch_size=32,
    learning_rate=1e-3,
    use_scheduler=True,
    use_checkpoint=True,
    early_stopping_patience=20
)

print(f"Final validation loss: {history['val_loss'][-1]:.4f}")
```

### train_epoch

Train the model for one epoch.

::: renalprog.modeling.train.train_epoch

### evaluate_model

Evaluate model on validation/test data.

::: renalprog.modeling.train.evaluate_model

### train_vae_with_postprocessing

Train VAE and post-processing network together.

::: renalprog.modeling.train.train_vae_with_postprocessing

## Utility Functions

### create_dataloader

Create PyTorch DataLoader from numpy arrays.

::: renalprog.modeling.train.create_dataloader

### frange_cycle_linear

Generate cyclical annealing schedule for KL divergence.

::: renalprog.modeling.train.frange_cycle_linear

**Example Usage:**

```python
from renalprog.modeling.train import frange_cycle_linear

# Create annealing schedule
schedule = frange_cycle_linear(
    n_iter=1000,
    start=0.0,
    stop=1.0,
    n_cycle=4,
    ratio=0.5
)

# Use in training loop
for i, beta in enumerate(schedule):
    loss = reconstruction_loss + beta * kl_loss
```

## Post-Processing Network

### NetworkReconstruction

Neural network for refining VAE reconstructions.

::: renalprog.modeling.train.NetworkReconstruction

### train_reconstruction_network

Train post-processing network.

::: renalprog.modeling.train.train_reconstruction_network

## See Also

- [Training API](training.md) - Complete training pipeline
- [Prediction API](prediction.md) - Using trained models
- [Configuration](config.md) - Model hyperparameters

