# Prediction API

Functions for applying trained VAE models to generate latent representations and trajectories.

## Overview

This module provides:

- Apply trained VAE to encode data
- Generate disease progression trajectories
- Evaluate reconstruction quality
- Patient connectivity analysis
- Latent space interpolation

## Core Prediction Functions

### apply_vae

Apply trained VAE model to encode gene expression data into latent space.

::: renalprog.modeling.predict.apply_vae

**Example Usage:**

```python
import torch
import pandas as pd
from pathlib import Path
from renalprog.modeling.train import VAE
from renalprog.modeling.predict import apply_vae

# Load model
model = VAE(input_dim=20000, mid_dim=1024, features=128)
model.load_state_dict(torch.load("models/my_vae/best_model.pt"))

# Load test data
test_expr = pd.read_csv("data/interim/split/test_expression.tsv", sep="\t", index_col=0)

# Apply VAE
results = apply_vae(
    model=model,
    data=test_expr.values,
    device='cuda',
    batch_size=32
)

latent = results['latent']  # Latent representations
reconstructed = results['reconstructed']  # Reconstructed expression
print(f"Latent space shape: {latent.shape}")
```



### interpolate_latent_linear

Linear interpolation between latent representations.

::: renalprog.modeling.predict.interpolate_latent_linear

### interpolate_latent_spherical

Spherical (SLERP) interpolation between latent representations.

::: renalprog.modeling.predict.interpolate_latent_spherical

**Example:**

```python
from renalprog.modeling.predict import interpolate_latent_spherical
import numpy as np

z_start = np.random.randn(10, 128)  # 10 samples, 128 latent dims
z_end = np.random.randn(10, 128)

# Spherical interpolation (better for normalized spaces)
trajectory = interpolate_latent_spherical(z_start, z_end, n_steps=50)
# Shape: (10, 50, 128)
```

## Reconstruction Evaluation

### evaluate_reconstruction

Comprehensive evaluation of VAE reconstruction quality.

::: renalprog.modeling.predict.evaluate_reconstruction

**Example Usage:**

```python
from renalprog.modeling.predict import evaluate_reconstruction

# Evaluate reconstruction quality
metrics = evaluate_reconstruction(
    model=model,
    original_data=test_expr.values,
    device='cuda',
    output_dir=Path("reports/reconstruction_eval")
)

print(f"MSE: {metrics['mse']:.4f}")
print(f"Pearson R: {metrics['pearson_mean']:.4f}")
print(f"Cosine similarity: {metrics['cosine_mean']:.4f}")
```

## Quality Metrics

### diagnostic_metrics

Calculate diagnostic metrics for model evaluation.

::: renalprog.modeling.predict.diagnostic_metrics

### quality_metrics

Calculate quality metrics for generated trajectories.

::: renalprog.modeling.predict.quality_metrics

## Trajectory Classification

### classify_trajectories

Classify disease progression trajectories as progressing vs. non-progressing.

::: renalprog.modeling.predict.classify_trajectories

**Example Usage:**

```python
from renalprog.modeling.predict import classify_trajectories

# Train classifier on trajectories
classifier, metrics = classify_trajectories(
    trajectories=trajectory_data,
    labels=progression_labels,
    output_dir=Path("models/trajectory_classifier")
)

print(f"Classification accuracy: {metrics['accuracy']:.3f}")
print(f"AUC-ROC: {metrics['auc_roc']:.3f}")
```

## Network Analysis

### build_trajectory_network

Build network graph of patient trajectories.

::: renalprog.modeling.predict.build_trajectory_network

### link_patients_closest

Link patients using closest neighbor strategy.

::: renalprog.modeling.predict.link_patients_closest

### link_patients_random

Link patients randomly (for control/comparison).

::: renalprog.modeling.predict.link_patients_random

## Dynamic Analysis

### dynamic_enrichment_analysis

Perform pathway enrichment along trajectory timepoints.

::: renalprog.modeling.predict.dynamic_enrichment_analysis

### calculate_all_possible_transitions

Calculate all possible patient-to-patient transitions.

::: renalprog.modeling.predict.calculate_all_possible_transitions

## Metadata

### get_metadata

Extract metadata from model directory.

::: renalprog.modeling.predict.get_metadata

## Complete Example

```python
import torch
import pandas as pd
from pathlib import Path
from sklearn.preprocessing import MinMaxScaler
from renalprog.modeling.train import VAE, NetworkReconstruction
from renalprog.modeling.predict import (
    apply_vae,
    generate_trajectory_data,
    evaluate_reconstruction
)
from renalprog.plots import plot_trajectory

# Load model and data
vae_model = VAE(input_dim=20000, mid_dim=1024, features=128)
vae_model.load_state_dict(torch.load("models/my_vae/best_model.pt"))

recnet_model = NetworkReconstruction([20000, 5000, 20000])
recnet_model.load_state_dict(torch.load("models/my_vae/recnet_model.pt"))

train_expr = pd.read_csv("data/interim/split/train_expression.tsv", sep="\t", index_col=0)
test_expr = pd.read_csv("data/interim/split/test_expression.tsv", sep="\t", index_col=0)
clinical = pd.read_csv("data/interim/split/test_clinical.tsv", sep="\t", index_col=0)

# Encode data
train_results = apply_vae(vae_model, train_expr.values, device='cuda')
test_results = apply_vae(vae_model, test_expr.values, device='cuda')

# Prepare scaler (use same as training)
scaler = MinMaxScaler()
scaler.fit(train_expr.values.T)  # Fit on training data

# Define a trajectory (list of patient IDs)
trajectory_patients = ['TCGA-A1-001', 'TCGA-A2-002', 'TCGA-A3-003']

# Generate trajectory data
trajectory_data = generate_trajectory_data(
    vae_model=vae_model,
    recnet_model=recnet_model,
    trajectory=trajectory_patients,
    gene_data=test_expr.T,  # Transpose: genes as rows, patients as columns
    n_timepoints=50,
    interpolation_method='spherical',
    device='cuda',
    scaler=scaler
)

# Evaluate reconstruction
metrics = evaluate_reconstruction(
    model=vae_model,
    original_data=test_expr.values,
    device='cuda',
    output_dir=Path("reports/reconstruction")
)

print(f"Generated trajectory with {len(trajectory_data)} time points")
print(f"Reconstruction MSE: {metrics['mse']:.4f}")
```

## See Also

- [Training API](training.md) - Train VAE models
- [Trajectories API](trajectories.md) - Trajectory analysis
- [Models API](models.md) - VAE architectures

