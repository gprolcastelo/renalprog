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

## Trajectory Generation

### generate_trajectories

Generate disease progression trajectories by interpolating in latent space.

::: renalprog.modeling.predict.generate_trajectories

**Example Usage:**

```python
from renalprog.modeling.predict import generate_trajectories

# Generate trajectories from early to late stage
trajectories = generate_trajectories(
    model=model,
    start_data=early_stage_samples.values,
    end_data=late_stage_samples.values,
    n_steps=50,
    interpolation='spherical',
    device='cuda'
)

# trajectories shape: (n_samples, n_steps, n_genes)
print(f"Generated {trajectories.shape[0]} trajectories")
print(f"Each with {trajectories.shape[1]} steps")
```

### create_patient_connections

Create optimal patient pairings for trajectory generation.

::: renalprog.modeling.predict.create_patient_connections

**Example Usage:**

```python
from renalprog.modeling.predict import create_patient_connections

# Find optimal patient connections
connections = create_patient_connections(
    latent_early=early_latent,
    latent_late=late_latent,
    method='closest',  # or 'random'
    output_path=Path("data/processed/patient_connections.csv")
)

print(f"Created {len(connections)} patient pairs")
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
from renalprog.modeling.train import VAE
from renalprog.modeling.predict import (
    apply_vae,
    create_patient_connections,
    generate_trajectories,
    evaluate_reconstruction
)
from renalprog.plots import plot_latent_space, plot_trajectory

# Load model and data
model = VAE(input_dim=20000, mid_dim=1024, features=128)
model.load_state_dict(torch.load("models/my_vae/best_model.pt"))

train_expr = pd.read_csv("data/interim/split/train_expression.tsv", sep="\t", index_col=0)
test_expr = pd.read_csv("data/interim/split/test_expression.tsv", sep="\t", index_col=0)
clinical = pd.read_csv("data/interim/split/test_clinical.tsv", sep="\t", index_col=0)

# Encode data
train_results = apply_vae(model, train_expr.values, device='cuda')
test_results = apply_vae(model, test_expr.values, device='cuda')

# Visualize latent space
plot_latent_space(
    latent=test_results['latent'],
    labels=clinical['stage'],
    output_path=Path("reports/figures/latent_space.png")
)

# Create patient connections
early_mask = clinical['stage'] == 'early'
late_mask = clinical['stage'] == 'late'

connections = create_patient_connections(
    latent_early=test_results['latent'][early_mask],
    latent_late=test_results['latent'][late_mask],
    method='closest',
    output_path=Path("data/processed/connections.csv")
)

# Generate trajectories
trajectories = generate_trajectories(
    model=model,
    start_data=test_expr.values[early_mask],
    end_data=test_expr.values[late_mask],
    n_steps=50,
    interpolation='spherical',
    device='cuda'
)

# Evaluate reconstruction
metrics = evaluate_reconstruction(
    model=model,
    original_data=test_expr.values,
    device='cuda',
    output_dir=Path("reports/reconstruction")
)

print(f"Generated {len(trajectories)} trajectories")
print(f"Reconstruction MSE: {metrics['mse']:.4f}")
```

## See Also

- [Training API](training.md) - Train VAE models
- [Trajectories API](trajectories.md) - Trajectory analysis
- [Models API](models.md) - VAE architectures

