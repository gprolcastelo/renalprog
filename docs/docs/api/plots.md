# Plots API

Visualization functions for gene expression analysis, model training, and results presentation.

## Overview

The plots module provides publication-quality visualization for:

- Training history and loss curves
- Latent space representations
- Gene expression heatmaps
- Trajectories and pathways
- Confusion matrices
- Enrichment results

## Core Plotting Functions

### save_plot

Utility function for saving plots with consistent formatting.

::: renalprog.plots.save_plot

**Example Usage:**

```python
from renalprog.plots import save_plot
import matplotlib.pyplot as plt
from pathlib import Path

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot([1, 2, 3], [4, 5, 6])
ax.set_title("My Plot")

save_plot(
    fig=fig,
    output_path=Path("reports/figures/my_plot.png"),
    dpi=300,
    bbox_inches='tight'
)
```

## Training Visualization

### plot_training_history

Visualize VAE training progress.

::: renalprog.plots.plot_training_history

**Example:**

```python
from renalprog.plots import plot_training_history
from pathlib import Path

# After training
history, model, checkpoints = train_vae(...)

# Plot training curves
plot_training_history(
    history=history,
    output_path=Path("reports/figures/training_history.png"),
    title="VAE Training Progress"
)
```

### plot_reconstruction_losses

Compare reconstruction losses across samples.

::: renalprog.plots.plot_reconstruction_losses


### plot_umap_plotly

Interactive UMAP visualization.

::: renalprog.plots.plot_umap_plotly

**Example:**

```python
from renalprog.plots import plot_umap_plotly

# Create interactive plot
fig = plot_umap_plotly(
    latent=latent,
    labels=clinical['stage'],
    sample_names=clinical.index.tolist(),
    title="Interactive Latent Space"
)

# Save as HTML
fig.write_html("reports/figures/latent_space_interactive.html")
```

## Trajectory Visualization

### plot_trajectory

Visualize disease progression trajectory.

::: renalprog.plots.plot_trajectory

**Example:**

```python
from renalprog.plots import plot_trajectory

# Plot single trajectory
plot_trajectory(
    trajectory=trajectories[0],  # Shape: (n_steps, n_genes)
    feature_names=selected_genes,
    output_path=Path("reports/figures/trajectory_001.png"),
    title="Disease Progression Trajectory",
    highlight_genes=['TP53', 'VEGFA', 'HIF1A']
)
```

## PCA Visualization

### plot_pca_variance

Visualize PCA variance explained.

::: renalprog.plots.plot_pca_variance

**Example:**

```python
from renalprog.plots import plot_pca_variance
from sklearn.decomposition import PCA

# Perform PCA
pca = PCA(n_components=50)
pca.fit(expression_data)

# Plot variance explained
plot_pca_variance(
    pca=pca,
    output_path=Path("reports/figures/pca_variance.png"),
    n_components=20
)
```

## Complete Visualization Workflow

```python
import torch
import pandas as pd
from pathlib import Path
from renalprog.modeling.train import VAE, train_vae
from renalprog.modeling.predict import apply_vae, generate_trajectories
from renalprog.plots import (
    plot_training_history,
    plot_trajectory,
    plot_umap_plotly
)

# Create output directory
output_dir = Path("reports/figures")
output_dir.mkdir(parents=True, exist_ok=True)

# 1. Load data
train_expr = pd.read_csv("data/interim/split/train_expression.tsv", sep="\t", index_col=0)
test_expr = pd.read_csv("data/interim/split/test_expression.tsv", sep="\t", index_col=0)
clinical = pd.read_csv("data/interim/split/test_clinical.tsv", sep="\t", index_col=0)

# 2. Train model and plot history
history, model, checkpoints = train_vae(
    train_data=train_expr.values,
    val_data=test_expr.values,
    input_dim=train_expr.shape[1],
    mid_dim=1024,
    features=128,
    output_dir=Path("models/my_vae"),
    n_epochs=100
)

plot_training_history(
    history=history,
    output_path=output_dir / "training_history.png"
)

# 3. Encode to latent space and visualize
results = apply_vae(model, test_expr.values, device='cuda')

plot_umap_plotly(
    latent=results['latent'],
    labels=clinical['stage'],
    sample_names=clinical.index.tolist(),
    title="Interactive Latent Space"
).write_html(output_dir / "latent_space_interactive.html")


# 4. Generate and plot trajectories
early_mask = clinical['stage'] == 'early'
late_mask = clinical['stage'] == 'late'

trajectories = generate_trajectories(
    model=model,
    start_data=test_expr.values[early_mask],
    end_data=test_expr.values[late_mask],
    n_steps=50,
    device='cuda'
)

# Plot first trajectory
plot_trajectory(
    trajectory=trajectories[0],
    feature_names=top_genes.tolist(),
    output_path=output_dir / "trajectory_001.png",
    title="Disease Progression Trajectory"
)

print(f"All figures saved to {output_dir}")
```

## See Also

- [Training API](training.md) - Generate training history
- [Prediction API](prediction.md) - Generate predictions to plot
- [Trajectories API](trajectories.md) - Generate trajectories
- [Complete Pipeline Tutorial](../tutorials/complete-pipeline.md)

