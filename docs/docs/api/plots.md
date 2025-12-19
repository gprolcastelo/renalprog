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

### plot_loss_landscape

Visualize loss landscape.

::: renalprog.plots.plot_loss_landscape

## Latent Space Visualization

### plot_latent_space

2D visualization of latent representations.

::: renalprog.plots.plot_latent_space

**Example:**

```python
from renalprog.plots import plot_latent_space
from renalprog.modeling.predict import apply_vae
import pandas as pd
from pathlib import Path

# Get latent representations
results = apply_vae(model, data, device='cuda')
latent = results['latent']

# Load labels
clinical = pd.read_csv("data/interim/split/test_clinical.tsv", sep="\t", index_col=0)

# Plot latent space colored by stage
plot_latent_space(
    latent=latent,
    labels=clinical['stage'],
    output_path=Path("reports/figures/latent_space_by_stage.png"),
    method='umap',  # or 'pca', 'tsne'
    title="Latent Space Representation"
)
```

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

## Gene Expression Visualization

### plot_gene_expression_heatmap

Heatmap of gene expression patterns.

::: renalprog.plots.plot_gene_expression_heatmap

**Example:**

```python
from renalprog.plots import plot_gene_expression_heatmap
import pandas as pd

# Load expression data
expr = pd.read_csv("data/interim/split/test_expression.tsv", sep="\t", index_col=0)

# Select top variable genes
var = expr.var(axis=0).sort_values(ascending=False)
top_genes = var.head(50).index

# Plot heatmap
plot_gene_expression_heatmap(
    expression=expr[top_genes],
    row_labels=expr.index.tolist(),
    col_labels=top_genes.tolist(),
    output_path=Path("reports/figures/expression_heatmap.png"),
    cluster_rows=True,
    cluster_cols=True
)
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

## Enrichment Visualization

### plot_enrichment_results

Visualize pathway enrichment results.

::: renalprog.plots.plot_enrichment_results

**Example:**

```python
from renalprog.plots import plot_enrichment_results
import pandas as pd

# Load enrichment results
enrichment = pd.read_csv("reports/enrichment/combined_results.csv")

# Plot top pathways
plot_enrichment_results(
    enrichment_df=enrichment,
    output_path=Path("reports/figures/enrichment_barplot.png"),
    top_n=20,
    metric='fdr',
    title="Top Enriched Pathways"
)
```

## Classification Visualization

### plot_confusion_matrix

Visualize classification performance.

::: renalprog.plots.plot_confusion_matrix

**Example:**

```python
from renalprog.plots import plot_confusion_matrix
from sklearn.metrics import confusion_matrix

# After classification
y_true = test_labels
y_pred = classifier.predict(test_features)
cm = confusion_matrix(y_true, y_pred)

plot_confusion_matrix(
    confusion_matrix=cm,
    class_names=['Non-progressing', 'Progressing'],
    output_path=Path("reports/figures/confusion_matrix.png"),
    normalize=True
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
    plot_latent_space,
    plot_gene_expression_heatmap,
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

plot_latent_space(
    latent=results['latent'],
    labels=clinical['stage'],
    output_path=output_dir / "latent_space_umap.png",
    method='umap'
)

plot_umap_plotly(
    latent=results['latent'],
    labels=clinical['stage'],
    sample_names=clinical.index.tolist(),
    title="Interactive Latent Space"
).write_html(output_dir / "latent_space_interactive.html")

# 4. Plot expression heatmap
var = test_expr.var(axis=0).sort_values(ascending=False)
top_genes = var.head(50).index

plot_gene_expression_heatmap(
    expression=test_expr[top_genes],
    row_labels=test_expr.index.tolist(),
    col_labels=top_genes.tolist(),
    output_path=output_dir / "expression_heatmap.png",
    cluster_rows=True,
    cluster_cols=True
)

# 5. Generate and plot trajectories
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

## Customization

All plotting functions accept matplotlib/plotly parameters for customization:

```python
from renalprog.plots import plot_latent_space

plot_latent_space(
    latent=latent,
    labels=labels,
    output_path=output_path,
    figsize=(12, 8),         # Custom figure size
    cmap='viridis',          # Custom colormap
    alpha=0.6,               # Transparency
    s=50,                    # Point size
    title="Custom Title",
    xlabel="Component 1",
    ylabel="Component 2"
)
```

## Publication-Quality Figures

For publication:

```python
import matplotlib.pyplot as plt
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'sans-serif',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12
})

# Then use plotting functions
plot_latent_space(...)
```

## See Also

- [Training API](training.md) - Generate training history
- [Prediction API](prediction.md) - Generate predictions to plot
- [Trajectories API](trajectories.md) - Generate trajectories
- [Complete Pipeline Tutorial](../tutorials/complete-pipeline.md)

