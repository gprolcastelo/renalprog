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
from sklearn.preprocessing import MinMaxScaler
from renalprog.modeling.train import VAE, NetworkReconstruction, train_vae
from renalprog.modeling.predict import apply_vae, generate_trajectory_data
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


# 4. Generate and plot trajectory
# Prepare scaler
scaler = MinMaxScaler()
scaler.fit(train_expr.values.T)

# Define trajectory (patient IDs in progression order)
trajectory_patients = ['TCGA-A1-001', 'TCGA-A2-002', 'TCGA-A3-003']

trajectory_data = generate_trajectory_data(
    vae_model=model,
    recnet_model=None,  # Optional reconstruction network
    trajectory=trajectory_patients,
    gene_data=test_expr.T,  # Transpose: genes as rows, patients as columns
    n_timepoints=50,
    device='cuda',
    scaler=scaler
)

# Get top variable genes for visualization
top_genes = trajectory_data.var(axis=0).nlargest(20).index

# Plot trajectory
plot_trajectory(
    trajectory=trajectory_data[top_genes],
    feature_names=top_genes.tolist(),
    output_path=output_dir / "trajectory_001.png",
    title="Disease Progression Trajectory"
)

print(f"All figures saved to {output_dir}")
```

## Classification Visualization

### plot_metrics

Visualize classification metrics across multiple model runs.

::: renalprog.plots.plot_metrics

**Example Usage:**

```python
from renalprog.plots import plot_metrics
import pandas as pd

# Classification results from multiple runs
results_df = pd.DataFrame({
    'Accuracy': [0.85, 0.87, 0.86, 0.88, 0.84],
    'Precision': [0.83, 0.85, 0.84, 0.86, 0.82],
    'Recall': [0.82, 0.84, 0.83, 0.85, 0.81],
    'F1-Score': [0.825, 0.845, 0.835, 0.855, 0.815],
    "Cohen's Kappa": [0.70, 0.74, 0.72, 0.76, 0.68]
})

# Create boxplot visualization
fig = plot_metrics(
    results_df,
    save_path="reports/figures/classification_metrics",
    title="XGBoost Classification Metrics"
)
```

### plot_trajectory_classification

Visualize classification probabilities along synthetic trajectories.

::: renalprog.plots.plot_trajectory_classification

**Example Usage:**

```python
from renalprog.plots import plot_trajectory_classification
import pandas as pd

# Load trajectory classification predictions
df_predictions = pd.read_csv("data/processed/trajectory_classifications.csv")

# Visualize how probabilities change over time
fig = plot_trajectory_classification(
    df_predictions,
    save_path="reports/figures/trajectory_classification",
    traj_type="test_to_test",
    n_timepoints=50,
    title="Disease Progression Classification"
)
```

**Complete Classification Workflow:**

```python
from renalprog.plots import plot_metrics, plot_trajectory_classification
from renalprog.modeling.train import classification_benchmark
from pathlib import Path

# Train multiple classifiers
results = []
models = []
for seed in range(10):
    metrics, model = classification_benchmark(
        X_train, y_train, X_test, y_test,
        seed=seed
    )
    results.append(metrics)
    models.append(model)

# Visualize metrics
results_df = pd.DataFrame(results)
plot_metrics(
    results_df,
    save_path="reports/figures/classification_metrics"
)

# Apply best model to trajectories
best_idx = results_df["Cohen's Kappa"].idxmax()
best_model = models[best_idx]

predictions = classify_trajectories(best_model, trajectory_data)

# Visualize trajectory classifications
plot_trajectory_classification(
    predictions,
    save_path="reports/figures/trajectory_classification",
    traj_type="train_to_train",
    n_timepoints=50
)
```

## See Also

- [Training API](training.md) - Generate training history
- [Prediction API](prediction.md) - Generate predictions to plot
- [Trajectories API](trajectories.md) - Generate trajectories
- [Step 5: Classification Tutorial](../tutorials/step5-classification.md) - Full classification workflow
- [Complete Pipeline Tutorial](../tutorials/complete-pipeline.md)

