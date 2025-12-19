# Quick Start Tutorial

Get started with `renalprog` in 10 minutes! This tutorial demonstrates the core functionality using example data.

## Installation

If you haven't installed `renalprog` yet:

```bash
# Clone repository
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog

# Create environment with mamba/conda
mamba env create -f environment.yml
mamba activate renalprog

# Install package
pip install -e .
```

See the [Installation Guide](../getting-started.md) for detailed instructions.

## Quick Example

Let's run a minimal example that demonstrates the complete pipeline:

### 1. Import Required Modules

```python
from renalprog import config, dataset, features, modeling, plots
from renalprog.modeling import VAE, generate_trajectories
import pandas as pd
import numpy as np
import torch
```

### 2. Load Example Data

```python
# Load preprocessed KIRC data
rnaseq_path = 'data/processed/rnaseq_maha.csv'
clinical_path = 'data/processed/clinical.csv'

# Load RNA-seq data
rnaseq = pd.read_csv(rnaseq_path, index_col=0)
clinical = pd.read_csv(clinical_path, index_col=0)

print(f"RNA-seq data shape: {rnaseq.shape}")
print(f"Samples: {rnaseq.shape[0]}, Genes: {rnaseq.shape[1]}")
```

**Expected output:**
```
RNA-seq data shape: (498, 5000)
Samples: 498, Genes: 5000
```

### 3. Prepare Data for VAE

```python
# Select important genes (optional - for faster training)
n_genes = 1000  # Use top 1000 most variable genes
gene_vars = rnaseq.var(axis=0)
top_genes = gene_vars.nlargest(n_genes).index
rnaseq_subset = rnaseq[top_genes]

# Convert to PyTorch tensor
X = torch.FloatTensor(rnaseq_subset.values)
print(f"Training data shape: {X.shape}")
```

### 4. Train a Simple VAE

```python
# Initialize VAE
vae = VAE(
    input_dim=X.shape[1],
    latent_dim=10,
    hidden_dims=[256, 128],
    seed=42
)

# Train
vae.fit(
    X,
    epochs=50,
    batch_size=32,
    learning_rate=1e-3,
    verbose=True
)
```

**Expected output:**
```
Epoch [10/50], Loss: 1234.56
Epoch [20/50], Loss: 987.65
Epoch [30/50], Loss: 765.43
Epoch [40/50], Loss: 654.32
Epoch [50/50], Loss: 598.21
Training complete!
```

!!! tip "Training Time"
    On a modern CPU, this should take 2-5 minutes. With a GPU, it takes under 1 minute.

### 5. Visualize Latent Space

```python
# Encode data to latent space
with torch.no_grad():
    latent = vae.encode(X).numpy()

# Create DataFrame with latent coordinates and stages
latent_df = pd.DataFrame(
    latent,
    index=rnaseq_subset.index,
    columns=[f'Z{i}' for i in range(latent.shape[1])]
)
latent_df['stage'] = clinical.loc[latent_df.index, 'ajcc_pathologic_tumor_stage']

# Plot with UMAP
from renalprog.plots import plot_latent_space

plot_latent_space(
    latent_df,
    color_by='stage',
    method='umap',
    save_path='reports/figures/quickstart_latent_umap.png'
)
```

### 6. Generate a Trajectory

```python
# Select early and late stage samples
early_samples = clinical[clinical['stage_binary'] == 'Early'].index
late_samples = clinical[clinical['stage_binary'] == 'Late'].index

# Get latent representations
early_latent = latent[rnaseq_subset.index.isin(early_samples)]
late_latent = latent[rnaseq_subset.index.isin(late_samples)]

# Pick random start and end points
start_point = early_latent[0]
end_point = late_latent[0]

# Generate trajectory (linear interpolation in latent space)
n_steps = 10
trajectory_latent = np.linspace(start_point, end_point, n_steps)

# Decode back to gene expression
trajectory_latent_tensor = torch.FloatTensor(trajectory_latent)
with torch.no_grad():
    trajectory_expression = vae.decode(trajectory_latent_tensor).numpy()

# Create DataFrame
trajectory_df = pd.DataFrame(
    trajectory_expression,
    columns=top_genes,
    index=[f'timepoint_{i}' for i in range(n_steps)]
)

print(f"Generated trajectory with {n_steps} timepoints")
print(trajectory_df.head())
```

### 7. Identify Changing Genes

```python
# Calculate fold change along trajectory
first_tp = trajectory_df.iloc[0]
last_tp = trajectory_df.iloc[-1]

# Log2 fold change
fc = np.log2((last_tp + 1) / (first_tp + 1))
fc_sorted = fc.abs().sort_values(ascending=False)

print("Top 10 most changing genes:")
print(fc_sorted.head(10))
```

### 8. Quick Classification

```python
from renalprog.modeling.classification import train_stage_classifier

# Prepare labels
y = (clinical.loc[rnaseq_subset.index, 'stage_binary'] == 'Late').astype(int)

# Train classifier
clf, metrics, shap_values = train_stage_classifier(
    X=rnaseq_subset.values,
    y=y.values,
    feature_names=top_genes.tolist(),
    n_folds=5
)

print(f"Classification Accuracy: {metrics['test_accuracy']:.3f}")
print(f"ROC AUC: {metrics['test_roc_auc']:.3f}")
```

**Expected output:**
```
Classification Accuracy: 0.892
ROC AUC: 0.945
```

### 9. Save Results

```python
# Save trajectory
trajectory_df.to_csv('data/interim/quickstart_trajectory.csv')

# Save model
vae.save('models/quickstart_vae.pt')

# Save latent representation
latent_df.to_csv('data/interim/quickstart_latent.csv')

print("Results saved!")
```

## Complete Script

Here's the complete script you can run:

```python title="quickstart_example.py"
"""Quick start example for renalprog."""

from renalprog import dataset, features, modeling, plots
from renalprog.modeling import VAE
import pandas as pd
import numpy as np
import torch
import os

# Create output directories
os.makedirs('reports/figures', exist_ok=True)
os.makedirs('data/interim', exist_ok=True)
os.makedirs('models', exist_ok=True)

# 1. Load data
print("Loading data...")
rnaseq = pd.read_csv('data/processed/rnaseq_maha.csv', index_col=0)
clinical = pd.read_csv('data/processed/clinical.csv', index_col=0)

# 2. Prepare data
print("Preparing data...")
n_genes = 1000
gene_vars = rnaseq.var(axis=0)
top_genes = gene_vars.nlargest(n_genes).index
rnaseq_subset = rnaseq[top_genes]
X = torch.FloatTensor(rnaseq_subset.values)

# 3. Train VAE
print("Training VAE...")
vae = VAE(
    input_dim=X.shape[1],
    latent_dim=10,
    hidden_dims=[256, 128],
    seed=42
)
vae.fit(X, epochs=50, batch_size=32, learning_rate=1e-3, verbose=True)

# 4. Get latent representations
print("Encoding to latent space...")
with torch.no_grad():
    latent = vae.encode(X).numpy()

latent_df = pd.DataFrame(
    latent,
    index=rnaseq_subset.index,
    columns=[f'Z{i}' for i in range(latent.shape[1])]
)
latent_df['stage'] = clinical.loc[latent_df.index, 'ajcc_pathologic_tumor_stage']

# 5. Visualize
print("Creating visualizations...")
plots.plot_latent_space(
    latent_df,
    color_by='stage',
    method='umap',
    save_path='reports/figures/quickstart_latent_umap.png'
)

# 6. Generate trajectory
print("Generating trajectory...")
early_samples = clinical[clinical['stage_binary'] == 'Early'].index
late_samples = clinical[clinical['stage_binary'] == 'Late'].index

early_latent = latent[rnaseq_subset.index.isin(early_samples)]
late_latent = latent[rnaseq_subset.index.isin(late_samples)]

start_point = early_latent[0]
end_point = late_latent[0]

n_steps = 10
trajectory_latent = np.linspace(start_point, end_point, n_steps)
trajectory_latent_tensor = torch.FloatTensor(trajectory_latent)

with torch.no_grad():
    trajectory_expression = vae.decode(trajectory_latent_tensor).numpy()

trajectory_df = pd.DataFrame(
    trajectory_expression,
    columns=top_genes,
    index=[f'timepoint_{i}' for i in range(n_steps)]
)

# 7. Save results
print("Saving results...")
trajectory_df.to_csv('data/interim/quickstart_trajectory.csv')
vae.save('models/quickstart_vae.pt')
latent_df.to_csv('data/interim/quickstart_latent.csv')

print("Quick start complete! Check reports/figures/ for visualizations.")
```

Run it with:
```bash
python quickstart_example.py
```

## What You've Learned

In this quick start, you've:

- ✅ Loaded and prepared gene expression data
- ✅ Trained a VAE to learn latent representations
- ✅ Visualized the latent space
- ✅ Generated a synthetic progression trajectory
- ✅ Identified genes that change along the trajectory
- ✅ Built a classification model

## Next Steps

Now that you've seen the basics, dive deeper:

1. **[Complete Pipeline](complete-pipeline.md)**: Full analysis workflow from raw data
2. **[VAE Training](step2-vae-training.md)**: Advanced model configuration
3. **[Trajectory Generation](step4-trajectories.md)**: Multiple trajectory types
4. **[Enrichment Analysis](step6-enrichment.md)**: Pathway-level interpretation

## Troubleshooting

### Out of Memory

If you get OOM errors:
```python
# Reduce number of genes
n_genes = 500

# Reduce batch size
vae.fit(X, batch_size=16, ...)

# Use smaller hidden dimensions
vae = VAE(input_dim=X.shape[1], hidden_dims=[128, 64], ...)
```

### Slow Training

To speed up training:
```python
# Use GPU if available
device = 'cuda' if torch.cuda.is_available() else 'cpu'
vae = VAE(input_dim=X.shape[1], device=device)

# Reduce epochs for testing
vae.fit(X, epochs=20, ...)
```

### Import Errors

Make sure you've installed the package:
```bash
pip install -e .
```

And activated the environment:
```bash
mamba activate renalprog
```

