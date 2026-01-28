# Visualization Tutorial

This guide shows you how to use the visualization functions in `renalprog` to create plots for your analysis.

## Overview

The `renalprog` package provides two main modules for visualization:

1. **`plots` module**: General-purpose plotting functions using Plotly (interactive plots)
2. **`enrichment` module**: Specialized heatmaps for pathway enrichment analysis using Matplotlib

!!! tip "Interactive vs Static"
    - **Plotly plots** (`plots` module): Interactive HTML plots with hover tooltips, zoom, and pan
    - **Matplotlib plots** (`enrichment` module): Static publication-ready heatmaps
    
    Both save in multiple formats: HTML, PNG, PDF, and SVG

!!! warning "Dependencies"
    All functions using plotly require the `kaleido` package for static image export and Google Chrome.

!!! note "Plotly outputs"
    The `plots` module functions save interactive plots in HTML, as well as static images (PNG, PDF, SVG).

## Plots Module

The `plots` module provides interactive visualizations using Plotly.

### Available Functions

| Function                       | Purpose | Output                                                            |
|--------------------------------|---------|-------------------------------------------------------------------|
| `save_plot()`                  | Helper function to save Plotly figures | Save in multiple formats (.html, .png, .pdf, and .svg by default) |
| `plot_training_history()`      | VAE training progress | Loss curves over epochs                                           |
| `plot_reconstruction_losses()` | Reconstruction network training | Train/test loss curves                                            |
| `plot_trajectory()`            | Gene expression trajectories | Line plots across timepoints                                      |
| `plot_pca_variance()`          | PCA variance explained | Bar chart of PCs' variances                                       |
| `plot_umap_plotly()`           | Dimensionality reduction | 2D/3D UMAP scatter plot                                           |
| `plot_enrichment_results()`    | Enrichment scores | Bar chart of pathway scores                                       |

### Importing the functions

```python
from renalprog.plots import (
    plot_training_history,
    plot_reconstruction_losses,
    plot_trajectory,
    plot_umap_plotly,
    
)
from renalprog.enrichment import plot_enrichment_results
from pathlib import Path
```

---

### 1. Training History Visualization

Visualize VAE training progress with train and validation losses.

```python
# After training a VAE
from renalprog.modeling.train import train_vae_with_postprocessing

# ... training code ...
vae_model, network, vae_history, reconstruction_history = train_vae_with_postprocessing(
    X_train=X_train,
    X_test=X_test,
    vae_config=vae_config,
    # ... other parameters ...
)

# Plot VAE training history
plot_training_history(
    history=vae_history,
    save_path=Path('reports/figures/vae_training.png'),
    title='VAE Training History',
    log_scale=False  # Set to True for log scale
)
```

**Expected Output:**
- Interactive plot showing train/validation loss over epochs
- Separate curves for total loss, reconstruction loss, and KL divergence
- Saved in HTML, PNG, PDF, and SVG formats

---

### 2. Reconstruction Network Training

Visualize the training of the reconstruction network (post-processing).

```python
# Plot reconstruction network training
plot_reconstruction_losses(
    loss_train=reconstruction_history["train_loss"],
    loss_test=reconstruction_history["test_loss"],
    save_path=Path('reports/figures/reconstruction_training.png'),
    title='Reconstruction Network Training',
    show_best_epoch=True  # Highlight epoch with lowest validation loss
)
```

**Features:**
- Shows train vs test loss
- Optionally marks the best epoch
- Helps identify overfitting

---

### 3. Gene Expression Trajectories

Plot how gene expression changes along disease progression trajectories.

```python
import pandas as pd

# Load trajectory data (timepoints × genes)
trajectory_df = pd.read_csv('trajectory_expression.csv', index_col=0)

# Plot multiple genes
genes_to_plot = ['TP53', 'VEGFA', 'HIF1A', 'VHL']

plot_trajectory(
    trajectory_df=trajectory_df,
    genes=genes_to_plot,
    save_path=Path('reports/figures/gene_trajectories.png'),
    title='Gene Expression Trajectories',
    normalize=True,  # Normalize each gene to [0, 1]
    show_markers=True,
    colormap='Viridis'
)
```

**Output:**
- Line plot with one trace per gene
- X-axis: Timepoints (pseudo-time from early to late)
- Y-axis: Expression level
- (HTML) Interactive hover showing exact values

---

### 4. UMAP Visualization

Create 2D or 3D UMAP plots to visualize high-dimensional data.

```python
import pandas as pd

# Load expression data and clinical metadata
data = pd.read_csv('preprocessed_rnaseq.csv', index_col=0)  # samples × genes
clinical = pd.read_csv('clinical_data.csv', index_col=0)    # samples × metadata

# Ensure data is samples × genes (transpose if needed)
if data.shape[0] > data.shape[1]:
    data = data.T

# Create UMAP plot colored by disease stage
plot_umap_plotly(
    data=data,
    clinical=clinical,
    colors_dict={'early': '#6495ed', 'late': '#b22222'},
    n_components=2,  # 2D plot (use 3 for 3D)
    save_fig=True,
    save_as='reports/figures/umap_by_stage',
    seed=2023,
    title='UMAP: Original Data by Stage',
    show=True
)
```

**Parameters:**
- `data`: Gene expression matrix (samples × genes)
- `clinical`: Clinical metadata with 'stage' column (values: 'early' or 'late')
- `colors_dict`: Mapping of stage to color
- `n_components`: 2 for 2D plot, 3 for 3D plot
- `seed`: Random seed for reproducibility

**Output:**
- Interactive scatter plot
- Hover shows sample ID and stage
- Can zoom, pan, and (in 3D) rotate


---

## Enrichment Module Heatmaps

The `enrichment` module provides specialized functions for creating publication-quality pathway enrichment heatmaps.

### Main Functions

#### 1. `generate_pathway_heatmap()`

Creates multiple heatmaps showing pathway regulation across disease progression.

**Purpose**: Visualize how biological pathways are regulated over pseudo-time from early to late disease stages.

**What it does**:
- Aggregates enrichment scores across all trajectories
- Creates 4-5 heatmaps:
  1. Top 50 most changing pathways
  2. Top 50 upregulated pathways
  3. Top 50 downregulated pathways
  4. High-level Reactome pathways
  5. Literature-curated pathways (optional)

#### 2. `plot_heatmap_regulation()` (internal helper)

Creates individual heatmap with custom styling.

---

### Complete Example: Pathway Enrichment Heatmaps

```python
from renalprog.enrichment import generate_pathway_heatmap
import pandas as pd
from pathlib import Path

# ============================================================================
# Step 1: Load Enrichment Data
# ============================================================================
# Load GSEA results from trajectory analysis
# Expected columns: [Patient, Idx, Transition, NAME, ES, NES, FDR q-val]
enrichment_df = pd.read_csv('trajectory_enrichment_results.csv')

print(f"Enrichment data shape: {enrichment_df.shape}")
print(f"Columns: {enrichment_df.columns.tolist()}")
print(f"Number of pathways: {enrichment_df['NAME'].nunique()}")
print(f"Number of timepoints: {enrichment_df['Idx'].nunique()}")

# ============================================================================
# Step 2: Generate Heatmaps
# ============================================================================
output_dir = Path('reports/figures/pathway_heatmaps')

heatmap_data, figures = generate_pathway_heatmap(
    enrichment_df=enrichment_df,
    output_dir=output_dir,
    fdr_threshold=0.05,      # Only include significant pathways
    colorbar=True,           # Show colorbar
    legend=False,            # Hide legend (optional)
    yticks_fontsize=12,      # Font size for pathway names
    show=False               # Don't display plots (just save)
)

print(f"\nGenerated {len(figures)} heatmaps:")
for name in figures.keys():
    print(f"  - {name}")

# ============================================================================
# Step 3: Examine Results
# ============================================================================
# The heatmap_data DataFrame contains summed NES values
print(f"\nHeatmap data shape: {heatmap_data.shape}")
print(f"Pathways (rows): {heatmap_data.shape[0]}")
print(f"Timepoints (columns): {heatmap_data.shape[1]}")

# View top pathways at first and last timepoint
first_timepoint = heatmap_data.iloc[:, 0].sort_values(ascending=False)
last_timepoint = heatmap_data.iloc[:, -1].sort_values(ascending=False)

print("\nTop 5 pathways at early stage:")
print(first_timepoint.head())

print("\nTop 5 pathways at late stage:")
print(last_timepoint.head())
```

---

### Understanding the Output

#### Heatmap Structure

```
Rows: Pathway names (e.g., "Immune System", "DNA Repair")
Columns: Pseudo-time (early → late)
Values: Sum of NES (Normalized Enrichment Score)
  - Positive (red): Pathway upregulated
  - Negative (blue): Pathway downregulated
  - Zero (white): No change
The values can be any other metric deemed appropriate.
```

#### Color Scheme

The heatmaps use a **diverging colormap** centered at zero:

- **Red**: Upregulated pathways (positive NES)
- **White**: No change (NES ≈ 0)
- **Blue**: Downregulated pathways (negative NES)

The color scale is **symmetric** around zero, making it easy to identify regulation direction.

---

### Output Files

After running `generate_pathway_heatmap()`, you'll find:

```
reports/figures/pathway_heatmaps/
├── heatmap_top50_changing.png          # Top 50 most dynamic pathways
├── heatmap_top50_upregulated.png       # Top 50 upregulated pathways
├── heatmap_top50_downregulated.png     # Top 50 downregulated pathways
├── heatmap_selected_high_level.png     # Reactome high-level pathways
└── heatmap_selected_literature.png     # Literature-curated pathways (if present)
```

Each heatmap is saved in **PNG** format at high resolution (300 DPI).

---

### Customization Options

#### Adjust Significance Threshold

```python
# More stringent: only FDR < 0.01
heatmap_data, figures = generate_pathway_heatmap(
    enrichment_df=enrichment_df,
    output_dir='reports/figures/',
    fdr_threshold=0.01,  # More stringent
    # ... other parameters
)
```

#### Change Visual Appearance

```python
# Larger font for pathway names
heatmap_data, figures = generate_pathway_heatmap(
    enrichment_df=enrichment_df,
    output_dir='reports/figures/',
    yticks_fontsize=14,  # Larger font
    colorbar=True,       # Show colorbar
    legend=True,         # Show legend
    show=True            # Display plots interactively
)
```

#### Access Individual Figures

```python
# Generate heatmaps
heatmap_data, figures = generate_pathway_heatmap(...)

# Access specific figure
fig_changing = figures['heatmap_top50_changing']
fig_upregulated = figures['heatmap_top50_upregulated']

# Further customize with matplotlib
import matplotlib.pyplot as plt

fig_changing.suptitle('My Custom Title', fontsize=20)
fig_changing.savefig('custom_heatmap.pdf', dpi=300, bbox_inches='tight')
plt.close(fig_changing)
```

---

### Understanding the Heatmap Data

The returned `heatmap_data` DataFrame can be used for further analysis:

```python
# Generate heatmaps
heatmap_data, figures = generate_pathway_heatmap(...)

# Analyze pathway dynamics
# Calculate change from early to late
pathway_changes = heatmap_data.iloc[:, -1] - heatmap_data.iloc[:, 0]
pathway_changes_sorted = pathway_changes.sort_values(ascending=False)

print("Most increasing pathways:")
print(pathway_changes_sorted.head(10))

print("\nMost decreasing pathways:")
print(pathway_changes_sorted.tail(10))

# Find pathways active throughout
mean_nes = heatmap_data.mean(axis=1)
std_nes = heatmap_data.std(axis=1)

consistently_high = mean_nes[mean_nes > 0.5].sort_values(ascending=False)
print("\nConsistently upregulated pathways:")
print(consistently_high.head(10))

# Export for further analysis
heatmap_data.to_csv('pathway_nes_matrix.csv')
```

---

## See Also

- [API Documentation - Plots](../api/plots.md)
- [API Documentation - Enrichment](../api/enrichment.md)
- [Enrichment Analysis Tutorial](step6-enrichment.md)

---
