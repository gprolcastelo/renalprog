# Trajectory Generation Pipeline

This document describes the synthetic cancer progression trajectory generation pipeline in renalprog.

## Overview

The trajectory generation pipeline creates synthetic cancer progression paths by interpolating between early-stage and late-stage cancer samples in the VAE's latent space. This allows us to model intermediate states of cancer progression that may not be directly observable in the data.

## Biological Motivation

Cancer progression is a continuous process, but clinical data typically captures discrete snapshots (stages I-IV). Synthetic trajectories help us:

1. **Model Intermediate States**: Generate hypothetical intermediate stages between observed samples
2. **Identify Progression Markers**: Find genes and pathways that change during progression
3. **Understand Dynamics**: Study the temporal ordering of molecular events
4. **Generate Hypotheses**: Predict which pathways are activated during progression

## Pipeline Steps

### 1. Patient Connection Creation

Identify pairs of patients to connect:

```python
from renalprog.modeling.trajectories import create_patient_connections

connections = create_patient_connections(
    clinical_df=clinical_data,
    source_stage="early",  # Stage I-II
    target_stage="late",   # Stage III-IV
    method="all_to_all"    # Connect every early to every late
)

# Each connection: (source_patient_id, target_patient_id)
```

Connection strategies:
- **all_to_all**: Every early patient → every late patient
- **nearest_neighbor**: Early patients → nearest late patients in latent space
- **matched**: Match on clinical covariates (age, sex, etc.)

### 2. Latent Space Interpolation

Generate intermediate points in VAE latent space:

```python
from renalprog.modeling.predict import generate_trajectories

trajectories = generate_trajectories(
    model=vae_model,
    source_samples=early_stage_data,
    target_samples=late_stage_data,
    n_steps=50,           # Number of intermediate timepoints
    method="linear"       # Interpolation method
)

# Output shape: (n_trajectories, n_steps, n_genes)
```

Interpolation methods:
- **linear**: Straight line in latent space (default)
- **spherical**: Interpolate along geodesic on hypersphere
- **bezier**: Smooth curve using Bezier interpolation

### 3. Gene Expression Reconstruction

Decode latent trajectories to gene expression space:

```python
from renalprog.modeling.predict import reconstruct_from_latent

# For each timepoint, decode latent vector to gene expression
gene_expression = []
for t in range(n_steps):
    z_t = trajectories[:, t, :]  # Latent vectors at time t
    x_t = vae_model.decode(z_t)  # Reconstruct gene expression
    gene_expression.append(x_t)

# Shape: (n_trajectories, n_steps, n_genes)
gene_expression = np.stack(gene_expression, axis=1)
```

### 4. Quality Control

Validate trajectory quality:

```python
from renalprog.modeling.trajectories import validate_trajectories

quality_metrics = validate_trajectories(
    trajectories=gene_expression,
    early_samples=X_early,
    late_samples=X_late,
    vae_model=vae_model
)

# Metrics:
# - reconstruction_error: How well VAE reconstructs samples
# - monotonicity: Do progression markers increase/decrease consistently?
# - smoothness: Are trajectories smooth or erratic?
# - stage_separation: Do endpoints match early/late distributions?
```

### 5. Control Trajectory Generation

Generate negative controls for comparison:

```python
from renalprog.modeling.trajectories import generate_control_trajectories

# Random noise trajectories
noise_controls = generate_control_trajectories(
    n_trajectories=len(trajectories),
    n_steps=n_steps,
    n_genes=n_genes,
    method="gaussian_noise",
    noise_scale=1.0
)

# Shuffled gene trajectories
shuffled_controls = generate_control_trajectories(
    trajectories=trajectories,
    method="shuffle_genes"
)
```

Control types:
- **gaussian_noise**: Random Gaussian noise
- **shuffle_genes**: Shuffle gene order independently
- **shuffle_time**: Randomize temporal ordering
- **reverse**: Reverse early→late to late→early

## Output Files

The trajectory pipeline generates:

```
data/interim/trajectories/
├── early_to_late/
│   ├── trajectory_metadata.csv      # Patient pairs, stages
│   ├── latent_trajectories.npy      # Latent space paths
│   ├── gene_expression.npy          # Reconstructed expression
│   └── quality_metrics.json         # Validation metrics
├── controls/
│   ├── noise_controls.npy           # Random noise
│   └── shuffled_controls.npy        # Shuffled genes
└── visualizations/
    ├── latent_space_trajectories.png
    ├── sample_trajectories.png
    └── progression_heatmap.png
```

## Usage Example

### Complete Trajectory Workflow

```python
from renalprog.modeling.predict import generate_trajectories
from renalprog.modeling.trajectories import (
    create_patient_connections,
    validate_trajectories,
    generate_control_trajectories
)
from renalprog.config import get_dated_dir, INTERIM_DATA_DIR
import numpy as np

# 1. Load trained VAE
from renalprog.modeling.train import load_vae
vae = load_vae("models/20251216_VAE_KIRC/best_model.pth")

# 2. Load data
X_early = np.load("data/interim/samples/early_stage.npy")
X_late = np.load("data/interim/samples/late_stage.npy")

# 3. Create connections
connections = create_patient_connections(
    clinical_df=clinical,
    source_stage="early",
    target_stage="late",
    method="all_to_all"
)

# 4. Generate trajectories
trajectories = generate_trajectories(
    model=vae,
    source_samples=X_early,
    target_samples=X_late,
    n_steps=50,
    method="linear"
)

# 5. Validate quality
metrics = validate_trajectories(
    trajectories=trajectories,
    early_samples=X_early,
    late_samples=X_late,
    vae_model=vae
)

# 6. Generate controls
controls = generate_control_trajectories(
    n_trajectories=len(trajectories),
    n_steps=50,
    n_genes=trajectories.shape[-1],
    method="gaussian_noise"
)

# 7. Save results
output_dir = get_dated_dir(INTERIM_DATA_DIR, "synthetic_trajectories")
np.save(f"{output_dir}/trajectories.npy", trajectories)
np.save(f"{output_dir}/controls.npy", controls)

print(f"Generated {len(trajectories)} trajectories")
print(f"Quality metrics: {metrics}")
```

### Script Usage

```bash
# Run trajectory generation pipeline
python scripts/pipeline_steps/4_trajectories.py \
  --vae_model models/20251216_VAE_KIRC/best_model.pth \
  --data_dir data/interim/20251216_train_test_split \
  --output_dir data/interim/20251216_synthetic_trajectories \
  --n_steps 50 \
  --method linear \
  --generate_controls
```

## Visualization

### Latent Space Trajectories

Visualize trajectories in 2D latent space (using UMAP/t-SNE):

```python
from renalprog.plots import plot_latent_trajectories

fig = plot_latent_trajectories(
    latent_trajectories=z_trajectories,
    early_latent=z_early,
    late_latent=z_late,
    reduction_method="umap"
)
fig.savefig("reports/figures/latent_trajectories.png")
```

### Gene Expression Heatmaps

Show how genes change along trajectories:

```python
from renalprog.plots import plot_trajectory_heatmap

fig = plot_trajectory_heatmap(
    trajectories=trajectories[:10],  # First 10 trajectories
    gene_names=important_genes,
    n_genes_show=50
)
fig.savefig("reports/figures/trajectory_heatmap.png")
```

### Individual Trajectory Plots

Plot selected genes over time:

```python
from renalprog.plots import plot_individual_trajectory

genes_of_interest = ["VHL", "PBRM1", "SETD2", "BAP1"]

fig = plot_individual_trajectory(
    trajectory=trajectories[0],
    gene_names=all_genes,
    highlight_genes=genes_of_interest
)
fig.savefig("reports/figures/individual_trajectory.png")
```

## Advanced Topics

### Custom Interpolation Methods

Implement custom interpolation:

```python
from renalprog.modeling.trajectories import BaseInterpolator

class CustomInterpolator(BaseInterpolator):
    def interpolate(self, z_start, z_end, n_steps):
        """Custom interpolation between latent vectors."""
        # Your custom interpolation logic
        # For example, cubic spline:
        from scipy.interpolate import CubicSpline
        
        t = np.linspace(0, 1, n_steps)
        cs = CubicSpline([0, 1], np.vstack([z_start, z_end]))
        return cs(t)

# Use custom interpolator
trajectories = generate_trajectories(
    model=vae,
    source_samples=X_early,
    target_samples=X_late,
    n_steps=50,
    interpolator=CustomInterpolator()
)
```

### Trajectory Clustering

Identify different progression patterns:

```python
from sklearn.cluster import KMeans
from renalprog.modeling.trajectories import extract_trajectory_features

# Extract features from trajectories
features = extract_trajectory_features(
    trajectories,
    method="slope"  # Gene expression slopes
)

# Cluster trajectories
kmeans = KMeans(n_clusters=3, random_state=42)
clusters = kmeans.fit_predict(features)

# Analyze each cluster
for i in range(3):
    cluster_trajs = trajectories[clusters == i]
    print(f"Cluster {i}: {len(cluster_trajs)} trajectories")
```

### Multi-Stage Trajectories

Generate trajectories across multiple stages:

```python
# Generate trajectories: Stage I → Stage II → Stage III → Stage IV
stages = ["stage_i", "stage_ii", "stage_iii", "stage_iv"]
multi_stage_trajs = []

for i in range(len(stages) - 1):
    trajs = generate_trajectories(
        model=vae,
        source_samples=samples[stages[i]],
        target_samples=samples[stages[i+1]],
        n_steps=20,
        method="linear"
    )
    multi_stage_trajs.append(trajs)

# Concatenate all segments
full_trajectories = np.concatenate(multi_stage_trajs, axis=1)
```

## Interpretation Guidelines

### What Do Trajectories Represent?

- **Not Real Patients**: Trajectories are hypothetical, model-generated paths
- **Biological Plausibility**: Validated by enrichment analysis and clinical markers
- **Statistical Averages**: Represent general trends, not individual variation
- **Hypothesis Generation**: Use to identify candidates for experimental validation

### Common Pitfalls

1. **Over-interpretation**: Don't assume every trajectory detail is biologically meaningful
2. **Model Artifacts**: VAE may introduce interpolation artifacts
3. **Missing Biology**: Model can't capture unknown mechanisms
4. **Data Bias**: Trajectories reflect training data distributions

### Best Practices

1. **Generate Controls**: Always compare to random/shuffled controls
2. **Validate Endpoints**: Ensure trajectories start/end at correct stages
3. **Check Monotonicity**: Progression markers should change consistently
4. **Cross-validate**: Compare to independent validation datasets
5. **Biological Validation**: Prioritize findings that match known biology

## Performance Considerations

### Memory Usage

For large datasets:
```python
# Generate trajectories in batches
batch_size = 100
all_trajectories = []

for i in range(0, len(connections), batch_size):
    batch_connections = connections[i:i+batch_size]
    batch_trajs = generate_trajectories(
        model=vae,
        source_samples=X_early[batch_connections[:, 0]],
        target_samples=X_late[batch_connections[:, 1]],
        n_steps=50
    )
    all_trajectories.append(batch_trajs)
    
trajectories = np.concatenate(all_trajectories, axis=0)
```

### Computational Speed

- Linear interpolation: Fast (~1-10 ms per trajectory)
- Spherical interpolation: Medium (~10-50 ms per trajectory)
- Bezier interpolation: Slow (~50-200 ms per trajectory)

For 1000 trajectories with 50 steps:
- Linear: ~10 seconds
- Spherical: ~1 minute
- Bezier: ~5 minutes

## See Also

- [VAE Training](./tutorials/step2-vae-training.md) - Train the VAE model
- [Classification](./CLASSIFICATION.md) - Classify trajectory timepoints
- [Enrichment Analysis](./ENRICHMENT_ANALYSIS.md) - Find enriched pathways
- [API Reference](./api/trajectories.md) - Detailed API documentation

## References

1. White, T. (2016). Sampling Generative Networks. arXiv preprint arXiv:1609.04468.

2. Shoemake, K. (1985). Animating rotation with quaternion curves. ACM SIGGRAPH Computer Graphics, 19(3), 245-254.

3. Way, G. P., & Greene, C. S. (2018). Extracting a biologically relevant latent space from cancer transcriptomes with variational autoencoders. Pacific Symposium on Biocomputing, 23, 80-91.

4. Ding, J., et al. (2018). Interpretable dimensionality reduction of single cell transcriptome data with deep generative models. Nature Communications, 9(1), 1-13.

