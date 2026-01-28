# Trajectories API

Functions for analyzing disease progression trajectories.

## Overview

The trajectories module provides analysis tools for:

- Trajectory network construction
- Patient connectivity analysis
- Temporal pathway enrichment
- Transition probability calculation
- Trajectory visualization

## Network Construction

### build_trajectory_network

Build directed graph of patient transitions.

::: renalprog.modeling.predict.build_trajectory_network

**Example Usage:**

```python
from renalprog.modeling.predict import build_trajectory_network
import pandas as pd
from pathlib import Path

# Load patient connections
connections = pd.read_csv("data/processed/patient_connections.csv")

# Build network
network = build_trajectory_network(
    connections=connections,
    output_path=Path("data/processed/trajectory_network.graphml")
)

print(f"Network has {network.number_of_nodes()} nodes")
print(f"Network has {network.number_of_edges()} edges")
```

### generate_trajectory_data

Generate complete trajectory dataset with metadata.

::: renalprog.modeling.predict.generate_trajectory_data

## Patient Connectivity

### link_patients_closest

Link patients using closest latent space neighbors.

::: renalprog.modeling.predict.link_patients_closest

**Example:**

```python
from renalprog.modeling.predict import link_patients_closest
import numpy as np

early_latent = np.random.randn(100, 128)
late_latent = np.random.randn(80, 128)

connections = link_patients_closest(
    latent_early=early_latent,
    latent_late=late_latent,
    patient_ids_early=['E001', 'E002', ...],
    patient_ids_late=['L001', 'L002', ...]
)

# Returns DataFrame with columns: early_patient, late_patient, distance
```

### link_patients_random

Link patients randomly (control method).

::: renalprog.modeling.predict.link_patients_random

## Transition Analysis

### calculate_all_possible_transitions

Calculate metrics for all possible patient transitions.

::: renalprog.modeling.predict.calculate_all_possible_transitions

**Example Usage:**

```python
from renalprog.modeling.predict import calculate_all_possible_transitions

# Calculate all transitions
transitions = calculate_all_possible_transitions(
    latent_early=early_latent,
    latent_late=late_latent,
    patient_ids_early=early_ids,
    patient_ids_late=late_ids,
    output_dir=Path("data/processed/transitions")
)

# Analyze transition patterns
print(transitions.describe())
```

## Dynamic Enrichment

### dynamic_enrichment_analysis

Perform pathway enrichment at each trajectory timepoint.

::: renalprog.modeling.predict.dynamic_enrichment_analysis

**Example Usage:**

```python
from renalprog.modeling.predict import dynamic_enrichment_analysis
from pathlib import Path

# Analyze pathway dynamics along trajectories
enrichment_results = dynamic_enrichment_analysis(
    trajectories=trajectory_gene_expression,  # Shape: (n_traj, n_steps, n_genes)
    gene_names=gene_list,
    pathway_file=Path("data/external/ReactomePathways.gmt"),
    output_dir=Path("reports/dynamic_enrichment")
)

# Results contain enrichment at each timepoint
for timepoint, results in enrichment_results.items():
    print(f"Timepoint {timepoint}: {len(results)} enriched pathways")
```

## Interpolation Methods

### interpolate_latent_linear

Linear interpolation between points.

::: renalprog.modeling.predict.interpolate_latent_linear

### interpolate_latent_spherical

Spherical interpolation (SLERP) for normalized spaces.

::: renalprog.modeling.predict.interpolate_latent_spherical

**Comparison:**

```python
from renalprog.modeling.predict import (
    interpolate_latent_linear,
    interpolate_latent_spherical
)
import numpy as np

z_start = np.random.randn(1, 128)
z_end = np.random.randn(1, 128)

# Linear interpolation
traj_linear = interpolate_latent_linear(z_start, z_end, n_steps=50)

# Spherical interpolation (preserves norm better)
traj_spherical = interpolate_latent_spherical(z_start, z_end, n_steps=50)

# Spherical is preferred for normalized latent spaces
```

## Visualization

### plot_trajectory

Visualize individual trajectory.

::: renalprog.plots.plot_trajectory

**Example:**

```python
from renalprog.plots import plot_trajectory
from pathlib import Path

# Plot single trajectory
plot_trajectory(
    trajectory=trajectory_data[0],  # Shape: (n_steps, n_features)
    feature_names=selected_genes,
    output_path=Path("reports/figures/trajectory_example.png"),
    title="Disease Progression Trajectory"
)
```

## See Also

- [Prediction API](prediction.md) - Apply trained models
- [Plots API](plots.md) - Visualization functions
- [Complete Pipeline Tutorial](../tutorials/complete-pipeline.md)

