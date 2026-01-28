# Configuration API

The `config` module provides centralized configuration management for the RenalProg pipeline, including paths, hyperparameters, and project-wide constants.

## Overview

The configuration module is organized into several components:

- **Path Management**: Centralized directory structure for data, models, and outputs
- **Preprocessing Configuration**: Parameters for data preprocessing
- **VAE Configuration**: Hyperparameters for VAE model training
- **Trajectory Configuration**: Settings for trajectory generation
- **Classification Configuration**: Parameters for survival classification
- **Enrichment Configuration**: Settings for pathway enrichment analysis


## Path Structure

The module defines a comprehensive directory structure:

```python
from renalprog.config import PATHS

# Access common paths
data_dir = PATHS['data']
models_dir = PATHS['models']
figures_dir = PATHS['figures']
```

### Available Paths

| Key | Description |
|-----|-------------|
| `root` | Project root directory |
| `data` | Main data directory |
| `raw` | Raw, immutable data |
| `interim` | Intermediate processed data |
| `processed` | Final, canonical data sets |
| `external` | External data sources (pathways, gene lists) |
| `models` | Trained models and checkpoints |
| `reports` | Analysis reports |
| `figures` | Generated figures and plots |
| `notebooks` | Jupyter notebooks |
| `references` | Reference materials |
| `scripts` | Pipeline scripts |

## Configuration Classes

### PreprocessingConfig

Configuration for data preprocessing steps.

::: renalprog.config.PreprocessingConfig
    options:
      show_root_heading: true
      members:
        - mean_threshold
        - var_threshold
        - min_sample_fraction
        - outlier_alpha
        - test_size
        - random_state

**Example Usage:**

```python
from renalprog.config import PreprocessingConfig

# Access default preprocessing parameters
config = PreprocessingConfig()
print(f"Mean threshold: {config.mean_threshold}")
print(f"Outlier alpha: {config.outlier_alpha}")
print(f"Test split ratio: {config.test_size}")
```

### VAEConfig

Configuration for VAE model architecture and training.

::: renalprog.config.VAEConfig
    options:
      show_root_heading: true

**Example Usage:**

```python
from renalprog.config import VAEConfig

# Use default VAE configuration
config = VAEConfig()

# Or customize for your experiment
custom_config = VAEConfig(
    mid_dim=256,
    latent_dim=10,
    learning_rate=0.0001,
    num_epochs=500,
    batch_size=64
)
```

### TrajectoryConfig

Configuration for trajectory generation from trained VAE models.

::: renalprog.config.TrajectoryConfig
    options:
      show_root_heading: true

**Example Usage:**

```python
from renalprog.config import TrajectoryConfig

config = TrajectoryConfig(
    num_trajectories=5000,
    trajectory_length=100,
    noise_scale=0.1
)
```

### ClassificationConfig

Configuration for survival classification models.

::: renalprog.config.ClassificationConfig
    options:
      show_root_heading: true

### EnrichmentConfig

Configuration for pathway enrichment analysis.

::: renalprog.config.EnrichmentConfig
    options:
      show_root_heading: true

## Utility Functions

### get_dated_dir

Create a dated directory for organizing time-stamped outputs.

::: renalprog.config.get_dated_dir

**Example:**

```python
from renalprog.config import get_dated_dir, PATHS

# Create a dated directory for today's model outputs
model_dir = get_dated_dir(PATHS['models'], prefix='VAE_KIRC')
# Returns: models/20251218_VAE_KIRC/
```

## KIRC-Specific Paths

The module provides a `KIRCPaths` class for managing KIRC dataset-specific file locations:

```python
from renalprog.config import KIRCPaths

# Access KIRC-specific data paths
rnaseq = KIRCPaths.RNASEQ_RAW
clinical = KIRCPaths.CLINICAL_RAW
pathways = KIRCPaths.REACTOME_PATHWAYS
```

## Best Practices

### Using Configuration in Scripts

Always import configuration at the module level:

```python
from renalprog.config import PATHS, VAEConfig, PreprocessingConfig

def main():
    # Access paths
    data_dir = PATHS['interim']
    
    # Use configuration objects
    vae_config = VAEConfig()
    preproc_config = PreprocessingConfig()
```

### Creating Custom Configurations

For experiments, create configuration variants:

```python
from renalprog.config import VAEConfig

# Base configuration
base_config = VAEConfig()

# Experiment variations
configs = {
    'small': VAEConfig(mid_dim=128, latent_dim=5),
    'medium': VAEConfig(mid_dim=256, latent_dim=10),
    'large': VAEConfig(mid_dim=512, latent_dim=20)
}
```

### Saving Configuration

Always save configuration with trained models:

```python
from renalprog.modeling.checkpointing import save_model_config

config = VAEConfig(experiment_name='my_experiment')
save_model_config(config, output_dir)
```

## See Also

- [Dataset API](dataset.md) - Data loading and preprocessing
- [Training API](training.md) - Model training with configuration