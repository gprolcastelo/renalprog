# API Reference

Complete API reference for the `renalprog` package. All functions, classes, and modules are documented here.

## Package Structure

```
renalprog/
├── __init__.py              # Package initialization
├── config.py                # Configuration and paths
├── dataset.py               # Data loading and splitting
├── features.py              # Preprocessing and feature engineering
├── plots.py                 # Visualization functions
├── modeling/                # Machine learning models
│   ├── __init__.py
│   ├── models.py           # VAE architectures
│   ├── train.py            # Training functions
│   ├── predict.py          # Prediction and inference
│   ├── trajectories.py     # Trajectory generation
│   ├── classification.py   # XGBoost classification
│   └── enrichment.py       # GSEA integration
└── utils/                   # Utility functions
    ├── __init__.py
    ├── logging.py          # Logging configuration
    └── helpers.py          # Helper functions
```

## Quick Links

### Core Modules
- [Configuration](config.md) - Paths and settings
- [Dataset](dataset.md) - Data loading and processing
- [Features](features.md) - Preprocessing and filtering

### Modeling
- [VAE Models](models.md) - Variational autoencoders
- [Training](training.md) - Model training
- [Prediction](prediction.md) - Inference and generation
- [Trajectories](trajectories.md) - Synthetic progression
- [Classification](classification.md) - Stage prediction
- [Enrichment](enrichment.md) - Pathway analysis

### Visualization
- [Plots](plots.md) - All plotting functions

### Utilities
- [Utils](utils.md) - Helper functions and logging

## Installation

```bash
pip install renalprog
```

Or for development:

```bash
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog
pip install -e .
```

## Basic Usage

### Import the Package

```python
import renalprog
from renalprog import config, dataset, features, modeling, plots
```

### Load Data

```python
# Load preprocessed data
rnaseq = dataset.load_data('data/processed/rnaseq_maha.csv')
clinical = dataset.load_data('data/processed/clinical.csv')

# Create train/test split
X_train, X_test, y_train, y_test = dataset.create_train_test_split(
    rnaseq_path='data/processed/rnaseq_maha.csv',
    clinical_path='data/processed/clinical.csv',
    test_size=0.2,
    seed=2023
)
```

### Train a Model

```python
from renalprog.modeling import VAE, train_vae

# Initialize VAE
vae = VAE(
    input_dim=X_train.shape[1],
    latent_dim=256,
    hidden_dims=[512, 256]
)

# Train
vae.fit(X_train, epochs=100, batch_size=32)

# Save
vae.save('models/my_vae.pt')
```

### Generate Trajectories

```python
from renalprog.modeling import generate_trajectory_data
from sklearn.preprocessing import MinMaxScaler

# Prepare trajectory (list of patient IDs in progression order)
trajectory = ['TCGA-A1-001', 'TCGA-A2-002', 'TCGA-A3-003']

# Use the same scaler that was used during VAE training
scaler = MinMaxScaler()
scaler.fit(X_train.T)  # Fit on training data (genes as features)

# Generate interpolated gene expression along trajectory
trajectory_data = generate_trajectory_data(
    vae_model=vae,
    recnet_model=reconstruction_network,  # Optional, can be None
    trajectory=trajectory,
    gene_data=gene_expression_df,  # DataFrame with genes as rows, patients as columns
    n_timepoints=50,
    interpolation_method='linear',  # or 'spherical'
    device='cpu',
    scaler=scaler
)
```

### Run Classification

```python
from renalprog.modeling.classification import train_stage_classifier

clf, metrics, shap_values = train_stage_classifier(
    X=X_train,
    y=y_train,
    feature_names=gene_names,
    n_folds=5
)
```

### Visualize

```python
from renalprog import plots

# Plot training history
plots.plot_training_history(
    history,
    save_path='figures/training.png'
)
```

## Module Details

### renalprog.config

Configuration module with paths and constants.

**Key Components**:

- `PATHS`: Dictionary of all project paths
- `VAEConfig`: VAE hyperparameter configuration
- `get_dated_dir()`: Create dated output directories
- `KIRCPaths`: KIRC-specific data paths

**Example**:
```python
from renalprog.config import PATHS, VAEConfig

# Access paths
data_dir = PATHS['data']
models_dir = PATHS['models']

# Configure VAE
config = VAEConfig()
config.LATENT_DIM = 256
config.EPOCHS = 600
```

See: [Configuration Reference](config.md)

---

### renalprog.dataset

Data loading, processing, and splitting utilities.

**Key Functions**:

- `download_data()`: Download TCGA data from Xena
- `process_downloaded_data()`: Process raw TCGA files
- `load_data()`: Load CSV files
- `create_train_test_split()`: Stratified train/test split

**Example**:
```python
from renalprog import dataset

# Download data
rnaseq, clinical, pheno = dataset.download_data(
    destination='data/raw',
    cancer_type='KIRC'
)

# Create split
X_train, X_test, y_train, y_test = dataset.create_train_test_split(
    rnaseq_path='data/processed/rnaseq.csv',
    clinical_path='data/processed/clinical.csv',
    test_size=0.2,
    seed=2023
)
```

See: [Dataset Reference](dataset.md)

---

### renalprog.features

Preprocessing and feature engineering.

**Key Functions**:

- `preprocess_rnaseq()`: Filter and normalize gene expression
- `filter_low_expression()`: Remove lowly expressed genes
- `detect_outliers()`: Mahalanobis distance outlier detection
- `normalize()`: Various normalization methods

**Example**:
```python
from renalprog import features

# Preprocess RNA-seq data
processed, info = features.preprocess_rnaseq(
    data=rnaseq,
    filter_expression=True,
    detect_outliers=True,
    mean_threshold=0.5,
    var_threshold=0.5,
    alpha=0.05
)
```

See: [Features Reference](features.md)

---

### renalprog.modeling

Machine learning models and training.

**Submodules**:

- `models`: VAE architectures
- `train`: Training functions
- `predict`: Inference and generation
- `trajectories`: Synthetic trajectory generation
- `classification`: XGBoost stage classification
- `enrichment`: GSEA integration

**Example**:
```python
from renalprog.modeling import VAE, train_vae, generate_trajectory_data

# Train VAE
vae, history = train_vae(
    X_train=X_train,
    X_test=X_test,
    config=vae_config
)

# Generate trajectory data
trajectory_data = generate_trajectory_data(
    vae_model=vae,
    recnet_model=None,  # Optional reconstruction network
    trajectory=['patient1', 'patient2', 'patient3'],
    gene_data=gene_expression_df,
    n_timepoints=50
)
```

See:
- [Models Reference](models.md)
- [Training Reference](training.md)
- [Trajectories Reference](trajectories.md)
- [Classification Reference](classification.md)
- [Enrichment Reference](enrichment.md)

---

### renalprog.plots

Visualization functions for all analysis steps.

**Key Functions**:

- `plot_training_history()`: Loss curves
- `plot_reconstruction()`: Original vs reconstructed
- `plot_trajectories()`: Trajectory visualization
- `plot_metrics()`: Classification metrics boxplots
- `plot_trajectory_classification()`: Trajectory classification over time
- `plot_enrichment_heatmap()`: Pathway enrichment

See: [Plots Reference](plots.md)

---

### renalprog.utils

Utility functions and helpers.

**Key Functions**:

- `configure_logging()`: Set up logging
- `set_seed()`: Set random seeds for reproducibility
- `Timer`: Context manager for timing code
- `check_file_exists()`: File validation

**Example**:
```python
from renalprog.utils import configure_logging, set_seed, Timer

# Configure logging
configure_logging(level='INFO')

# Set random seed
set_seed(2023)

# Time code execution
with Timer("Training"):
    vae.fit(X_train, epochs=100)
```

See: [Utils Reference](utils.md)

---

## Advanced Topics

### Custom VAE Architectures

```python
from renalprog.modeling.models import BaseVAE
import torch.nn as nn

class CustomVAE(BaseVAE):
    def __init__(self, input_dim, latent_dim):
        super().__init__(input_dim, latent_dim)
        
        # Custom encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 1024),
            nn.BatchNorm1d(1024),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(1024, 512),
            nn.ReLU(),
        )
        
        self.fc_mu = nn.Linear(512, latent_dim)
        self.fc_logvar = nn.Linear(512, latent_dim)
        
        # Custom decoder
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 512),
            nn.ReLU(),
            nn.Linear(512, 1024),
            nn.ReLU(),
            nn.Linear(1024, input_dim)
        )
```

### Custom Preprocessing

```python
from renalprog.features import preprocess_rnaseq

def custom_preprocess(data, **kwargs):
    """Custom preprocessing pipeline."""
    
    # Step 1: Log transform
    data_log = np.log2(data + 1)
    
    # Step 2: Quantile normalization
    data_norm = quantile_normalize(data_log)
    
    # Step 3: Standard preprocessing
    data_processed, info = preprocess_rnaseq(
        data_norm,
        **kwargs
    )
    
    return data_processed, info
```

### Batch Processing

```python
from renalprog.modeling import generate_trajectory_data
from sklearn.preprocessing import MinMaxScaler
import glob

# Process multiple experiments
for experiment_dir in glob.glob('data/interim/experiment_*'):
    vae = VAE.load(f'{experiment_dir}/vae_model.pt')
    recnet = NetworkReconstruction.load(f'{experiment_dir}/recnet_model.pt')
    
    # Load the scaler used during training
    scaler = MinMaxScaler()
    scaler.fit(X_train.T)
    
    trajectory_data = generate_trajectory_data(
        vae_model=vae,
        recnet_model=recnet,
        trajectory=['early_patient', 'mid_patient', 'late_patient'],
        gene_data=gene_df,
        n_timepoints=50,
        scaler=scaler
    )
    
    trajectory_data.to_csv(f'{experiment_dir}/trajectory_data.csv')
```

## Type Hints

The package uses type hints throughout:

```python
from typing import Tuple, Optional, Dict, List
import pandas as pd
import numpy as np
import torch

def preprocess_rnaseq(
    data: pd.DataFrame,
    filter_expression: bool = True,
    mean_threshold: float = 0.5,
    var_threshold: float = 0.5,
    detect_outliers: bool = True,
    alpha: float = 0.05,
    seed: Optional[int] = None
) -> Tuple[pd.DataFrame, Dict[str, any]]:
    """Preprocess RNA-seq data."""
    pass
```

## Error Handling

Common exceptions:

```python
from renalprog.dataset import load_data

try:
    data = load_data('nonexistent_file.csv')
except FileNotFoundError as e:
    print(f"File not found: {e}")
except pd.errors.EmptyDataError:
    print("File is empty")
except Exception as e:
    print(f"Unexpected error: {e}")
```

## Performance

### GPU Acceleration

```python
# Use GPU if available
device = 'cuda' if torch.cuda.is_available() else 'cpu'

vae = VAE(
    input_dim=5000,
    latent_dim=256,
    device=device
)
```

### Parallel Processing

```python
from renalprog.enrichment import run_gsea_parallel

# Use multiple cores
run_gsea_parallel(
    deseq_dir='data/processed/deseq',
    n_threads=8
)
```

## Testing

Run the test suite:

```bash
# All tests
pytest

# Specific module
pytest tests/test_dataset.py

# With coverage
pytest --cov=renalprog --cov-report=html
```

## Contributing

See [Contributing Guidelines](../contributing/guidelines.md) for development information.

## Version History

- **v0.1.0** (2024-12): Initial release
  - VAE training pipeline
  - Trajectory generation
  - XGBoost classification
  - GSEA integration

## License

Apache License 2.0. See [LICENSE](../license.md).

## Citation

```bibtex
@software{renalprog2024,
  author = {Prol-Castelo, Guillermo},
  title = {renalprog: Cancer Progression Forecasting with Generative AI},
  year = {2024},
  url = {https://github.com/gprolcastelo/renalprog}
}
```

