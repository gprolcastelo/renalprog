# Utils API

Utility functions used across the RenalProg package.

## Overview

The utils module provides:

- Random seed setting for reproducibility
- Logging configuration
- Device selection (CPU/GPU)
- Data preprocessing utilities
- Helper functions

## Core Utilities

### set_seed

Set random seed for reproducibility across all libraries.

::: renalprog.utils.set_seed

**Example Usage:**

```python
from renalprog.utils import set_seed

# Set seed at the start of your script
set_seed(42)

# All random operations will be reproducible
import numpy as np
import torch

print(np.random.rand(5))  # Same output every time
print(torch.randn(5))     # Same output every time
```

**Libraries Affected:**

- `random` (Python standard library)
- `numpy`
- `torch` (CPU)
- `torch.cuda` (GPU)
- Sets `torch.backends.cudnn.deterministic = True`

### configure_logging

Configure logging with scientific output formatting.

::: renalprog.utils.configure_logging

**Example Usage:**

```python
from renalprog.utils import configure_logging
import logging

# Basic configuration (recommended)
configure_logging()

# Now logging works throughout the package
import renalprog.dataset as dataset
dataset.download_data()  # Will show progress logs

# Debug mode with timestamps
configure_logging(
    level=logging.DEBUG,
    format_string="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
)

# Custom file logging
file_handler = logging.FileHandler("pipeline.log")
configure_logging(handlers=[file_handler])
```

**Logging Levels:**

```python
import logging

configure_logging(level=logging.DEBUG)    # Show everything
configure_logging(level=logging.INFO)     # Default - show info and above
configure_logging(level=logging.WARNING)  # Show only warnings and errors
configure_logging(level=logging.ERROR)    # Show only errors
```

### get_device

Get appropriate device for PyTorch computation.

::: renalprog.utils.get_device

**Example Usage:**

```python
from renalprog.utils import get_device
import torch

# Automatically select CUDA if available
device = get_device()
print(f"Using device: {device}")  # cuda:0 or cpu

# Force CPU usage
device = get_device(force_cpu=True)
print(f"Using device: {device}")  # cpu

# Use device in model
model = VAE(input_dim=20000, mid_dim=1024, features=128)
model = model.to(device)

data = torch.randn(32, 20000).to(device)
output = model(data)
```

## Complete Script Template

Here's a template for a complete analysis script using all utilities:

```python
#!/usr/bin/env python3
"""
Complete analysis pipeline for RenalProg.

This script demonstrates proper usage of utility functions for
reproducibility and logging.
"""

import logging
from pathlib import Path
import pandas as pd
import torch

from renalprog.utils import set_seed, configure_logging, get_device
from renalprog.dataset import download_data, process_downloaded_data, create_train_test_split
from renalprog.features import preprocess_rnaseq
from renalprog.modeling.train import train_vae
from renalprog.modeling.predict import apply_vae, generate_trajectories
from renalprog.plots import plot_training_history

# ============================================================================
# Configuration
# ============================================================================

# Reproducibility
SEED = 42
set_seed(SEED)

# Logging
configure_logging(
    level=logging.INFO,
    format_string="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# Device
device = get_device()
logger.info(f"Using device: {device}")

# Paths
BASE_DIR = Path(".")
DATA_DIR = BASE_DIR / "data"
MODEL_DIR = BASE_DIR / "models" / "my_experiment"
OUTPUT_DIR = BASE_DIR / "reports"

# Ensure directories exist
for dir_path in [DATA_DIR, MODEL_DIR, OUTPUT_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Main Pipeline
# ============================================================================

def main():
    """Run complete analysis pipeline."""
    
    logger.info("Starting RenalProg analysis pipeline")
    
    # Step 1: Download data
    logger.info("Step 1: Downloading TCGA data")
    rnaseq_path, clinical_path, phenotype_path = download_data(
        destination=DATA_DIR / "raw"
    )
    
    # Step 2: Process for KIRC
    logger.info("Step 2: Processing KIRC data")
    rnaseq, clinical, phenotype = process_downloaded_data(
        rnaseq_path=rnaseq_path,
        clinical_path=clinical_path,
        phenotype_path=phenotype_path,
        cancer_type="KIRC",
        output_dir=DATA_DIR / "raw"
    )
    
    # Step 3: Preprocess
    logger.info("Step 3: Preprocessing gene expression")
    rnaseq_preprocessed = preprocess_rnaseq(
        rnaseq=rnaseq,
        output_dir=DATA_DIR / "interim" / "preprocessed"
    )
    
    # Step 4: Train/test split
    logger.info("Step 4: Creating train/test split")
    create_train_test_split(
        rnaseq=rnaseq_preprocessed,
        clinical=clinical,
        test_size=0.2,
        random_state=SEED,
        output_dir=DATA_DIR / "interim" / "split"
    )
    
    # Step 5: Load split data
    logger.info("Step 5: Loading split data")
    train_expr = pd.read_csv(
        DATA_DIR / "interim" / "split" / "train_expression.tsv",
        sep="\t", index_col=0
    )
    test_expr = pd.read_csv(
        DATA_DIR / "interim" / "split" / "test_expression.tsv",
        sep="\t", index_col=0
    )
    
    # Step 6: Train VAE
    logger.info("Step 6: Training VAE")
    history, model, checkpoints = train_vae(
        train_data=train_expr.values,
        val_data=test_expr.values,
        input_dim=train_expr.shape[1],
        mid_dim=1024,
        features=128,
        output_dir=MODEL_DIR,
        n_epochs=100,
        batch_size=32,
        learning_rate=1e-3,
        device=device,
        use_scheduler=True,
        early_stopping_patience=20
    )
    
    # Step 7: Plot training history
    logger.info("Step 7: Visualizing training")
    plot_training_history(
        history=history,
        output_path=OUTPUT_DIR / "figures" / "training_history.png"
    )
    
    # Step 8: Generate latent representations
    logger.info("Step 8: Generating latent representations")
    results = apply_vae(
        model=model,
        data=test_expr.values,
        device=device
    )
    
    # Step 9: Visualize latent space
    logger.info("Step 9: Visualizing latent space")
    clinical_test = pd.read_csv(
        DATA_DIR / "interim" / "split" / "test_clinical.tsv",
        sep="\t", index_col=0
    )
    
    
    logger.info("Pipeline completed successfully!")
    logger.info(f"Results saved to {OUTPUT_DIR}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        raise
```

## Best Practices

### 1. Always Set Seed First

```python
from renalprog.utils import set_seed

# First line of your script
set_seed(42)
```

### 2. Configure Logging Early

```python
from renalprog.utils import configure_logging

# Second line of your script
configure_logging()
```

### 3. Check Device Availability

```python
from renalprog.utils import get_device

device = get_device()
if device.type == 'cuda':
    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
else:
    print("Running on CPU")
```

### 4. Handle Errors Gracefully

```python
import logging
from renalprog.utils import configure_logging

configure_logging()
logger = logging.getLogger(__name__)

try:
    # Your code here
    results = some_function()
except Exception as e:
    logger.error(f"Analysis failed: {e}", exc_info=True)
    raise
```

### 5. Use Context Managers

```python
import torch
from renalprog.utils import get_device

device = get_device()

# Inference mode (faster, uses less memory)
model.eval()
with torch.no_grad():
    output = model(data.to(device))
```

## Environment Variables

You can control behavior via environment variables:

```bash
# Force CPU usage
export CUDA_VISIBLE_DEVICES=""
python my_script.py

# Use specific GPU
export CUDA_VISIBLE_DEVICES=1
python my_script.py

# Limit threads
export OMP_NUM_THREADS=4
python my_script.py
```

## Reproducibility Checklist

For fully reproducible results:

- ✅ Set random seed with `set_seed()`
- ✅ Use fixed `random_state` parameters
- ✅ Set `torch.backends.cudnn.deterministic = True`
- ✅ Document package versions
- ✅ Save configuration parameters
- ✅ Version control code
- ✅ Track data provenance

```python
from renalprog.utils import set_seed
import torch
import numpy as np
import pandas as pd
import json

# Reproducibility
set_seed(42)

# Save configuration
config = {
    'seed': 42,
    'torch_version': torch.__version__,
    'numpy_version': np.__version__,
    'pandas_version': pd.__version__,
    'cuda_available': torch.cuda.is_available(),
    'cudnn_deterministic': torch.backends.cudnn.deterministic,
    'cudnn_benchmark': torch.backends.cudnn.benchmark
}

with open('config.json', 'w') as f:
    json.dump(config, f, indent=2)
```

## See Also

- [Configuration API](config.md) - Project configuration
- [Complete Pipeline Tutorial](../tutorials/complete-pipeline.md)
- [Contributing Guide](../contributing/guidelines.md)

