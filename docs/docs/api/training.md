# Training API

Complete training pipeline for VAE models with checkpointing and monitoring.

## Overview

The training module provides high-level functions for:

- Complete VAE training workflow
- Automatic checkpointing
- Training history visualization
- Early stopping
- Learning rate scheduling

## Main Training Function

### train_vae

The main training function that orchestrates the entire VAE training process.

::: renalprog.modeling.train.train_vae

## Key Features

### Automatic Checkpointing

The training function automatically saves:

- Model state dict
- Optimizer state
- Training history
- Configuration parameters

Checkpoints are saved when:
- Validation loss improves (best model)
- At regular intervals (every `checkpoint_every` epochs)
- After training completes (final model)

### Early Stopping

Training stops automatically if validation loss doesn't improve for `early_stopping_patience` epochs. This prevents overfitting and saves computation time.

### Learning Rate Scheduling

When `use_scheduler=True`, the learning rate is reduced when validation loss plateaus:

```python
scheduler = ReduceLROnPlateau(
    optimizer,
    mode='min',
    factor=0.5,
    patience=10,
    verbose=True
)
```

### Cyclical KL Annealing

The KL divergence weight β is gradually increased using cyclical annealing to prevent posterior collapse:

::: renalprog.modeling.train.frange_cycle_linear

## Training History

The training function returns a dictionary with:

| Key | Description |
|-----|-------------|
| `train_loss` | Training loss per epoch |
| `val_loss` | Validation loss per epoch |
| `train_recon` | Training reconstruction loss per epoch |
| `val_recon` | Validation reconstruction loss per epoch |
| `train_kl` | Training KL divergence per epoch |
| `val_kl` | Validation KL divergence per epoch |
| `learning_rates` | Learning rate per epoch |

## Complete Example

```python
import pandas as pd
from pathlib import Path
from renalprog.modeling.train import train_vae
from renalprog.plots import plot_training_history
from renalprog.utils import set_seed, configure_logging

# Configure
configure_logging()
set_seed(42)

# Load data
train_expr = pd.read_csv("data/interim/split/train_expression.tsv", sep="\t", index_col=0)
test_expr = pd.read_csv("data/interim/split/test_expression.tsv", sep="\t", index_col=0)

# Train VAE
history, best_model, checkpoints = train_vae(
    train_data=train_expr.values,
    val_data=test_expr.values,
    input_dim=train_expr.shape[1],
    mid_dim=1024,
    features=128,
    output_dir=Path("models/my_vae"),
    n_epochs=100,
    batch_size=32,
    learning_rate=1e-3,
    use_scheduler=True,
    use_checkpoint=True,
    checkpoint_every=10,
    early_stopping_patience=20,
    device='cuda'
)

# Plot results
plot_training_history(
    history,
    output_path=Path("reports/figures/training_history.png")
)

# Load best model for inference
best_model.eval()
import torch
with torch.no_grad():
    reconstruction, mu, log_var, z = best_model(
        torch.FloatTensor(test_expr.values).to(device)
    )
    
print(f"Best validation loss: {min(history['val_loss']):.4f}")
print(f"Final learning rate: {history['learning_rates'][-1]:.6f}")
```

## Checkpointing API

For manual checkpoint management:

::: renalprog.modeling.checkpointing

### save_checkpoint

Save model checkpoint with metadata.

::: renalprog.modeling.checkpointing.save_model_config

### load_checkpoint

Load model from checkpoint.

::: renalprog.modeling.checkpointing.load_model_config

## See Also

- [Models API](models.md) - VAE architectures
- [Prediction API](prediction.md) - Using trained models
- [Configuration](config.md) - Training hyperparameters
- [Complete Pipeline Tutorial](../tutorials/complete-pipeline.md)

