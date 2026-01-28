# Performance Optimization

The VAE model architecture can easily be tuned using [Ax](https://ax.dev/) or [Optuna](https://optuna.org/).

## Hyperparameter Optimization with Ax

This example shows how to optimize VAE hyperparameters using Adaptive Experimentation Platform (Ax).

### Minimal Working Example

```python
"""
Minimal example for VAE hyperparameter optimization using Ax.

This script optimizes:
- Latent dimension size
- Middle layer dimension
- Learning rate
- Beta (KL weight)
"""

from ax.service.ax_client import AxClient
from ax.service.utils.instantiation import ObjectiveProperties
import torch
import numpy as np
from pathlib import Path

from renalprog import dataset
from renalprog.config import VAEConfig
from renalprog.modeling.train import train_vae

# ============================================================================
# 1. Setup: Load Data
# ============================================================================
print("Loading data...")

# Load preprocessed data
X_train, X_test, y_train, y_test, _, _ = dataset.create_train_test_split(
    rnaseq_path=Path('data/interim/preprocessed_KIRC_data/preprocessed_rnaseq.csv'),
    clinical_path=Path('data/interim/preprocessed_KIRC_data/clinical_data.csv'),
    test_size=0.2,
    seed=2023,
    output_dir=Path('data/interim/train_test_split')
)

input_dim = X_train.shape[1]
print(f"Input dimension: {input_dim}")
print(f"Training samples: {X_train.shape[0]}")

# ============================================================================
# 2. Define Evaluation Function
# ============================================================================
def evaluate_vae(parameterization):
    """
    Train VAE with given hyperparameters and return validation loss.
    
    Args:
        parameterization: Dict with hyperparameters from Ax
        
    Returns:
        Dict with 'val_loss' metric
    """
    # Extract hyperparameters
    latent_dim = parameterization['latent_dim']
    mid_dim = parameterization['mid_dim']
    learning_rate = parameterization['learning_rate']
    beta_ratio = parameterization['beta_ratio']
    
    print(f"\nTrying: latent_dim={latent_dim}, mid_dim={mid_dim}, "
          f"lr={learning_rate:.4f}, beta={beta_ratio:.2f}")
    
    # Configure VAE
    vae_config = VAEConfig()
    vae_config.INPUT_DIM = input_dim
    vae_config.LATENT_DIM = latent_dim
    vae_config.MID_DIM = mid_dim
    vae_config.LEARNING_RATE = learning_rate
    vae_config.BETA_RATIO = 0.5
    vae_config.BETA_CYCLES = 3  # Single cycle for speed
    vae_config.EPOCHS = 200 * vae_config.BETA_CYCLES  
    vae_config.BATCH_SIZE = 32
    
    # Train VAE
    try:
        vae_model, history = train_vae(
            X_train=X_train,
            X_test=X_test,
            config=vae_config,
            save_dir=None,  # Don't save intermediate models
            force_cpu=False # Running on GPU is recommended
        )
        
        # Get average validation loss for the last 20 epochs
        val_loss = np.mean(history['val_loss'][-20:])
        
        print(f"  → Validation loss: {val_loss:.4f}")
        
        return {'val_loss': (val_loss, 0.0)}  # (mean, sem)
        
    except Exception as e:
        print(f"  → Training failed: {e}")
        return {'val_loss': (float('inf'), 0.0)}

# ============================================================================
# 3. Setup Ax Client
# ============================================================================
ax_client = AxClient()

ax_client.create_experiment(
    name="vae_optimization",
    parameters=[
        {
            "name": "latent_dim",
            "type": "range",
            "bounds": [64, 512],
            "value_type": "int",
            "log_scale": True,  # Search in log space
        },
        {
            "name": "mid_dim",
            "type": "range",
            "bounds": [256, 2048],
            "value_type": "int",
            "log_scale": True,
        },
        {
            "name": "learning_rate",
            "type": "range",
            "bounds": [1e-4, 1e-2],
            "value_type": "float",
            "log_scale": True,
        },
    ],
    objectives={
        "val_loss": ObjectiveProperties(minimize=True)
    },
)

# ============================================================================
# 4. Run Optimization
# ============================================================================
print("\n" + "="*80)
print("Starting Bayesian Optimization")
print("="*80)

n_trials = 20  # Number of configurations to try

for trial_idx in range(n_trials):
    print(f"\n{'='*80}")
    print(f"Trial {trial_idx + 1}/{n_trials}")
    print(f"{'='*80}")
    
    # Get next parameters to try
    parameters, trial_index = ax_client.get_next_trial()
    
    # Evaluate
    result = evaluate_vae(parameters)
    
    # Report results back to Ax
    ax_client.complete_trial(trial_index=trial_index, raw_data=result)

# ============================================================================
# 5. Get Best Configuration
# ============================================================================
print("\n" + "="*80)
print("OPTIMIZATION COMPLETE")
print("="*80)

best_parameters, metrics = ax_client.get_best_parameters()

print("\nBest hyperparameters found:")
for param, value in best_parameters.items():
    print(f"  {param}: {value}")

print(f"\nBest validation loss: {metrics[0]['val_loss']:.4f}")

# ============================================================================
# 6. Train Final Model with Best Configuration
# ============================================================================
print("\n" + "="*80)
print("Training final model with best hyperparameters...")
print("="*80)

final_config = VAEConfig()
final_config.INPUT_DIM = input_dim
final_config.LATENT_DIM = best_parameters['latent_dim']
final_config.MID_DIM = best_parameters['mid_dim']
final_config.LEARNING_RATE = best_parameters['learning_rate']
final_config.BETA_RATIO = best_parameters['beta_ratio']
final_config.EPOCHS = 600  # Full training
final_config.BATCH_SIZE = 8
final_config.BETA_CYCLES = 3

final_model, final_history = train_vae(
    X_train=X_train,
    X_test=X_test,
    config=final_config,
    save_dir=Path('models/optimized_vae'),
    force_cpu=False
)

print(f"\n Final model saved to: models/optimized_vae/")
print(f" Final validation loss: {final_history['val_loss'][-1]:.4f}")
```

### See Also

- [Ax Documentation](https://ax.dev/)
- [Bayesian Optimization in Pytorch](https://botorch.org/)
- [VAE Configuration API](../api/config.md)
- [VAE Training API](../api/models.md)

