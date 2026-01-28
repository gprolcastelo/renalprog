# Step 2: VAE Training Pipeline

This guide explains how to train Variational Autoencoder (VAE) models and Reconstruction Networks for the RenalProg pipeline.

## Overview

The VAE training pipeline performs the following steps:

1. **Data Loading**: Load preprocessed RNAseq data from Step 1
2. **Train/Test Split**: Create stratified train/test splits
3. **VAE Training**: Train a Variational Autoencoder with cyclic β-annealing
4. **Reconstruction Network Training**: Train a postprocessing network to refine VAE outputs
5. **Visualization**: Plot training histories and save models

The pipeline trains two sequential models:
- **VAE**: Learns a compressed latent representation of gene expression data
- **Reconstruction Network**: Refines the VAE's reconstructions back to original gene space

## Prerequisites

Before running the VAE training pipeline, ensure you have:

- **Preprocessed data**: Completed Step 1 data processing
- **Python environment**: With PyTorch, numpy, pandas installed
- **Sufficient compute**: GPU recommended but not required (can use CPU with `use_cpu=True`)


## Usage

### Basic Usage

```bash
python scripts/pipeline_steps/2_models.py
```

This will:
- Load preprocessed KIRC data
- Create train/test split (80/20)
- Train VAE with default hyperparameters
- Train Reconstruction Network
- Save models and training plots

### Customizing Parameters

Edit the script to modify:
- **Cancer type**: Change `cancer_type = "KIRC"` to `"BRCA"` or other types
- **Data paths**: Update `path_rnaseq` and `path_clinical` for your data
- **Model architecture**: Modify VAE and Reconstruction Network dimensions
- **Training parameters**: Adjust epochs, batch size, learning rates, etc.

## Configuration Parameters

### Data Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cancer_type` | `"KIRC"` | Cancer type identifier |
| `path_rnaseq` | `"data/interim/preprocessed_KIRC/preprocessed_rnaseq.csv"` | Path to preprocessed RNAseq data |
| `path_clinical` | `"data/interim/preprocessed_KIRC/clinical_data.csv"` | Path to clinical data |
| `test_size` | `0.2` | Fraction of data to use for testing (20%) |

### VAE Architecture

| Parameter | Default | Description |
|-----------|---------|-------------|
| `INPUT_DIM` | `X_train.shape[1]` | Number of input features (genes), auto-detected |
| `MID_DIM` | `512` | Dimension of middle hidden layer |
| `LATENT_DIM` | `256` | Dimension of latent space |
| `BETA_CYCLES` | `3` | Number of β-annealing cycles |
| `EPOCHS` | `600` | Total training epochs (200 per cycle × 3 cycles) |
| `BETA_RATIO` | `0.5` | Fraction of each cycle spent increasing β |
| `BATCH_SIZE` | `8` | Training batch size |

### Reconstruction Network Architecture

| Parameter | Default | Description |
|-----------|---------|-------------|
| `recnet_dims` | `[INPUT_DIM, 3512, 824, 3731, INPUT_DIM]` | Layer dimensions for reconstruction network |
| `recnet_lr` | `1e-4` | Learning rate for reconstruction network |
| `recnet_epochs` | `1000` | Number of training epochs |
| `batch_recnet` | `8` | Batch size for reconstruction network |

### Hardware

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_cpu` | `True` | Force CPU usage even if GPU is available |

## Training Process

### Step 1: Train/Test Split

The pipeline creates a stratified split preserving the distribution of cancer stages:

```python
X_train, X_test, y_train, y_test = create_train_test_split(
    rnaseq_path=path_rnaseq,
    clinical_path=path_clinical,
    test_size=0.2,
    seed=2023
)
```

Output saved to: `data/interim/YYYYMMDD_train_test_split/`

### Step 2: VAE Training

The VAE uses **cyclic β-annealing** to balance reconstruction and regularization:

- **β-annealing**: Gradually increases β from 0 to 1 during first half of each cycle
- **Purpose**: Prevents posterior collapse and improves latent space quality
- **Cycles**: 3 cycles of 200 epochs each (600 total epochs)

**VAE Architecture:**
```
Input (genes) → Encoder → Latent (256D) → Decoder → Reconstruction (genes)
                  ↓
              Sampling
```

### Step 3: Reconstruction Network Training

After VAE training, a postprocessing network refines the reconstructions:

- **Purpose**: Correct systematic reconstruction errors from VAE
- **Architecture**: Deep feedforward network with hidden layers
- **Training**: 1000 epochs on VAE-encoded then decoded data

## Output Files

### Model Files

```
models/YYYYMMDD_models_KIRC/
├── vae/
│   ├── final_model.pth                    # Trained VAE weights
│   ├── config.json                        # VAE architecture configuration
│   └── vae_training_history.png           # Training loss plots
├── reconstruction_network.pth             # Trained reconstruction network
├── network_dims.csv                       # Network layer dimensions
└── reconstruction_network_history.png     # Reconstruction training plots
```

### Train/Test Split

```
data/interim/YYYYMMDD_train_test_split/
├── X_train.csv                            # Training gene expression data
├── X_test.csv                             # Test gene expression data
├── y_train.csv                            # Training labels (stages)
├── y_test.csv                             # Test labels (stages)
├── train_indices.csv                      # Training sample indices
└── test_indices.csv                       # Test sample indices
```

## Model Architecture Details

### VAE Configuration

The default VAE architecture for KIRC (~8500 genes):

```
Encoder:
  Input: 8516 genes
  ↓
  Hidden: 512 dimensions
  ↓
  Latent: 256 dimensions (μ and σ)

Decoder:
  Latent: 256 dimensions
  ↓
  Hidden: 512 dimensions
  ↓
  Output: 8516 genes
```

### Reconstruction Network

The default reconstruction network:

```
Input: 8516 genes (from VAE decoder)
  ↓
Layer 1: 3512 neurons
  ↓
Layer 2: 824 neurons
  ↓
Layer 3: 3731 neurons
  ↓
Output: 8516 genes (refined reconstruction)
```

## Training Strategies

### β-Annealing

The pipeline uses cyclic β-annealing:

1. **Warmup Phase** (first 50% of cycle):
   - β increases from 0 to 1
   - Model focuses on reconstruction first, then regularization

2. **Full Training** (second 50% of cycle):
   - β = 1 (standard VAE loss)
   - Model balances reconstruction and regularization

3. **Repeat** for multiple cycles (default: 3 cycles)

### Loss Functions

**VAE Loss:**
```
Total Loss = Reconstruction Loss + β × KL Divergence
```

**Reconstruction Network Loss:**
```
MSE Loss = Mean Squared Error between refined output and original data
```

## Complete Example Workflow

### Default Training (KIRC)

```bash
# Train with default parameters
python scripts/pipeline_steps/2_models.py
```

### Custom Training (BRCA)

Edit the script:

```python
# Change cancer type
cancer_type = "BRCA"

# Update paths
path_rnaseq = "data/interim/preprocessed_BRCA_data/preprocessed_rnaseq.csv"
path_clinical = "data/interim/preprocessed_BRCA_data/clinical_data.csv"

# Adjust VAE architecture for different gene count
vae_config.INPUT_DIM = X_train.shape[1]  # Auto-detects gene count
vae_config.MID_DIM = 1024  # Larger middle layer for BRCA
vae_config.LATENT_DIM = 512  # Larger latent space

# Adjust reconstruction network
recnet_dims = [X_train.shape[1], 4096, 1024, 4096, X_train.shape[1]]
```

## Monitoring Training

### Training Progress

The script logs training progress:

```
[INFO] Epoch 1/600: Loss=1234.56, Recon=1200.00, KL=34.56, β=0.01
[INFO] Epoch 50/600: Loss=456.78, Recon=400.00, KL=56.78, β=0.50
[INFO] Epoch 100/600: Loss=234.56, Recon=180.00, KL=54.56, β=1.00
...
```

### Training Plots

After training, inspect the generated plots:

1. **VAE Training History**: Shows total loss, reconstruction loss, and KL divergence over epochs
2. **Reconstruction Network History**: Shows training and test MSE loss

## Verification

After training, verify the models:

```python
import torch
import pandas as pd
from renalprog.modeling.train import VAE

# Load trained VAE
model_dir = "models/YYYYMMDD_models_KIRC"
checkpoint = torch.load(f"{model_dir}/vae/final_model.pth")

# Load test data
X_test = pd.read_csv("data/interim/YYYYMMDD_train_test_split/X_test.csv", index_col=0)

# Test reconstruction
# ... (see Step 3 for full reconstruction testing)
```

## Next Steps

After completing VAE training:

1. **Verify models**: Check training plots and final losses
2. **Test reconstruction**: Evaluate on test set
3. **Proceed to Step 3**: Check reconstruction quality and generate synthetic data

```bash
python scripts/pipeline_steps/3_check_reconstruction.py
```

## Additional Resources

- [VAE Tutorial](https://arxiv.org/abs/1312.6114) - Original VAE paper
- [β-VAE](https://openreview.net/forum?id=Sy2fzU9gl) - Disentangled representations
- [PyTorch VAE Examples](https://github.com/pytorch/examples/tree/master/vae)
- [Understanding VAE Loss](https://stats.stackexchange.com/questions/332179/how-to-weight-kld-loss-vs-reconstruction-loss-in-variational-auto-encoder)
