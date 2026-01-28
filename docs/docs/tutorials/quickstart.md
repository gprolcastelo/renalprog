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

See the [Installation Guide](installation.md) for detailed instructions.

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

You can either download preprocessed data from Zenodo or use your own data.

#### Option A: Download Preprocessed Data from Zenodo (Recommended)

```python
from pathlib import Path

# Download preprocessed KIRC data from Zenodo
rnaseq, clinical = dataset.download_preprocessed_from_zenodo(
    rnaseq_url='https://zenodo.org/records/17987300/files/Static_KIRC.csv?download=1',
    clinical_url='https://zenodo.org/records/17987300/files/nodes_metadata.csv?download=1',
    output_dir=Path('data/interim/preprocessed_KIRC_data')
)

print(f"RNA-seq data shape: {rnaseq.shape}")
print(f"Genes: {rnaseq.shape[0]:,}, Samples: {rnaseq.shape[1]:,}")
print(f"Clinical data shape: {clinical.shape}")
```

Expected output:
```
================================================================================
Downloading preprocessed data from Zenodo
================================================================================

[DOWNLOADING RNASEQ DATA]
  Source URL: https://zenodo.org/records/17987300/files/Static_KIRC.csv?download=1
  Destination: data/interim/preprocessed_KIRC_data/preprocessed_rnaseq.csv
Downloading RNAseq: 100.0%
  ✓ Successfully downloaded RNAseq data
  Shape: 8,516 genes × 530 samples

[DOWNLOADING CLINICAL DATA]
  Source URL: https://zenodo.org/records/17987300/files/nodes_metadata.csv?download=1
  Destination: data/interim/preprocessed_KIRC_data/clinical_data.csv
Downloading clinical: 100.0%
  ✓ Successfully downloaded clinical data
  Shape: 530 samples × 4 features

[DATA VALIDATION]
  RNAseq samples: 530
  Clinical samples: 530
  Common samples: 530

================================================================================
Preprocessed data download completed successfully
================================================================================
```

!!! tip "Available Datasets"
    Preprocessed data is available for:
    
    - **KIRC** (Kidney Renal Clear Cell Carcinoma): [Zenodo Record 17987300](https://zenodo.org/records/17987300)
    - **BRCA** (Breast Cancer): [Zenodo Record 17986123](https://zenodo.org/records/17986123)

#### Option B: Load Your Own Preprocessed Data

```python
# Load from local files
rnaseq = pd.read_csv('data/processed/rnaseq_maha.csv', index_col=0)
clinical = pd.read_csv('data/processed/clinical.csv', index_col=0)

print(f"RNA-seq data shape: {rnaseq.shape}")
print(f"Samples: {rnaseq.shape[0]}, Genes: {rnaseq.shape[1]}")
```

### 3. Train VAE

First, create train/test splits:

```python
from renalprog.config import VAEConfig, INTERIM_DATA_DIR, MODELS_DIR, get_dated_dir
from renalprog.modeling.train import train_vae_with_postprocessing
from pathlib import Path

# Create train/test splits (80/20)
traintest_dir = INTERIM_DATA_DIR / "train_test_split"
X_train, X_test, y_train, y_test, _, _ = dataset.create_train_test_split(
    rnaseq_path=Path('data/interim/preprocessed_KIRC_data/preprocessed_rnaseq.csv'),
    clinical_path=Path('data/interim/preprocessed_KIRC_data/clinical_data.csv'),
    test_size=0.2,
    seed=2023,
    output_dir=traintest_dir,
    stage_column="stage"
)

print(f"Training samples: {X_train.shape[0]}")
print(f"Test samples: {X_test.shape[0]}")
print(f"Number of genes: {X_train.shape[1]}")
```

Now configure and train the VAE: the given parameters are intended for the quick start testing. 
Increase the number of epochs for better performance in real applications.

```python
# Configure VAE
vae_config = VAEConfig()
vae_config.INPUT_DIM = X_train.shape[1]  # Number of genes
vae_config.MID_DIM = 128
vae_config.LATENT_DIM = 16
vae_config.BETA_CYCLES = 1
vae_config.EPOCHS = 100 * vae_config.BETA_CYCLES  # 600 epochs total
vae_config.BETA_RATIO = 0.5
vae_config.BATCH_SIZE = 8

# Configure postprocessing Reconstruction Network
recnet_dims = [X_train.shape[1], 3_512, 824, 3_731, X_train.shape[1]]

# Create output directory
model_dir = MODELS_DIR / "quickstart_vae"
model_dir.mkdir(parents=True, exist_ok=True)

# Train VAE with postprocessing network
vae_model, network, vae_history, reconstruction_history = train_vae_with_postprocessing(
    X_train=X_train,
    X_test=X_test,
    vae_config=vae_config,
    reconstruction_network_dims=recnet_dims,
    reconstruction_epochs=1_000,
    reconstruction_lr=1e-4,
    batch_size_reconstruction=8,
    save_dir=model_dir,
    force_cpu=True,  # Set to False to use GPU if available
)

print(f"Training complete! Models saved to {model_dir}")
```

**Expected output:**
```
[INFO] VAE Configuration:
[INFO]   Input dim: 8516
[INFO]   Mid dim: 128
[INFO]   Latent dim: 16
[INFO]   Beta cycles: 1
[INFO]   Total epochs: 100
[INFO]   Batch size: 8
[INFO] Reconstruction Network dims: [8516, 3512, 824, 3731, 8516]
[INFO] Models will be saved to: models/quickstart_vae
[INFO] 
Starting training...
[INFO] Starting full VAE + postprocessing pipeline
[INFO] Step 1: Training VAE
[INFO] Saved config: models/quickstart_vae/vae/config.json
[INFO] Using device: cpu
[INFO] Model: VAE(input_dim=8516, mid_dim=128, latent_dim=16)
[INFO] Parameters: 2,195,044
[INFO] ModelCheckpointer initialized: models/quickstart_vae/vae
[INFO] Monitoring: val_loss (min)
[INFO] Using cyclical beta annealing: 0.0 -> 1.0 over 1 cycles
[INFO] Starting training for 100 epochs
Epochs:   0%|                                                                                                                                         | 0/100 [00:00<?, ?it/s, train_loss=1667.9040, val_loss=1361.7740, beta=0.000][INFO] Epoch 1/100 - train_loss: 1667.9040, val_loss: 1361.7740                                                                                                                                                                     
Epochs:   9%|███████████▊                                                                                                                       | 9/100 [00:09<01:24,  1.08it/s, train_loss=552.9183, val_loss=710.7239, beta=0.180][INFO] Epoch 10/100 - train_loss: 552.9183, val_loss: 710.7239                                                                                                                                                                      
Epochs:  19%|████████████████████████▋                                                                                                         | 19/100 [00:22<02:07,  1.58s/it, train_loss=556.0336, val_loss=709.1441, beta=0.380][INFO] Epoch 20/100 - train_loss: 556.0336, val_loss: 709.1441                                                                                                                                                                      
Epochs:  29%|█████████████████████████████████████▋                                                                                            | 29/100 [00:37<01:43,  1.46s/it, train_loss=526.8585, val_loss=669.9171, beta=0.580][INFO] Epoch 30/100 - train_loss: 526.8585, val_loss: 669.9171                                                                                                                                                                      
Epochs:  39%|██████████████████████████████████████████████████▋                                                                               | 39/100 [00:52<01:30,  1.48s/it, train_loss=462.7390, val_loss=596.3372, beta=0.780][INFO] Epoch 40/100 - train_loss: 462.7390, val_loss: 596.3372                                                                                                                                                                      
Epochs:  49%|███████████████████████████████████████████████████████████████▋                                                                  | 49/100 [01:06<01:10,  1.38s/it, train_loss=389.8898, val_loss=518.6038, beta=0.980][INFO] Epoch 50/100 - train_loss: 389.8898, val_loss: 518.6038                                                                                                                                                                      
Epochs:  59%|████████████████████████████████████████████████████████████████████████████▋                                                     | 59/100 [01:19<00:53,  1.30s/it, train_loss=335.6482, val_loss=460.9847, beta=1.000][INFO] Epoch 60/100 - train_loss: 335.6482, val_loss: 460.9847                                                                                                                                                                      
Epochs:  69%|█████████████████████████████████████████████████████████████████████████████████████████▋                                        | 69/100 [01:31<00:40,  1.31s/it, train_loss=307.4428, val_loss=433.4705, beta=1.000][INFO] Epoch 70/100 - train_loss: 307.4428, val_loss: 433.4705                                                                                                                                                                      
Epochs:  79%|██████████████████████████████████████████████████████████████████████████████████████████████████████▋                           | 79/100 [01:44<00:25,  1.22s/it, train_loss=293.1740, val_loss=419.2723, beta=1.000][INFO] Epoch 80/100 - train_loss: 293.1740, val_loss: 419.2723                                                                                                                                                                      
Epochs:  89%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████▋              | 89/100 [01:57<00:14,  1.29s/it, train_loss=278.1044, val_loss=404.8974, beta=1.000][INFO] Epoch 90/100 - train_loss: 278.1044, val_loss: 404.8974                                                                                                                                                                      
Epochs:  99%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▋ | 99/100 [02:10<00:01,  1.28s/it, train_loss=261.5659, val_loss=388.2671, beta=1.000][INFO] Epoch 100/100 - train_loss: 261.5659, val_loss: 388.2671                                                                                                                                                                     
Epochs: 100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 100/100 [02:10<00:00,  1.31s/it, train_loss=261.5659, val_loss=388.2671, beta=1.000]
[INFO] Saved checkpoint: models/quickstart_vae/vae/final_model.pth
[INFO] Step 3: Training reconstruction network
[INFO] Training reconstruction network for 10 epochs
Reconstruction Network Training: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 10/10 [05:38<00:00, 33.87s/it, train_loss=53.913388, test_loss=54.212544]
[INFO] Reconstruction network training complete        
```

!!! tip "Training Time"
    With these settings (100 VAE epochs + 10 reconstruction epochs, and the rest of the configuration parameters), CPU training should take 6-7 minutes. Using a GPU will significantly speed up training.

!!! info "Two-Stage Training"
    This pipeline trains:
    
    1. **VAE**: Learns compressed latent representation
    2. **Reconstruction Network**: Post-processes VAE output for better reconstruction
    
    Both models are saved and can be used for downstream analyses. The pipeline avoids data leakage by training moth models on training data only.

### 4. Visualize Training Progress

```python
from renalprog.plots import plot_training_history, plot_reconstruction_losses

# Plot VAE training history
plot_training_history(
    history=vae_history,
    save_path=model_dir / "vae_training_history.png",
    title="VAE Training History"
)

# Plot reconstruction network training history
plot_reconstruction_losses(
    loss_train=reconstruction_history["train_loss"],
    loss_test=reconstruction_history["test_loss"],
    save_path=model_dir / "reconstruction_training_history.png",
    title="Reconstruction Network Training History"
)

print("Training plots saved!")
```

## Full quick start code

You can find the complete code for this quick start example in the `scripts/quickstart_vae.py` script in the repository.
Easily run it from the command line:

```bash
python scripts/quickstart_vae.py
``` 


## Next Steps

Congratulations! You've successfully:

✅ Downloaded preprocessed data from Zenodo  
✅ Split data into train/test sets  
✅ Trained a VAE with postprocessing network  
✅ Visualized training progress  

Now that you've completed the quick start, explore the full tutorials:

### 📚 Recommended Tutorials

1. **[Complete Pipeline](complete-pipeline.md)** - Full analysis workflow from data processing to results
2. **[Step 3: VAE Reconstruction Check](step3-reconstruction.md)** - Evaluate VAE reconstruction quality
3. **[Step 4: Generate Trajectories](step4-trajectories.md)** - Create disease progression trajectories
4. **[Step 5: Classification](step5-classification.md)** - Train classifiers on gene expression
5. **[Step 6: Enrichment Analysis](step6-enrichment.md)** - Perform pathway enrichment analysis

## 📜 Citation

If you use this pipeline in your research, please cite this work.
For details, see the [citation page](../citation.md) file.