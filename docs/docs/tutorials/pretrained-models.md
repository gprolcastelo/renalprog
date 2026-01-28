# Using Pretrained Models

This guide shows you how to use the pretrained VAE models from Hugging Face to reproduce the results from the paper without training models from scratch.
The pretrained models were trained on the preprocessed TCGA data available in this repository (see [Preprocessing Tutorial](step1-data-processing.md) for details).

## Overview

Pre-trained models are available on Hugging Face Hub and include:

- **VAE models**: Variational Autoencoders trained on TCGA data
- **Reconstruction Networks**: Post-processing networks that refine VAE outputs
- **Configuration files**: Model architectures and hyperparameters

!!! info "Available Cancer Types"
    Pretrained models are available for:
    
    - **KIRC** (Kidney Renal Clear Cell Carcinoma)
    - **BRCA** (Breast Invasive Carcinoma)

## Quick Start

### Using the Pipeline Script

The easiest way to use pretrained models is with the `3_check_reconstruction.py` script:

```bash
# Download and use KIRC pretrained models
python scripts/pipeline_steps/3_check_reconstruction.py \
    --hf_models \
    --cancer_type KIRC

# Or for BRCA
python scripts/pipeline_steps/3_check_reconstruction.py \
    --hf_models \
    --cancer_type BRCA
```

This will:
1. Download the VAE model and configuration from Hugging Face
2. Download the Reconstruction Network model
3. Load your preprocessed data
4. Generate reconstructions
5. Create UMAP visualizations comparing original vs reconstructed data

## Manual Usage in Python

### Step 1: Install Hugging Face Hub

Check the [official documentation](https://huggingface.co/docs/huggingface_hub/installation) to install the `huggingface-hub` library.
The simplest way is via pip:


```bash
pip install huggingface-hub
```

### Step 2: Download and Load Models

```python
import huggingface_hub as hf
import torch
import json
from pathlib import Path
from renalprog.modeling.train import VAE, NetworkReconstruction
from renalprog.config import MODELS_DIR

# Set cancer type
cancer_type = 'KIRC'  # or 'BRCA'

# Create local directory for pretrained models
model_dir = MODELS_DIR / "pretrained" / cancer_type
model_dir.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Download VAE Configuration
# ============================================================================
print(f"Downloading VAE config for {cancer_type}...")
vae_config_path = hf.hf_hub_download(
    repo_id="gprolcastelo/evenflow_models",
    filename=f"{cancer_type}/config.json",
    local_dir=model_dir.parent
)

# Load configuration
with open(vae_config_path, 'r') as f:
    vae_config = json.load(f)

print(f"VAE Configuration: {vae_config}")

# ============================================================================
# Download and Load VAE Model
# ============================================================================
# Model filenames for each cancer type
vae_models = {
    'KIRC': "KIRC/20250321_VAE_idim8516_md512_feat256mse_relu.pth",
    'BRCA': "BRCA/20251209_VAE_idim8954_md1024_feat512mse_relu.pth"
}

print(f"Downloading VAE model for {cancer_type}...")
vae_model_path = hf.hf_hub_download(
    repo_id="gprolcastelo/evenflow_models",
    filename=vae_models[cancer_type],
    local_dir=model_dir.parent
)

# Initialize VAE
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model_vae = VAE(
    input_dim=vae_config['INPUT_DIM'],
    mid_dim=vae_config['MID_DIM'],
    features=vae_config['LATENT_DIM']
).to(device)

# Load weights
checkpoint = torch.load(vae_model_path, map_location=device, weights_only=False)
model_vae.load_state_dict(checkpoint)
model_vae.eval()

print(f"✓ VAE model loaded successfully!")

# ============================================================================
# Download and Load Reconstruction Network
# ============================================================================
print(f"Downloading Reconstruction Network for {cancer_type}...")

# Download network dimensions
network_dims_path = hf.hf_hub_download(
    repo_id="gprolcastelo/evenflow_models",
    filename=f"{cancer_type}/network_dims.csv",
    local_dir=model_dir.parent
)

# Load dimensions
import pandas as pd
network_dims = pd.read_csv(network_dims_path).values.tolist()[0]
print(f"Network dimensions: {network_dims}")

# Download model
recnet_model_path = hf.hf_hub_download(
    repo_id="gprolcastelo/evenflow_models",
    filename=f"{cancer_type}/network_reconstruction.pth",
    local_dir=model_dir.parent
)

# Initialize Reconstruction Network
model_recnet = NetworkReconstruction(layer_dims=network_dims).to(device)

# Load weights
checkpoint_recnet = torch.load(recnet_model_path, map_location=device, weights_only=False)
model_recnet.load_state_dict(checkpoint_recnet)
model_recnet.eval()

print(f"✓ Reconstruction Network loaded successfully!")
```

### Step 3: Use Models for Inference

```python
from renalprog.utils import apply_VAE
import pandas as pd

# Load your preprocessed data
data = pd.read_csv('data/interim/preprocessed_KIRC_data/preprocessed_rnaseq.csv', index_col=0)

# If data is genes × samples, transpose it to samples × genes
if data.shape[0] > data.shape[1]:
    data = data.T

print(f"Data shape: {data.shape}")

# ============================================================================
# Apply VAE
# ============================================================================
data_tensor = torch.tensor(data.values, dtype=torch.float32)

reconstruction_vae, _, _, latent, scaler = apply_VAE(
    data_tensor, 
    model_vae, 
    y=None
)

print(f"VAE reconstruction shape: {reconstruction_vae.shape}")
print(f"Latent representation shape: {latent.shape}")

# Convert to DataFrame
df_reconstruction_vae = pd.DataFrame(
    reconstruction_vae, 
    index=data.index, 
    columns=data.columns
)

df_latent = pd.DataFrame(
    latent, 
    index=data.index
)

# ============================================================================
# Apply Reconstruction Network (Post-processing)
# ============================================================================
rec_tensor = torch.tensor(reconstruction_vae, dtype=torch.float32).to(device)

with torch.no_grad():
    reconstruction_final = model_recnet(rec_tensor)

# Convert to DataFrame
df_reconstruction_final = pd.DataFrame(
    reconstruction_final.cpu().numpy(),
    index=data.index,
    columns=data.columns
)

print(f"Final reconstruction shape: {df_reconstruction_final.shape}")

# Save results
df_reconstruction_final.to_csv('reconstructed_data.csv')
df_latent.to_csv('latent_representation.csv')

print("✓ Reconstruction complete!")
```

## VAE Model Architecture Details

| Model | Input Dim | Mid Dim | Latent Dim | File                                                    |
|-------|-----------|---------|------------|---------------------------------------------------------|
| KIRC  | 8,516 | 512 | 256 | `KIRC/20250321_VAE_idim8516_md512_feat256mse_relu.pth`  |
| BRCA  | 8,954 | 1,024 | 512 | `BRCA/20251209_VAE_idim8954_md1024_feat512mse_relu.pth` |

## Hugging Face Repository

All pretrained models are hosted at:

🤗 **[gprolcastelo/evenflow_models](https://huggingface.co/gprolcastelo/evenflow_models)**

### Repository Structure

```
evenflow_models/
├── KIRC/
│   ├── config.json
│   ├── network_dims.csv
│   ├── network_reconstruction.pth
│   └── 20250321_VAE_idim8516_md512_feat256mse_relu.pth
└── BRCA/
    ├── config.json
    ├── network_dims.csv
    ├── network_reconstruction.pth
    └── 20251209_VAE_idim8954_md1024_feat512mse_relu.pth
```

## Complete Example: Reconstruction Validation

For a detailed example of using the pretrained models to validate reconstructions, refer to the [reconstruction tutorial](step3-reconstruction.md).
As a summary, when using `3_check_reconstruction.py`, you have several options:

```bash
# Basic usage with pretrained models
python scripts/pipeline_steps/3_check_reconstruction.py --hf_models --cancer_type KIRC

# Include SDMetrics evaluation (takes longer)
python scripts/pipeline_steps/3_check_reconstruction.py --hf_models --cancer_type KIRC --sdmetrics

# Use locally trained models instead
python scripts/pipeline_steps/3_check_reconstruction.py --cancer_type KIRC
```

### Arguments

| Argument | Description                              | Default |
|----------|------------------------------------------|---------|
| `--cancer_type` | Cancer type (KIRC or BRCA)               | KIRC |
| `--hf_models` | Load pretrained models from Hugging Face | False |
| `--sdmetrics` | Evaluate using SDMetrics (very slow)     | False |

## Output Files

When running the reconstruction check, you'll get:

```
reports/figures/YYYYMMDD_CANCER_umap_reconstruction/
├── preprocessed.html                          # UMAP of original data
├── VAE_output.html                            # UMAP of VAE reconstruction
├── recnet_output.html                         # UMAP of final reconstruction
├── preprocessed_and_vae.html                  # Comparison: original vs VAE
└── preprocessed_and_recnet.html               # Comparison: original vs final

models/pretrained/CANCER/
├── config.json                                # VAE configuration
├── network_dims.csv                           # Reconstruction network architecture
├── CANCER/
│   ├── 20250321_VAE_*.pth                    # VAE weights
│   └── network_reconstruction.pth             # Reconstruction network weights
```

## 📜 Citation

If you use these pretrained models in your research, please [cite](../citation.md)

## See Also

- [Reconstruction Tutorial](step3-reconstruction.md)
- [VAE Training Tutorial](step2-vae-training.md)
- [Complete Pipeline](complete-pipeline.md)
- [API Documentation - Models](../api/models.md)

