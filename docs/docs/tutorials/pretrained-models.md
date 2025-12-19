# Using Pretrained Models

This guide shows you how to use the pretrained VAE models provided with the RenalProg package to reproduce the results from the paper.

## Overview

The repository includes pretrained models for two cancer types:

- **KIRC** (Kidney Renal Clear Cell Carcinoma)
- **BRCA** (Breast Invasive Carcinoma)

Each pretrained model comes with:

- Trained VAE weights
- Reconstruction network weights
- Network architecture specifications

## Pretrained Model Structure

```
models/pretrained/
├── KIRC/
│   ├── 20250321_VAE_idim8516_md512_feat256mse_relu.pth  # VAE weights
│   ├── network_reconstruction.pth                       # Reconstruction network weights
│   └── network_dims.csv                                 # Architecture specifications
└── BRCA/
    ├── 20251209_VAE_idim8954_md1024_feat512mse_relu.pth
    ├── network_reconstruction.pth
    └── network_dims.csv
```

## Preprocessed Data

The preprocessed data for each cancer type is available in:

```
data/interim/
├── preprocessed_KIRC/
│   ├── clinical.csv    # Patient metadata (stages, demographics)
│   └── rnaseq.csv      # Preprocessed RNA-seq expression data
└── preprocessed_BRCA/
    ├── clinical.csv
    └── rnaseq.csv
```

### Data Format

**clinical.csv**: Patient metadata with columns:
- Index: Patient IDs (e.g., `TCGA-A3-3319-01`)
- `ajcc_pathologic_tumor_stage`: Cancer stage (Stage I, II, III, IV)
- Additional clinical features (varies by cancer type)

**rnaseq.csv**: Gene expression data
- Index: Gene symbols (e.g., `TSPAN6`, `TNMD`)
- Columns: Patient IDs
- Values: Log-transformed normalized gene expression

## Quick Start

### 1. Generate Trajectories from Pretrained Models

Use the provided script to generate patient trajectories:

```bash
# For KIRC
python scripts/pipeline_steps/use_pretrained_model.py \
    --cancer_type KIRC \
    --model_dir models/pretrained/KIRC \
    --data_dir data/interim/preprocessed_KIRC \
    --output_dir data/processed/trajectories_KIRC_pretrained

# For BRCA
python scripts/pipeline_steps/use_pretrained_model.py \
    --cancer_type BRCA \
    --model_dir models/pretrained/BRCA \
    --data_dir data/interim/preprocessed_BRCA \
    --output_dir data/processed/trajectories_BRCA_pretrained
```

### 2. Run Enrichment Analysis

After generating trajectories, run pathway enrichment analysis:

```bash
# For KIRC
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --trajectory_dir data/processed/trajectories_KIRC_pretrained/early_to_late/test_to_test \
    --output_dir data/processed/enrichment_KIRC_pretrained \
    --n_threads 8 \
    --gsea_path ./GSEA_4.3.2/gsea-cli.sh \
    --pathways_file data/external/ReactomePathways.gmt

# For BRCA
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --trajectory_dir data/processed/trajectories_BRCA_pretrained/early_to_late/test_to_test \
    --output_dir data/processed/enrichment_BRCA_pretrained \
    --n_threads 8 \
    --gsea_path ./GSEA_4.3.2/gsea-cli.sh \
    --pathways_file data/external/ReactomePathways.gmt
```

### 3. Generate Pathway Heatmaps

Create visualizations of pathway enrichment:

```bash
# For KIRC
python scripts/pipeline_steps/6b_generate_pathway_heatmap.py \
    --enrichment_file data/processed/enrichment_KIRC_pretrained/trajectory_enrichment.csv \
    --output_dir data/processed/enrichment_KIRC_pretrained \
    --fdr_threshold 0.05

# For BRCA
python scripts/pipeline_steps/6b_generate_pathway_heatmap.py \
    --enrichment_file data/processed/enrichment_BRCA_pretrained/trajectory_enrichment.csv \
    --output_dir data/processed/enrichment_BRCA_pretrained \
    --fdr_threshold 0.05
```

## Advanced Usage

### Loading Pretrained Models in Python

You can load and use the pretrained models directly in your Python code:

```python
import torch
import pandas as pd
from pathlib import Path
from renalprog.modeling.train import VAE, NetworkReconstruction
from renalprog.utils import get_device

# Configuration
cancer_type = "KIRC"
model_dir = Path("models/pretrained/KIRC")
data_dir = Path("data/interim/preprocessed_KIRC")

# Load network architecture
network_dims = pd.read_csv(model_dir / "network_dims.csv", index_col=0)
input_dim = int(network_dims['in_dim'].values[0])
layer1_dim = int(network_dims['layer1_dim'].values[0])
layer2_dim = int(network_dims['layer2_dim'].values[0])
layer3_dim = int(network_dims['layer3_dim'].values[0])

# Initialize VAE model
# For KIRC: mid_dim=512, latent_dim=256
# For BRCA: mid_dim=1024, latent_dim=512
device = get_device(force_cpu=False)

if cancer_type == "KIRC":
    vae_model = VAE(
        input_dim=input_dim,
        mid_dim=512,
        latent_dim=256,
        loss_fn='mse',
        activation='relu'
    ).to(device)
else:  # BRCA
    vae_model = VAE(
        input_dim=input_dim,
        mid_dim=1024,
        latent_dim=512,
        loss_fn='mse',
        activation='relu'
    ).to(device)

# Load VAE weights
vae_weights = list(model_dir.glob("*VAE_*.pth"))[0]
vae_model.load_state_dict(torch.load(vae_weights, map_location=device))
vae_model.eval()

# Initialize reconstruction network
reconstruction_network = NetworkReconstruction(
    layer_dims=[input_dim, layer1_dim, layer2_dim, layer3_dim, input_dim]
).to(device)

# Load reconstruction network weights
reconstruction_network.load_state_dict(
    torch.load(model_dir / "network_reconstruction.pth", map_location=device)
)
reconstruction_network.eval()

# Load data
rnaseq = pd.read_csv(data_dir / "rnaseq.csv", index_col=0)
clinical = pd.read_csv(data_dir / "clinical.csv", index_col=0)

# Generate latent representations
with torch.no_grad():
    data_tensor = torch.FloatTensor(rnaseq.T.values).to(device)
    _, mu, logvar = vae_model(data_tensor)
    latent_reps = mu.cpu().numpy()

print(f"Generated latent representations: {latent_reps.shape}")
```

### Generating Custom Trajectories

You can customize trajectory generation parameters:

```python
from renalprog.modeling.predict import (
    calculate_all_possible_transitions,
    link_patients_random,
    build_trajectory_network,
    generate_trajectory_data
)

# Calculate transitions
all_traj = calculate_all_possible_transitions(
    data=rnaseq,
    metadata_selection=clinical[['stage']],
    start_stage='early',
    end_stage='late',
    link_next=5,  # Number of nearest neighbors
    distance_metric='wasserstein'
)

# Build trajectory network
trajectory_network = build_trajectory_network(
    all_traj=all_traj,
    seed=2023
)

# Generate synthetic trajectories
generate_trajectory_data(
    data=rnaseq,
    metadata=clinical,
    trajectory_network=trajectory_network,
    vae_model=vae_model,
    reconstruction_network=reconstruction_network,
    n_timepoints=50,
    output_dir=output_dir / "trajectories",
    interpolation_method='linear',
    device=device
)
```

## Model Specifications

### KIRC Pretrained Model

- **VAE Architecture**:
  - Input dimension: 8,516 genes
  - Middle dimension: 512
  - Latent dimension: 256
  - Loss function: MSE
  - Activation: ReLU

- **Reconstruction Network**:
  - Architecture: [8,954, 3,104, 790, 4,027, 8,954]
  
- **Training Details**:
  - Beta-VAE with 3 cycles
  - 600 total epochs
  - Batch size: 8

### BRCA Pretrained Model

- **VAE Architecture**:
  - Input dimension: 8,954 genes
  - Middle dimension: 1,024
  - Latent dimension: 512
  - Loss function: MSE
  - Activation: ReLU

- **Reconstruction Network**:
  - Architecture: [8,954, 3,104, 790, 4,027, 8,954]

- **Training Details**:
  - Beta-VAE with 3 cycles
  - 600 total epochs
  - Batch size: 8

## Reproducing Paper Results

To reproduce the exact results from the paper:

1. **Use the pretrained models** (recommended for exact reproduction)
2. **Follow the complete pipeline**:

```bash
# 1. Generate trajectories
python scripts/pipeline_steps/use_pretrained_model.py \
    --cancer_type KIRC \
    --model_dir models/pretrained/KIRC \
    --data_dir data/interim/preprocessed_KIRC \
    --output_dir data/processed/paper_reproduction_KIRC

# 2. Run enrichment analysis
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --trajectory_dir data/processed/paper_reproduction_KIRC/early_to_late/test_to_test \
    --output_dir data/processed/enrichment_paper_KIRC \
    --n_threads 8 \
    --gsea_path ./GSEA_4.3.2/gsea-cli.sh \
    --pathways_file data/external/ReactomePathways.gmt

# 3. Generate heatmaps
python scripts/pipeline_steps/6b_generate_pathway_heatmap.py \
    --enrichment_file data/processed/enrichment_paper_KIRC/trajectory_enrichment.csv \
    --output_dir data/processed/enrichment_paper_KIRC \
    --fdr_threshold 0.05
```

## Troubleshooting

### GPU Memory Issues

If you encounter GPU memory issues, force CPU usage:

```python
device = get_device(force_cpu=True)
```

### Missing Dependencies

Ensure all dependencies are installed:

```bash
pip install -r requirements.txt
```

### Data Compatibility

Ensure your input data matches the expected format:
- Gene expression data should be log-transformed
- Clinical data should include `ajcc_pathologic_tumor_stage` column
- Patient IDs should match between rnaseq.csv and clinical.csv

## Next Steps

- [Complete Pipeline Tutorial](complete-pipeline.md)
- [Enrichment Analysis Guide](../ENRICHMENT_ANALYSIS.md)
- [API Reference](../api/index.md)

## Citation

If you use these pretrained models, please cite:

```bibtex
@article{renalprog2024,
  title={Deep Learning Analysis of Cancer Progression Trajectories},
  author={Your Name et al.},
  journal={In Preparation},
  year={2024},
  note={Preprint available soon}
}
```

