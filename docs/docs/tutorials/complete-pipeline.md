# Complete Pipeline Tutorial

This tutorial walks through the complete `renalprog` analysis pipeline from raw TCGA data to pathway enrichment results. This represents the full workflow used in the publication.

## Overview

The pipeline consists of six main steps:

1. **Data Processing**: Download and preprocess TCGA data
2. **VAE Training**: Train deep generative models  
3. **Reconstruction Validation**: Assess model quality
4. **Trajectory Generation**: Create synthetic progression paths
5. **Classification**: Stage prediction and biomarker discovery
6. **Enrichment Analysis**: Pathway-level interpretation

!!! note
    We recommend using HPC for steps 2 and 6 due to computational demands.
    Step 2 will benefit from GPU acceleration, while step 6 can be parallelized across multiple CPU cores.

## Prerequisites

Ensure you have completed the installation: see [Installation Guide](installation.md).

## Pipeline Execution

### Option 1: Run All Steps Sequentially

```bash
# Navigate to scripts directory
cd scripts/pipeline_steps

# Step 1: Data Processing (30 minutes)
python 1_data_processing.py

# Step 2: VAE Training (2-4 hours with GPU)
python 2_models.py

# Step 3: Reconstruction Check (10 minutes)
python 3_check_reconstruction.py

# Step 4: Trajectory Generation (30 minutes)
python 4_trajectories.py

# Step 5: Classification (1 hour)
python 5_classification.py
```

### Option 2: Run with Makefile

```bash
# Run entire pipeline
make pipeline

# Or run individual steps
make data
make models
make trajectories
make classification
```



## Detailed Step-by-Step Guide

### Step 1: Data Processing

**Script**: `scripts/pipeline_steps/1_data_processing.py`

**What it does**:
- Downloads TCGA KIRC RNA-seq and clinical data
- Filters low-expression genes
- Removes outliers using Mahalanobis distance
- Creates train/test splits

**Expected outputs**:
```
data/interim/preprocessed_KIRC_data/
├── preprocessed_rnaseq.csv         # Filtered gene expression (498 samples × ~5000 genes)
├── preprocessing_info.json         # Filter statistics
└── stages.csv                      # Cancer stage labels
```

**Validation**:
```python
import pandas as pd

rnaseq = pd.read_csv('data/interim/preprocessed_KIRC_data/preprocessed_rnaseq.csv', index_col=0)
print(f"Samples: {rnaseq.shape[0]}, Genes: {rnaseq.shape[1]}")
# Expected: Samples: 498, Genes: ~5000

# Check for missing values
assert rnaseq.isnull().sum().sum() == 0, "Missing values detected!"
print("✓ Data quality check passed")
```

---

### Step 2: VAE Training

**Script**: `scripts/pipeline_steps/2_models.py`

**What it does**:
- Trains Variational Autoencoder on gene expression
- Learns low-dimensional latent representation
- Trains reconstruction network for improved decoding
- Saves models and training histories

**Configuration** (in script):
```python
vae_config.INPUT_DIM = X_train.shape[1]  # Number of genes
vae_config.MID_DIM = 512                  # Hidden layer size
vae_config.LATENT_DIM = 256               # Latent space dimensionality
vae_config.EPOCHS = 600                   # Total epochs (3 cycles × 200)
vae_config.BATCH_SIZE = 8
```

**Expected outputs**:
```
models/20251217_models_KIRC/
├── vae/
│   ├── vae_model.pt                 # Trained VAE
│   ├── vae_config.json              # Model configuration
│   └── vae_training_history.png     # Loss curves
├── reconstruction/
│   ├── reconstruction_network.pt    # Post-processing network
│   └── reconstruction_network_history.png
└── training_data/
    ├── X_train.csv
    └── X_test.csv
```

**Validation**:
```python
import torch
from renalprog.modeling.models import VAE

# Load model
vae = VAE.load('models/20251217_models_KIRC/vae/vae_model.pt')

# Check reconstruction
X_test = pd.read_csv('models/20251217_models_KIRC/training_data/X_test.csv', index_col=0)
X_tensor = torch.FloatTensor(X_test.values)

with torch.no_grad():
    X_recon = vae(X_tensor)[0]
    mse = ((X_tensor - X_recon) ** 2).mean()
    print(f"Reconstruction MSE: {mse:.4f}")
    # Expected: MSE < 1.0 for well-trained model
```

---

### Step 3: Reconstruction Validation

**Script**: `scripts/pipeline_steps/3_check_reconstruction.py`

**What it does**:
- Visualizes latent space with UMAP/t-SNE
- Compares original vs. reconstructed gene expression
- Generates quality control plots

**Expected outputs**:
```
reports/figures/reconstruction/
├── latent_space_umap.png           # UMAP colored by stage
├── latent_space_tsne.png           # t-SNE colored by stage
├── reconstruction_scatter.png      # Original vs. reconstructed
└── gene_correlation.png            # Per-gene reconstruction quality
```

**Validation**:
```python
# Check that plots exist
import os

plot_dir = 'reports/figures/reconstruction/'
required_plots = [
    'latent_space_umap.png',
    'reconstruction_scatter.png'
]

for plot in required_plots:
    assert os.path.exists(os.path.join(plot_dir, plot)), f"Missing: {plot}"
print("✓ All validation plots generated")
```

---

### Step 4: Trajectory Generation

**Script**: `scripts/pipeline_steps/4_trajectories.py`

**What it does**:
- Generates synthetic patient trajectories
- Interpolates between early and late stage samples
- Decodes trajectories to gene expression space
- Creates multiple trajectory types

**Configuration**:
```python
n_trajectories = 500                # Number of trajectories to generate
n_timepoints = 20                   # Timepoints per trajectory
trajectory_types = [
    'early_to_late',                # Early stage → Late stage
    'stage1_to_stage4',             # Stage I → Stage IV
    'matched_pairs'                 # Patient-specific progressions
]
```

**Expected outputs**:
```
data/interim/20251217_synthetic_data/kirc/
├── early_to_late/
│   ├── TCGA-XXX-YYY_to_TCGA-AAA-BBB.csv
│   ├── TCGA-XXX-ZZZ_to_TCGA-CCC-DDD.csv
│   └── ...                         # 500 trajectory files
├── trajectory_metadata.csv         # Start/end points, stages
└── generation_params.json          # Generation parameters
```

**Validation**:
```python
import pandas as pd
import glob

# Count trajectories
trajectory_files = glob.glob('data/interim/*/kirc/early_to_late/*.csv')
print(f"Generated {len(trajectory_files)} trajectories")
# Expected: 500 trajectories

# Check trajectory format
traj = pd.read_csv(trajectory_files[0], index_col=0)
print(f"Trajectory shape: {traj.shape}")
# Expected: (20 timepoints, ~5000 genes)

assert traj.shape[0] == 20, "Wrong number of timepoints"
assert traj.isnull().sum().sum() == 0, "Missing values in trajectory"
print("✓ Trajectory validation passed")
```

---

### Step 5: Classification

**Script**: `scripts/pipeline_steps/5_classification.py`

**What it does**:
- Trains XGBoost classifier for stage prediction
- Calculates SHAP values for interpretability
- Identifies important gene signatures
- Evaluates performance with cross-validation

**Expected outputs**:
```
models/20251217_classification_kirc/
├── classifier.pkl                  # Trained XGBoost model
├── classification_metrics.json     # Accuracy, ROC AUC, etc.
├── shap_values.npy                 # SHAP importance values
├── important_genes.csv             # Top biomarker genes
└── figures/
    ├── confusion_matrix.png
    ├── roc_curve.png
    ├── shap_summary.png
    └── shap_waterfall.png
```

**Validation**:
```python
import json

# Check performance metrics
with open('models/20251217_classification_kirc/classification_metrics.json') as f:
    metrics = json.load(f)

print(f"Test Accuracy: {metrics['test_accuracy']:.3f}")
print(f"ROC AUC: {metrics['test_roc_auc']:.3f}")

# Expected performance for KIRC dataset
assert metrics['test_accuracy'] > 0.85, "Low accuracy!"
assert metrics['test_roc_auc'] > 0.90, "Low ROC AUC!"
print("✓ Classification performance acceptable")
```

---

### Step 6: Enrichment Analysis

**For the full tutorial, see**: [Enrichment Tutorial](step6-enrichment.md)

**Script**: `scripts/enrichment/pipeline.sh`

**What it does**:
- Calculates differential expression (DESeq-style) for each timepoint
- Runs GSEA on ranked gene lists
- Identifies enriched pathways along trajectories
- Combines results into single dataset

**Prerequisites**:
- GSEA CLI tool installed ([installation guide](../advanced/GSEA_INSTALLATION.md))
- ReactomePathways.gmt in `data/external/`


**Expected outputs**:
```
output_dir/
├── test_to_test/                                          # Synthetic trajectories
│   ├── early_to_late/                                     # Transition type
│   │   ├── patient1_to_patient2/                          # Patient trajectory
│   │   │   ├── patient1_to_patient_0.rnk                  # Ranked gene list for timepoint 0
│   │   │   ├── patient1_to_patient_1.rnk
│   │   │   ├── ...
│   │   │   ├── reports/                                   # GSEA output for all patients in directory
│   │   │   │   ├── patient1_to_patient_0.GseaPreranked.*  # GSEA output files
│   │   │   │   ├── gsea_report_for_na_pos_*.tsv           # Positive enrichment report
│   │   │   │   └── gsea_report_for_na_neg_*.tsv           # Negative enrichment report
│   │   ├──patient3_to_patient4/
│   │   ├── ...
│   └── gsea_commands_*.cmd                                # GSEA command files
├── full_gsea_reports_kirc.csv                             # Final combined results
└── heatmap_kirc_significantNES.csv                        # Significant pathways heatmap data
```

**Final result format**:
```csv
Patient,Idx,Transition,NAME,ES,NES,FDR q-val
TCGA-3Z-A93Z-01_to_TCGA-A3-A8OW-01,0,early_to_late,Cell Cycle,0.65,2.13,0.001
TCGA-3Z-A93Z-01_to_TCGA-A3-A8OW-01,0,early_to_late,DNA Repair,0.52,1.87,0.012
...
```

The Idx column represents the timepoint along the trajectory.

**Time**: approximate times for KIRC, after generating 50 timepoints of 279 trajectories:
- DESeq analysis: 8 hours, running in 10 nodes with 112 CPUs each (no multithreading), and assigning 21 CPUs per task, with 2GB per CPU.
- GSEA: 30 minutes, running in 5 nodes with 64 CPUs each (no multithreading), and assigning 1 CPUs per task, with 2GB per CPU.


---


## Citation

If you use this pipeline in your research, please cite this work.
For details, see the [citation page](../citation.md) file.

