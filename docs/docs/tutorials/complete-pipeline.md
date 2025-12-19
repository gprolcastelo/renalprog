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

**Total Time**: ~8-12 hours (varies with hardware)

**Requirements**: 
- 16+ GB RAM
- 50+ GB free disk space
- GPU recommended for steps 2-4
- R 4.0+ for step 6

## Prerequisites

Ensure you have completed the installation:

```bash
mamba env create -f environment.yml
mamba activate renalprog
pip install -e .
```

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

# Step 6: Enrichment Analysis (2-6 hours with 8 cores)
python 6_enrichment_analysis.py
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
make enrichment
```

### Option 3: Use Job Scripts (HPC)

If you're on an HPC cluster with SLURM:

```bash
# Submit all jobs
sbatch jobs/job_data
sbatch jobs/job_models
sbatch jobs/job_trajectories
sbatch jobs/job_classification
sbatch jobs/job_enrichment
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

**Time**: ~30 minutes  
**See**: [Data Processing Tutorial](step1-data-processing.md)

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

**Time**: 2-4 hours (GPU), 8-12 hours (CPU)  
**See**: [VAE Training Tutorial](step2-vae-training.md)

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

**Time**: ~10 minutes  
**See**: [Reconstruction Tutorial](step3-reconstruction.md)

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

**Time**: ~30 minutes  
**See**: [Trajectory Tutorial](step4-trajectories.md)

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

**Time**: ~1 hour  
**See**: [Classification Tutorial](step5-classification.md)

---

### Step 6: Enrichment Analysis

**Script**: `scripts/pipeline_steps/6_enrichment_analysis.py`

**What it does**:
- Calculates differential expression (DESeq-style) for each timepoint
- Runs GSEA on ranked gene lists
- Identifies enriched pathways along trajectories
- Combines results into single dataset

**Prerequisites**:
- GSEA CLI tool installed ([installation guide](../GSEA_INSTALLATION.md))
- ReactomePathways.gmt in `data/external/`

**Configuration**:
```python
n_threads = 8                       # Parallel processing threads
gsea_path = './GSEA_4.3.2/gsea-cli.sh'
pathways_file = 'data/external/ReactomePathways.gmt'
```

**Expected outputs**:
```
data/processed/20251217_enrichment/
├── deseq/
│   └── early_to_late/
│       ├── TCGA-XXX_to_YYY/
│       │   ├── patient_tp0_foldchange.rnk
│       │   ├── patient_tp1_foldchange.rnk
│       │   └── gsea_tp0/
│       │       ├── gsea_report_for_na_pos_*.tsv
│       │       └── gsea_report_for_na_neg_*.tsv
│       └── ...
└── trajectory_enrichment.csv       # Final combined results
```

**Final result format**:
```csv
Patient,Idx,Transition,NAME,ES,NES,FDR q-val
TCGA-3Z-A93Z-01_to_TCGA-A3-A8OW-01,0,early_to_late,Cell Cycle,0.65,2.13,0.001
TCGA-3Z-A93Z-01_to_TCGA-A3-A8OW-01,0,early_to_late,DNA Repair,0.52,1.87,0.012
...
```

**Validation**:
```python
import pandas as pd

# Load enrichment results
enrichment = pd.read_csv('data/processed/20251217_enrichment/trajectory_enrichment.csv')

print(f"Total enrichment results: {len(enrichment)}")
print(f"Unique patients: {enrichment['Patient'].nunique()}")
print(f"Unique pathways: {enrichment['NAME'].nunique()}")

# Check for significant pathways
sig_pathways = enrichment[enrichment['FDR q-val'] < 0.05]
print(f"Significant pathways (FDR < 0.05): {len(sig_pathways)}")

assert len(enrichment) > 0, "No enrichment results!"
assert enrichment['Patient'].nunique() == 500, "Missing trajectories!"
print("✓ Enrichment analysis complete")
```

**Time**: 2-6 hours (depends on CPU cores)  
**See**: [Enrichment Tutorial](step6-enrichment.md)

---

## Checking Pipeline Outputs

After running the complete pipeline, verify all outputs:

```bash
# Run validation script
python scripts/validate_pipeline_outputs.py
```

Or manually check key files:

```python
import os
import pandas as pd

# Check each step's outputs
checks = {
    'Step 1': 'data/interim/preprocessed_KIRC_data/preprocessed_rnaseq.csv',
    'Step 2': 'models/20251217_models_KIRC/vae/vae_model.pt',
    'Step 3': 'reports/figures/reconstruction/latent_space_umap.png',
    'Step 4': 'data/interim/20251217_synthetic_data/kirc/trajectory_metadata.csv',
    'Step 5': 'models/20251217_classification_kirc/classifier.pkl',
    'Step 6': 'data/processed/20251217_enrichment/trajectory_enrichment.csv',
}

for step, filepath in checks.items():
    exists = os.path.exists(filepath)
    status = "✓" if exists else "✗"
    print(f"{status} {step}: {filepath}")
```

## Expected Disk Usage

After completing the pipeline:

```
data/raw/                    ~5 GB    (TCGA downloads)
data/interim/                ~15 GB   (trajectories + intermediate files)
data/processed/              ~10 GB   (enrichment results)
models/                      ~1 GB    (trained models)
reports/                     ~100 MB  (figures)
-------------------------------------------
Total:                       ~31 GB
```

## Performance Benchmarks

Typical runtimes on different systems:

| System | Step 2 | Step 6 | Total | Config |
|--------|--------|--------|-------|--------|
| HPC Cluster | 1 hr | 1 hr | 3 hrs | 48 cores, A100 GPU |
| Desktop (2023) | 2 hrs | 2 hrs | 6 hrs | 16 cores, RTX 4080 |
| Desktop (2020) | 3 hrs | 4 hrs | 9 hrs | 8 cores, GTX 1080 |
| Laptop (2018) | 8 hrs | 6 hrs | 16 hrs | 4 cores, CPU only |

## Troubleshooting

### Pipeline Fails at Step X

Resume from the failing step:

```bash
# Fix the issue, then run from that step onward
python scripts/pipeline_steps/X_step_name.py
```

### Out of Disk Space

Clean up intermediate files:

```bash
# Remove GSEA temporary files
rm -rf data/processed/*/deseq/*/gsea_*

# Remove old dated directories
rm -rf data/interim/202412*  # Keep only latest
```

### Out of Memory

Reduce batch sizes and threads:

```python
# In 2_models.py
vae_config.BATCH_SIZE = 4  # Reduce from 8

# In 6_enrichment_analysis.py
n_threads = 2  # Reduce from 8
```

## Next Steps

After completing the pipeline:

1. **Analyze Results**: Explore the enrichment results
2. **Create Visualizations**: Generate publication figures
3. **Interpret Pathways**: Map significant pathways to biology
4. **Validate Findings**: Test predictions experimentally

See:
- [Visualization Tutorial](visualization.md)
- [Result Interpretation Guide](../reproducibility/results.md)
- [API Reference](../api/index.md) for custom analyses

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{renalprog2024,
  author = {Prol-Castelo, Guillermo},
  title = {renalprog: Cancer Progression Forecasting with Generative AI},
  year = {2024},
  url = {https://github.com/gprolcastelo/renalprog}
}
```

