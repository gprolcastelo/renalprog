# Step 3: Reconstruction Quality Check

This guide explains how to evaluate the reconstruction quality of trained VAE and Reconstruction Network models.

## Overview

The reconstruction quality check pipeline performs the following steps:

1. **Model Loading**: Load trained VAE and Reconstruction Network models
2. **Data Reconstruction**: Apply models to preprocessed data
3. **UMAP Visualization**: Create UMAP plots comparing original and reconstructed data
4. **Quality Metrics** (optional): Calculate reconstruction quality metrics using sdmetrics

The pipeline supports two model loading options:
- **Local models**: Load models trained in Step 2 from local disk
- **Hugging Face models**: Download and use pre-trained models from Hugging Face Hub

## Prerequisites

Before running the reconstruction check pipeline, ensure you have:

- **Trained models**: Either from Step 2 or available on Hugging Face
- **Preprocessed data**: From Step 1 data processing
- **Python environment**: With PyTorch, pandas, plotly installed
- **Optional**: sdmetrics package for quality metrics (only if using `--sdmetrics`)

## Usage

### Basic Usage (Local Models)

```bash
python scripts/pipeline_steps/3_check_reconstruction.py
```

This will:
- Load locally trained KIRC models
- Load preprocessed KIRC data
- Generate reconstructions with VAE and Reconstruction Network
- Create UMAP visualizations
- Save all outputs

### Using Hugging Face Pre-trained Models

```bash
python scripts/pipeline_steps/3_check_reconstruction.py --hf_models
```

This will:
- Download pre-trained KIRC models from Hugging Face
- Load preprocessed KIRC data
- Generate reconstructions
- Create UMAP visualizations

### For BRCA Cancer Type

```bash
# With local models
python scripts/pipeline_steps/3_check_reconstruction.py --cancer_type BRCA

# With Hugging Face models
python scripts/pipeline_steps/3_check_reconstruction.py --cancer_type BRCA --hf_models
```

### With Quality Metrics Evaluation

```bash
python scripts/pipeline_steps/3_check_reconstruction.py --sdmetrics
```

**Warning**: Using `--sdmetrics` can take a very long time (several hours) as it computes comprehensive quality metrics.

## Command-Line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--cancer_type` | String | 'KIRC' | Cancer type to process (KIRC or BRCA) |
| `--hf_models` | Flag | False | Load pre-trained models from Hugging Face Hub |
| `--sdmetrics` | Flag | False | Compute comprehensive quality metrics (very slow!) |

## Model Loading

### Local Models

When not using `--hf_models`, the script expects models in:

```
models/models_{CANCER_TYPE}/
├── vae/
│   ├── final_model.pth           # VAE weights
│   └── config.json               # VAE architecture
├── reconstruction_network.pth    # Reconstruction Network weights
└── network_dims.csv              # Network dimensions
```

### Hugging Face Models

When using `--hf_models`, the script downloads from:

**Repository**: `gpcastelo/evenflow_models`

**Files downloaded**:
- `{CANCER_TYPE}/config.json` - VAE configuration
- `{CANCER_TYPE}/20250321_VAE_idim8516_md512_feat256mse_relu.pth` (KIRC) - VAE weights
- `{CANCER_TYPE}/20251209_VAE_idim8954_md1024_feat512mse_relu.pth` (BRCA) - VAE weights
- `{CANCER_TYPE}/network_dims.csv` - Reconstruction Network dimensions
- `{CANCER_TYPE}/network_reconstruction.pth` - Reconstruction Network weights

**Note**: Hugging Face models are cached locally after first download.

## Processing Steps

### Step 1: Load Preprocessed Data

Loads data from Step 1:
- `data/interim/preprocessed_{CANCER_TYPE}_data/preprocessed_rnaseq.csv`
- `data/interim/preprocessed_{CANCER_TYPE}_data/clinical_data.csv`

### Step 2: Load Models

Either from local disk or Hugging Face, including:
- VAE model with configuration
- Reconstruction Network with layer dimensions

### Step 3: Generate Reconstructions

1. **VAE Reconstruction**:
   - Encode data to latent space
   - Decode back to gene expression space
   - Save latent representations

2. **Reconstruction Network Refinement**:
   - Take VAE output as input
   - Apply postprocessing network
   - Generate refined reconstructions

### Step 4: Create UMAP Visualizations

Generate multiple UMAP plots:

1. **Original Data**: UMAP of preprocessed data
2. **VAE Reconstruction**: UMAP of VAE output
3. **Reconstruction Network Output**: UMAP of postprocessed data
4. **Original vs VAE**: Comparison plot
5. **Original vs RecNet**: Comparison plot

### Step 5: Quality Metrics (Optional)

If `--sdmetrics` is enabled:
- Calculate Bland-Altman plots
- Compute Kolmogorov-Smirnov statistics
- Generate comprehensive quality reports

## Output Files

### UMAP Visualizations

```
reports/figures/YYYYMMDD_{CANCER_TYPE}_umap_reconstruction/
├── preprocessed.html                    # Original data UMAP
├── VAE_output.html                      # VAE reconstruction UMAP
├── reconstruction_network_output.html   # RecNet output UMAP
├── preprocessed_and_vae.html           # Comparison: original vs VAE
└── preprocessed_and_recnet.html        # Comparison: original vs RecNet
```

### Quality Metrics (if using --sdmetrics)

```
data/processed/YYYYMMDD_reconstruction_evaluation/
├── vae_reconstruction/
│   ├── bland_altman_plots/
│   ├── ks_statistics.csv
│   └── bland_altman_summary.csv
└── recnet_reconstruction/
    ├── bland_altman_plots/
    ├── ks_statistics.csv
    └── bland_altman_summary.csv
```

## Complete Example Workflow

### Quick Check with Local Models

```bash
# Step 1: Train models (from Step 2)
python scripts/pipeline_steps/2_models.py

# Step 2: Check reconstruction quality
python scripts/pipeline_steps/3_check_reconstruction.py
```

### Reproducible Check with Hugging Face Models

```bash
# Download and use pre-trained models
python scripts/pipeline_steps/3_check_reconstruction.py --hf_models --cancer_type KIRC
```

### Comprehensive Evaluation

```bash
# Full evaluation with quality metrics (takes several hours!)
python scripts/pipeline_steps/3_check_reconstruction.py --hf_models --sdmetrics
```

# Troubleshooting

### Hugging Face Download Fails

If model download from Hugging Face fails:
- Check your internet connection
- Verify Hugging Face Hub is accessible
- Try authenticating: `huggingface-cli login`
- Check repository is public: `https://huggingface.co/gpcastelo/evenflow_models`

### CUDA Out of Memory

If GPU memory is insufficient:
- The script uses `force_cpu=True` by default
- If you enabled GPU, reduce batch size or use CPU

### Dimension Mismatch

If you see dimension mismatch errors:
- Verify models match the cancer type
- Check preprocessed data has correct number of genes
- Ensure models were trained on same preprocessed data

### sdmetrics Takes Too Long

The `--sdmetrics` option can take several hours:
- Run in background: `nohup python ... --sdmetrics &`
- Use a smaller subset for testing
- Consider skipping if not needed for your analysis

## Next Steps

After verifying reconstruction quality:

1. **If reconstruction is good**: Proceed to Step 4 (trajectory generation)
2. **If reconstruction needs improvement**: Return to Step 2 and adjust training parameters

```bash
# Proceed to trajectory generation
python scripts/pipeline_steps/4_trajectories.py
```

## Additional Resources

- [UMAP Documentation](https://umap-learn.readthedocs.io/)
- [Bland-Altman Analysis](https://en.wikipedia.org/wiki/Bland%E2%80%93Altman_plot)
- [Kolmogorov-Smirnov Test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test)
- [SDMetrics Documentation](https://docs.sdv.dev/sdmetrics/)
- [Hugging Face Hub Documentation](https://huggingface.co/docs/hub/index)
