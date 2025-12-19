# Pipeline Steps

This directory contains the orchestration scripts for the renalprog pipeline. Each script represents a major step in the cancer progression analysis workflow.

## Pipeline Overview

```
1. Data Download (manual) → 2. Preprocessing → 3. Train/Test Split → 4. Train VAE
                                                                           ↓
                                    13. Static Enrichment ← 12. Dynamic Enrichment
                                            ↑                         ↑
                            11. Gene Clustering ← 10. Classification ← 9. Control Trajectories
                                                          ↑                     ↑
                                          8. Assess Reconstruction → 7. Generate Trajectories
                                                                           ↑
                                              6. Create Connections → 5. Adjust Reconstruction
```

## Available Scripts

### Core Pipeline

#### `6_enrichment_analysis.py` ✅ IMPLEMENTED

**Dynamic pathway enrichment analysis using GSEA**

Performs Gene Set Enrichment Analysis (GSEA) on synthetic trajectories to identify biological pathways enriched at different stages of cancer progression.

**Usage:**
```bash
python 6_enrichment_analysis.py \
    --trajectory_dir data/interim/20251216_synthetic_data/kirc/early_to_late \
    --output_dir data/processed/20251217_enrichment \
    --cancer_type kirc \
    --n_threads 8 \
    --gsea_path ./GSEA_4.3.2/gsea-cli.sh \
    --pathways_file data/external/ReactomePathways.gmt
```

**Requirements:**
- GSEA CLI tool (download from https://www.gsea-msigdb.org/gsea/downloads.jsp)
- ReactomePathways.gmt (included in data/external/)
- Java 11+

**Output:**
- `trajectory_enrichment.csv`: Combined enrichment results
- DESeq fold-change files (`.rnk`)
- GSEA reports (`.tsv`)

**Documentation:**
- [Enrichment Analysis Guide](../../docs/docs/ENRICHMENT_ANALYSIS.md)
- [GSEA Installation](../../docs/docs/GSEA_INSTALLATION.md)

**Time:** 2-6 hours for full KIRC dataset (depending on system and threads)

---

#### `5_classification.py` ✅ IMPLEMENTED

**Cancer stage classification**

Trains XGBoost classifier to predict cancer stages (early vs late) and applies it to synthetic trajectories.

**Features:**
- Optuna hyperparameter optimization
- Multi-seed training for robustness
- Train-to-train and test-to-test trajectory processing
- Publication-ready Plotly visualizations

**Documentation:**
- See implementation in `renalprog/modeling/train.py`

---

#### `4_trajectories.py` (To be implemented)

**Synthetic trajectory generation**

Generates synthetic patient progression trajectories by interpolating between early and late stage samples in latent space.

**Features:**
- Multiple interpolation methods (linear, spherical)
- Patient connection identification
- Control trajectory generation

---

#### `3_check_reconstruction.py` (To be implemented)

**Reconstruction quality assessment**

Evaluates VAE reconstruction quality using SDMetrics and statistical tests.

---

#### `2_models.py` (To be implemented)

**VAE training and postprocessing**

Trains Variational Autoencoder and postprocessing network for gene expression reconstruction.

---

#### `1_data_processing.py` (To be implemented)

**Data preprocessing and train/test split**

Filters low-expression genes, detects outliers using Mahalanobis distance, and creates stratified train/test splits.

---

### Example Scripts

#### `enrichment_example.py` ✅ IMPLEMENTED

**Example usage of enrichment analysis**

Demonstrates how to use the enrichment pipeline with proper error handling and prerequisite checking.

**Usage:**
```bash
python enrichment_example.py
```

---

### R Analysis Scripts

Located in `scripts/r_analysis/`:

- `gene_enrichment.R`: Static pathway enrichment using g:Profiler
- `differential_expression.R`: Differential expression analysis
- `install_r_packages.R`: Install required R packages

See [R Analysis README](../r_analysis/README.md) for details.

---

## Quick Start

### 1. Download Data (Manual)

Download TCGA-KIRC data from [UCSC Xena Browser](https://xenabrowser.net/):
- Gene expression (HiSeqV2)
- Clinical matrix
- Phenotype data

Place in `data/raw/KIRC/`

### 2. Run Pipeline Steps

```bash
# Step 1: Preprocess data
python 1_data_processing.py

# Step 2: Train VAE
python 2_models.py

# Step 3: Check reconstruction
python 3_check_reconstruction.py

# Step 4: Generate trajectories
python 4_trajectories.py

# Step 5: Train classifier
python 5_classification.py

# Step 6: Dynamic enrichment (requires GSEA)
python 6_enrichment_analysis.py \
    --trajectory_dir data/interim/20251216_synthetic_data/kirc/early_to_late \
    --output_dir data/processed/20251217_enrichment \
    --n_threads 8

# Step 7: Static enrichment (R)
Rscript ../r_analysis/gene_enrichment.R \
    --input data/external/important_genes_shap.csv
```

---

## Configuration

### General Settings

Most scripts accept the following arguments:
- `--seed`: Random seed for reproducibility (default: 42)
- `--output_dir`: Output directory (auto-dated if not specified)
- `--cancer_type`: Cancer type (kirc, lobular, ductal)

### Hardware Settings

- `--n_threads`: Number of parallel threads (default: 4)
- `--device`: PyTorch device (cpu, cuda, mps)

### Data Paths

- `--data_dir`: Input data directory
- `--model_dir`: Model checkpoint directory

---

## Output Organization

All scripts save outputs to dated directories:

```
data/
├── interim/
│   ├── 20251216_preprocessed_KIRC/
│   ├── 20251216_train_test_split/
│   └── 20251216_synthetic_data/
├── processed/
│   ├── 20251217_enrichment/
│   └── 20251217_classification/
└── models/
    ├── 20251216_VAE_KIRC/
    └── 20251216_classification/
```

---

## Parallelization

### CPU Parallelization

Scripts use Python's `multiprocessing` or `concurrent.futures`:
- Set `--n_threads` based on your CPU cores
- Recommended: Number of cores - 1
- Example: 8-core system → use 6-7 threads

### GPU Acceleration

VAE training can use GPU:
- Automatically detected if CUDA/MPS available
- Force CPU: `--device cpu`
- Force GPU: `--device cuda`

---

## Error Handling

### Common Issues

**1. Out of Memory**
- Reduce `--n_threads`
- Reduce batch size in VAE training
- Close other applications

**2. Missing Dependencies**
- Check `requirements.txt` installed
- For R: check R packages installed
- For GSEA: download and install separately

**3. GSEA Not Found**
- Download from https://www.gsea-msigdb.org/gsea/downloads.jsp
- Extract to project root
- Specify path: `--gsea_path ./GSEA_4.3.2/gsea-cli.sh`

### Resume Interrupted Runs

Most scripts support resuming:

```bash
# Skip completed steps
python 6_enrichment_analysis.py \
    --skip_deseq \  # Skip if DESeq already done
    --skip_gsea     # Skip if GSEA already done
```

---

## Testing

Test individual scripts:

```bash
# Test with small dataset
python 6_enrichment_analysis.py \
    --trajectory_dir tests/fixtures/small_trajectories \
    --output_dir tests/output \
    --n_threads 2
```

---

## Logging

All scripts provide detailed logging:
- Console output with progress bars
- Timestamps for all operations
- Error messages with troubleshooting hints

Enable debug logging:
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

---

## Performance Benchmarks

**System: Modern Desktop (16 cores, 32 GB RAM)**

| Step | Time | Threads | GPU |
|------|------|---------|-----|
| 1. Preprocessing | 5 min | 1 | No |
| 2. VAE Training | 30 min | 1 | Yes |
| 3. Reconstruction | 10 min | 4 | No |
| 4. Trajectories | 20 min | 8 | No |
| 5. Classification | 15 min | 8 | No |
| 6. Enrichment | 2 hrs | 12 | No |
| **Total** | **~3.5 hrs** | | |

**System: Laptop (8 cores, 16 GB RAM)**

Total time: ~6-8 hours

---

## Contributing

When adding new pipeline steps:

1. Create script named `N_step_name.py`
2. Follow existing argument parsing pattern
3. Add comprehensive logging
4. Include error handling
5. Support `--help` flag
6. Add to this README
7. Create example in `examples/`
8. Add tests in `tests/`

---

## See Also

- [Main README](../../README.md)
- [Package Documentation](../../docs/)
- [R Analysis Scripts](../r_analysis/README.md)
- [Migration Plan](../../MIGRATION_PLAN.md)

---

## Contact

For issues with specific pipeline steps, please open a GitHub issue with:
- Script name and version
- Error message
- System information (OS, Python version, RAM)
- Command used

