<div align="center">
  <img src="docs/docs/assets/images/kidneys.png" alt="renalprog logo" width="200"/>
</div>

# renalprog

**A Python package for simulating kidney cancer progression with synthetic data generation and machine learning.**

[![Cookiecutter Data Science](https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter)](https://cookiecutter-data-science.drivendata.org)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

<p align="center">
  <sub>Logo: <a href="https://www.flaticon.com/free-icons/kidneys" title="kidneys icons">Kidneys icons created by Smashicons - Flaticon</a></sub>
</p>

## Overview

`renalprog` is a comprehensive bioinformatics pipeline for analyzing kidney cancer (KIRC) progression using deep learning and pathway enrichment analysis. The package integrates Variational Autoencoders (VAEs) with differential expression analysis and gene set enrichment to model and interpret cancer progression trajectories.

### Key Features

- **Data Preprocessing**: Automated filtering of low-expression genes and robust outlier detection using Mahalanobis distance
- **Deep Learning Models**: Variational Autoencoder (VAE) and Conditional VAE (CVAE) implementations for learning latent representations
- **Trajectory Generation**: Generate synthetic patient trajectories between cancer stages
- **Stage Classification**: XGBoost-based classification of early vs. late stage cancer
- **Enrichment Analysis**: Integration with R-based DESeq2 and GSEA for pathway analysis
- **Visualization**: Comprehensive plotting functions for all analysis steps


### For New Users
- Follow the [Quick Start Tutorial](docs/docs/tutorials/quickstart.md)
- Explore [Step-by-Step Tutorials](docs/docs/tutorials/index.md)

### For Paper Reproducers
- See [Reproducibility Guide](docs/docs/reproducibility/index.md) for exact replication of published results

### For Contributors
- Check [Contributing Guide](CONTRIBUTING.md)
- See [Documentation TODO](DOCUMENTATION_TODO.md) for ways to help

## Installation

> 📖 **For detailed installation instructions**, see [INSTALLATION.md](INSTALLATION.md)

### Prerequisites

- Python 3.9 or higher
- R 4.0+ (for enrichment analysis) - **Can be installed via conda/mamba** (recommended)
- Conda or Mamba (recommended for environment management)
- CUDA-capable GPU (optional, for faster VAE training)

### Recommended: Using Mamba/Conda + uv

#### Quick Setup with environment.yml (Easiest)

```bash
# Clone repository
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog

# Create environment from file (includes Python, R, and all dependencies)
mamba env create -f environment.yml
mamba activate renalprog

# Install the package in editable mode
pip install -e .
```

#### Manual Setup

```bash
# Clone repository
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog

# Create environment with Python 3.9 AND R
mamba create -n renalprog "python==3.9" "r-base>=4.0"
mamba activate renalprog

# Install R packages via conda (recommended for reproducibility)
mamba install -c conda-forge r-gprofiler2 r-ggplot2 r-optparse

# Install uv for faster Python package management
pip install uv

# Install Python package
uv pip install -e .

# Install testing dependencies
uv pip install pytest pytest-cov
```

**Note**: Installing R via conda/mamba ensures all dependencies are managed in the same environment, improving reproducibility.

### Alternative: Using pip + venv

```bash
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -e .
```

### With Development Dependencies

```bash
# Using uv
uv pip install -e ".[dev]"

# Or using pip
pip install -e ".[dev]"
```

### Quick Setup Script

For automated setup (Linux/Mac):
```bash
chmod +x quick_setup.sh
./quick_setup.sh
```

### R Dependencies for Enrichment Analysis

The package includes R scripts for gene enrichment analysis. You can install R dependencies in several ways:

#### Option 1: Via Conda/Mamba (Recommended)

If using conda/mamba environment:

```bash
# Activate your environment
mamba activate renalprog

# Install R and packages (if not done during initial setup)
mamba install -c conda-forge r-base r-gprofiler2 r-ggplot2 r-optparse
```

#### Option 2: Via R's install.packages()

If you already have R installed system-wide or prefer CRAN:

```bash
# Install R packages using the provided script
Rscript scripts/r_analysis/install_r_packages.R
```

On Windows PowerShell:
```powershell
Rscript scripts\r_analysis\install_r_packages.R
```

**Required R packages:**
- `r-gprofiler2` / `gprofiler2` - Gene enrichment via g:Profiler API
- `r-ggplot2` / `ggplot2` - Visualization
- `r-optparse` / `optparse` - Command-line parsing

See [`scripts/r_analysis/README.md`](scripts/r_analysis/README.md) for detailed R setup and usage.

## Quick Start

### Using Pretrained Models (Recommended for Paper Reproduction)

The fastest way to reproduce paper results is using the provided pretrained models:

```bash
# 1. Generate trajectories using pretrained KIRC model
python scripts/pipeline_steps/use_pretrained_model.py \
    --cancer_type KIRC \
    --model_dir models/pretrained/KIRC \
    --data_dir data/interim/preprocessed_KIRC \
    --output_dir data/processed/trajectories_KIRC_pretrained

# 2. Run enrichment analysis
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --trajectory_dir data/processed/trajectories_KIRC_pretrained/early_to_late/test_to_test \
    --output_dir data/processed/enrichment_KIRC_pretrained \
    --n_threads 8 \
    --gsea_path ./GSEA_4.3.2/gsea-cli.sh \
    --pathways_file data/external/ReactomePathways.gmt

# 3. Generate pathway heatmaps
python scripts/pipeline_steps/6b_generate_pathway_heatmap.py \
    --enrichment_file data/processed/enrichment_KIRC_pretrained/trajectory_enrichment.csv \
    --output_dir data/processed/enrichment_KIRC_pretrained \
    --fdr_threshold 0.05
```

**Pretrained models available for**:
- KIRC (Kidney Renal Clear Cell Carcinoma)
- BRCA (Breast Invasive Carcinoma)

See [Pretrained Models Tutorial](docs/docs/tutorials/pretrained-models.md) for detailed documentation.

### Training from Scratch

#### 1. Data Download

Download TCGA-KIRC data from the [UCSC Xena Browser](https://xenabrowser.net/):

1. Navigate to TCGA Kidney Renal Clear Cell Carcinoma (KIRC)
2. Download:
   - Gene expression (HiSeqV2)
   - Clinical matrix
   - Phenotype data

#### 1. Data Download

Download TCGA-KIRC data from the [UCSC Xena Browser](https://xenabrowser.net/):

1. Navigate to TCGA Kidney Renal Clear Cell Carcinoma (KIRC)
2. Download:
   - Gene expression (HiSeqV2)
   - Clinical matrix
   - Phenotype data

Place downloaded files in `data/raw/KIRC/`.

#### 2. Preprocessing

```python
from renalprog import features, dataset
from renalprog.config import get_dated_dir, INTERIM_DATA_DIR

# Preprocess RNA-seq data
import pandas as pd

data = pd.read_csv("data/raw/KIRC/expression.csv", index_col=0)
processed_data, metadata = features.preprocess_rnaseq(
    data,
    filter_expression=True,
    detect_outliers=True
)

# Save preprocessed data
output_dir = get_dated_dir(INTERIM_DATA_DIR, "preprocessed_KIRC")
features.save_preprocessing_results(processed_data, metadata, output_dir)
```

#### 3. Create Train/Test Split

```python
from renalprog.dataset import create_train_test_split
from renalprog.config import get_dated_dir, INTERIM_DATA_DIR

output_dir = get_dated_dir(INTERIM_DATA_DIR, "train_test_split")

X_train, X_test, y_train, y_test, _, _ = create_train_test_split(
    rnaseq_path="data/interim/preprocessed_KIRC/rnaseq_preprocessed.csv",
    clinical_path="data/interim/preprocessed_KIRC/clinical.csv",
    test_size=0.2,
    seed=2023,
    output_dir=output_dir
)
```

#### 4. Train VAE

```python
from renalprog.modeling.train import train_vae
from renalprog.config import VAEConfig, get_dated_dir, MODELS_DIR

# Configure VAE
config = VAEConfig()
config.INPUT_DIM = X_train.shape[1]
config.MID_DIM = 1024
config.LATENT_DIM = 16
config.EPOCHS = 100

# Train model
model_dir = get_dated_dir(MODELS_DIR, "VAE_KIRC")
model, history = train_vae(
    X_train, X_test,
    y_train, y_test,
    config=config,
    save_dir=model_dir
)
```

#### 5. Generate Trajectories

```python
from renalprog.modeling.predict import generate_trajectories

trajectories = generate_trajectories(
    model=model,
    source_samples=early_stage_samples,
    target_samples=late_stage_samples,
    n_steps=50,
    method="linear"
)
```

## Pipeline Overview

The complete pipeline consists of 13 steps:

1. **Data Download** (manual): Download TCGA-KIRC data from Xena Browser
2. **Preprocessing**: Filter low-expression genes and detect outliers
3. **Train/Test Split**: Create stratified splits maintaining stage distribution
4. **Train VAE**: Learn latent representation of gene expression
5. **Adjust Reconstruction**: Train postprocessing network for better reconstruction
6. **Assess Reconstruction**: Evaluate quality using SDMetrics
7. **Create Connections**: Identify patient pairs for trajectory generation
8. **Generate Trajectories**: Create synthetic progression paths
9. **Generate Controls**: Create noise-based control trajectories
10. **Classification**: Train XGBoost classifier for stage prediction
11. **Apply Classifier**: Classify synthetic trajectories
12. **Dynamic Enrichment**: GSEA analysis on trajectory timepoints
13. **Static Enrichment**: Pathway enrichment on gene clusters (R script)

### Running the Full Pipeline

The complete pipeline can be run step-by-step:

```bash
# Steps 1-5: Python pipeline (data processing, VAE training, trajectory generation)
python scripts/pipeline_steps/1_data_processing.py
python scripts/pipeline_steps/2_models.py
python scripts/pipeline_steps/3_check_reconstruction.py
python scripts/pipeline_steps/4_trajectories.py
python scripts/pipeline_steps/5_classification.py

# Step 6: Dynamic enrichment analysis (requires GSEA)
python scripts/pipeline_steps/6_enrichment_analysis.py \
  --trajectory_dir data/interim/20251216_synthetic_data/kirc/early_to_late \
  --output_dir data/processed/20251217_enrichment \
  --n_threads 8 \
  --gsea_path ./GSEA_4.3.2/gsea-cli.sh

# Step 7: Static enrichment analysis (R)
Rscript scripts/r_analysis/gene_enrichment.R \
  --input data/external/important_genes_shap.csv \
  --sources "GO,KEGG,REAC,WP" \
  --top_n 20
```

See the [full pipeline documentation](docs/pipeline.md) for details.

### GSEA Installation (Required for Step 6)

For dynamic enrichment analysis, you need to install GSEA:

1. Download from: https://www.gsea-msigdb.org/gsea/downloads.jsp
2. Extract to project root (creates `GSEA_4.3.2/` directory)
3. See [GSEA Installation Guide](docs/docs/GSEA_INSTALLATION.md) for detailed instructions

## Project Structure

```
renalprog/
├── LICENSE                      # Apache 2.0 license
├── README.md                    # This file
├── pyproject.toml              # Package configuration
├── setup.cfg                    # Tool configurations
├── requirements.txt             # Python dependencies
├── Makefile                     # Convenience commands
├── data/
│   ├── external/               # Gene lists, pathway databases
│   ├── interim/                # Intermediate processed data
│   ├── processed/              # Final analysis outputs
│   └── raw/                    # Original TCGA data (not in git)
├── models/                     # Trained model checkpoints
├── notebooks/                  # Jupyter notebooks for exploration
├── reports/                    # Generated analysis reports
│   └── figures/               # Publication figures
├── scripts/                    # Standalone scripts
│   ├── r_analysis/            # R scripts (DESeq2, clusterProfiler)
│   └── pipeline_steps/        # Individual pipeline step scripts
├── tests/                      # Unit and integration tests
└── renalprog/                  # Main package
    ├── __init__.py
    ├── config.py               # Configuration and paths
    ├── dataset.py              # Data loading and splitting
    ├── features.py             # Preprocessing and feature engineering
    ├── plots.py                # Visualization functions
    ├── modeling/
    │   ├── __init__.py
    │   ├── train.py           # Model training
    │   └── predict.py         # Inference and trajectory generation
    └── utils/
        └── __init__.py        # Utility functions
```

## Testing

### Quick Start

Run the test suite:

```bash
pytest tests/
```

Run tests with coverage:

```bash
pytest tests/ --cov=renalprog --cov-report=html
```

The test suite includes a lightweight synthetic dataset that runs in <5 minutes on CPU.

### Test Categories

Tests are organized with markers for selective execution:

- `slow`: Tests that take longer to run (e.g., model training)
- `requires_gpu`: Tests requiring GPU acceleration
- `requires_heavy_deps`: Tests requiring heavy dependencies (SDV, full PyTorch, etc.)
- `integration`: Integration tests with real data/models

**Run fast tests only** (recommended for quick validation):
```bash
pytest tests/ -m "not slow"
```

**Run CI-compatible tests** (this is what runs in GitHub Actions):
```bash
pytest tests/ -m "not slow and not requires_heavy_deps"
```

## Documentation

📚 **Comprehensive documentation is available in the `docs/` directory.**

### Quick Links

- **[Installation Guide](INSTALLATION.md)** - Detailed installation instructions
- **[Quick Start Tutorial](docs/docs/tutorials/quickstart.md)** - Get started in 5 minutes
- **[Complete Pipeline Guide](docs/docs/reproducibility/pipeline.md)** - Reproduce paper results
- **[API Reference](docs/docs/api/index.md)** - Function and class documentation
- **[GSEA Setup Guide](docs/GSEA_SETUP_GUIDE.md)** - GSEA installation and configuration
- **[Citations](docs/CITATIONS.md)** - Proper citations for all tools and databases

### Documentation Contents

- **Getting Started**
  - Installation and setup
  - Data requirements
  - Quick start tutorial

- **Step-by-Step Tutorials**
  - [Step 1: Data Processing](docs/docs/tutorials/step1-data-processing.md)
  - [Step 2: VAE Training](docs/docs/tutorials/step2-vae-training.md)
  - [Step 3: Reconstruction Validation](docs/docs/tutorials/step3-reconstruction.md)
  - [Step 4: Trajectory Generation](docs/docs/tutorials/step4-trajectories.md)
  - [Step 5: Stage Classification](docs/docs/tutorials/step5-classification.md)
  - [Step 6: Enrichment Analysis](docs/docs/tutorials/step6-enrichment.md)
  - [Visualization Guide](docs/docs/tutorials/visualization.md)

- **Reproducibility**
  - System requirements
  - Data preparation
  - Running the pipeline
  - Expected results
  - Troubleshooting

- **API Reference**
  - Configuration (`renalprog.config`)
  - Dataset (`renalprog.dataset`)
  - Features (`renalprog.features`)
  - VAE Models (`renalprog.modeling.vae`)
  - Training (`renalprog.modeling.train`)
  - Prediction (`renalprog.modeling.predict`)
  - Trajectories (`renalprog.modeling.trajectories`)
  - Classification (`renalprog.modeling.classification`)
  - Enrichment (`renalprog.modeling.enrichment`)
  - Visualization (`renalprog.plots`)
  - Utilities (`renalprog.utils`)

- **Advanced Topics**
  - Architecture details
  - GSEA integration
  - R integration
  - Custom models
  - Performance tuning
  - GPU acceleration

- **Contributing**
  - Development guidelines
  - Code style
  - Testing

### Build Documentation Locally

The documentation uses MkDocs with Material theme:

```bash
# Install documentation dependencies
pip install -r docs/requirements.txt

# Serve documentation locally with live reload
cd docs
mkdocs serve

# Open http://127.0.0.1:8000/ in your browser
```

See [docs/BUILD.md](docs/BUILD.md) for complete documentation building instructions.

### Deploy Documentation

Deploy to GitHub Pages:

```bash
cd docs
mkdocs gh-deploy
```

**Note:** Documentation will be available at https://yourusername.github.io/renalprog/ after deployment.

## Troubleshooting

### Common Issues

#### GSEA Enrichment Analysis Errors

If you encounter errors during enrichment analysis, particularly "No GSEA results found to combine":

**Quick Diagnosis:**
```bash
python scripts/debug_gsea_results.py --deseq_dir data/processed/YYYYMMDD_enrichment/deseq
```

**Documentation:**
- 🚀 **[Quick Reference](docs/GSEA_QUICK_REFERENCE.md)** - Fast fixes for common issues
- 📖 **[Troubleshooting Guide](docs/GSEA_TROUBLESHOOTING.md)** - Comprehensive diagnostic steps
- 🔧 **[Implementation Details](docs/GSEA_ERROR_HANDLING_IMPROVEMENTS.md)** - Technical documentation

**Common causes:**
1. GSEA not installed or wrong path - see [Installation Guide](INSTALLATION.md#gsea-setup)
2. Java version incompatible (need 11+) - `java -version`
3. GSEA commands failed - check `failed_gsea_commands.txt`
4. Retry failed commands: `python scripts/retry_failed_gsea.py --failed_commands <file>`

#### Other Issues

For other issues, check:
- Installation problems: See [INSTALLATION.md](INSTALLATION.md)
- Pipeline errors: Check logs in `logs/` directory
- R integration issues: See [R Setup Guide](docs/R_SETUP_GUIDE.md)
- Performance issues: See documentation section on "Performance Tuning"

**Getting Help:**
1. Check existing [GitHub Issues](https://github.com/gprolcastelo/renalprog/issues)
2. Review documentation at https://gprolcastelo.github.io/renalprog/
3. Open a new issue with:
   - Error message and full traceback
   - Python and package versions
   - Steps to reproduce
   - Relevant log files

## Examples

See the `notebooks/` directory for example analyses:

- `01_preprocessing_example.ipynb`: Data preprocessing walkthrough
- `02_vae_training.ipynb`: VAE training and evaluation
- `03_trajectory_generation.ipynb`: Creating and analyzing trajectories
- `04_enrichment_analysis.ipynb`: Pathway enrichment analysis
- `paper_figures.ipynb`: Reproduce paper figures

## Citation

> ⚠️ **TEMPORARY CITATION**: This citation will be updated once the preprint/publication is available. Please check back regularly or watch the repository for updates.

If you use RenalProg in your research, please cite it as follows until the official publication is available:

```bibtex
@software{renalprog2024,
  title = {RenalProg: A Deep Learning Framework for Kidney Cancer Progression Modeling},
  author = {[Guillermo Prol-Castelo, Elina Syrri, Nikolaos Manginas, Vasileos Manginas, Nikos Katzouris, Davide Cirillo, George Paliouras, Alfonso Valencia]},
  year = {2025},
  url = {https://github.com/gprolcas/renalprog},
  note = {Preprint in preparation}
}
```

**Status**: Preprint in preparation | **Expected publication**: 2026

For the most up-to-date citation information, please visit the [Citation page](https://gprolcastelo.github.io/renalprog/citation/) in our documentation.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

1. Clone the repository
2. Install development dependencies: `pip install -e ".[dev]"`
3. Install pre-commit hooks: `pre-commit install`
4. Run tests to verify setup: `pytest tests/`

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- TCGA Research Network for KIRC data
- Barcelona Supercomputing Center (BSC-CNS)
- [EVENFLOW EU Project](https://evenflow-project.eu)
- All contributors and collaborators

## Contact

For questions and support, please open an issue on [GitHub](https://github.com/gprolcastelo/renalprog/issues).

---

**Note**: This package is under active development. APIs may change between versions.
