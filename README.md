<div align="center">
  <img src="docs/docs/assets/images/kidneys.png" alt="renalprog logo" width="200"/>
</div>

# renalprog

**A Python package for simulating kidney cancer progression with synthetic data generation and machine learning.**

[![Cookiecutter Data Science](https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter)](https://cookiecutter-data-science.drivendata.org)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

<p align="center">
  <sub>Logo: <a href="https://www.flaticon.com/free-icons/kidneys" title="kidneys icons">Kidneys icons created by Smashicons - Flaticon</a></sub>
</p>

## Overview

`renalprog` is a comprehensive bioinformatics pipeline for analyzing kidney cancer (KIRC) progression using deep learning and pathway enrichment analysis. The package integrates Variational Autoencoders (VAEs) with differential expression analysis and gene set enrichment to model and interpret cancer progression trajectories.

### Key Features

- **Data Preprocessing**: Automated filtering of low-expression genes and robust outlier detection using Mahalanobis distance
- **Deep Learning Models**: Variational Autoencoder (VAE), Conditional VAE (CVAE), and Autoencoder (AE) implementations for learning latent representations
- **Trajectory Generation**: Generate synthetic patient trajectories between cancer stages
- **Stage Classification**: XGBoost-based classification of early vs. late stage cancer
- **Enrichment Analysis**: Integration with pyDESeq2 and GSEA for pathway analysis
- **Visualization**: Comprehensive plotting functions for all analysis steps

## 📚 Documentation

Comprehensive documentation is available with guides for users, contributors, and reproducers:

- **[📖 Full Documentation](https://gprolcastelo.github.io/renalprog/)** - Complete documentation website
- **[🚀 Quick Start](https://gprolcastelo.github.io/renalprog/tutorials/quickstart/)** 

### For Paper Reproducers
- See the [Tutorial](https://gprolcastelo.github.io/renalprog/tutorials/) for replication of published results.

### For Contributors
- Check [Contributing Guide](https://gprolcastelo.github.io/renalprog/contributing/guidelines/)

## Installation

> 📖 **For detailed installation instructions**, see [Installation Guide](https://gprolcastelo.github.io/renalprog/tutorials/installation/

### Prerequisites

- Python 3.9 or higher
- R 4.0+ (for enrichment analysis) - **Can be installed via conda/mamba** (recommended)
- Conda or Mamba (recommended for environment management)
- CUDA-capable GPU (optional, for faster VAE training)

### Recommended: Mamba

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

**Note**: Installing R via conda/mamba ensures all dependencies are managed in the same environment.

### Quick Setup Script

For automated setup (Linux/Mac):
```bash
chmod +x quick_setup.sh
./quick_setup.sh
```

## Quick Start

### Using Pretrained Models (Recommended for Paper Reproduction)

The fastest way to reproduce paper results is using the provided pretrained models, available on HuggingFace: [gprolcastelo/evenflow_models](https://huggingface.co/gprolcastelo/evenflow_models)



See [Pretrained Models Tutorial](https://gprolcastelo.github.io/renalprog/tutorials/pretrained-models/) for detailed documentation.

## Pipeline Overview

The complete pipeline consists of these steps:

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

```

See the [full pipeline documentation](docs/pipeline.md) for details.

## Project Structure

```
renalprog/
├── LICENSE                      # Apache 2.0 license
├── README.md                    # This file
├── pyproject.toml               # Package configuration
├── setup.cfg                    # Tool configurations
├── requirements.txt             # Python dependencies
├── Makefile                     # Convenience commands
├── data/
│   ├── external/                # Gene lists, pathway databases
│   ├── interim/                 # Intermediate processed data
│   ├── processed/               # Final analysis outputs
│   └── raw/                     # Original TCGA data (not in git)
├── models/                      # Trained model checkpoints
├── notebooks/                   # Jupyter notebooks for exploration
├── reports/                     # Generated analysis reports
│   └── figures/                 # Publication figures
├── scripts/                     # Standalone scripts
│   ├── r_analysis/              # R scripts (DESeq2, clusterProfiler)
│   └── pipeline_steps/          # Individual pipeline step scripts
├── tests/                       # Unit and integration tests
└── renalprog/                   # Main package
    ├── __init__.py
    ├── config.py                # Configuration and paths
    ├── dataset.py               # Data loading and splitting
    ├── enrichment.py            # Enrichment analysis functions
    ├── features.py              # Preprocessing and feature engineering
    ├── plots.py                 # Visualization functions
    ├── modeling/
    │   ├── __init__.py
    │   ├── train.py             # Model training
    │   └── predict.py           # Inference and trajectory generation
    └── utils/
        └── __init__.py          # Utility functions
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

## Troubleshooting

Getting Help:

1. Check existing [GitHub Issues](https://github.com/gprolcastelo/renalprog/issues)
2. Review documentation at https://gprolcastelo.github.io/renalprog/
3. Open a new issue with:
   - Error message and full traceback
   - Python and package versions
   - Steps to reproduce
   - Relevant log files

## Citation

> ⚠️ **TEMPORARY CITATION**: This citation will be updated once the preprint/publication is available. Please check back regularly or watch the repository for updates.

If you use `renalprog` in your research, please cite it as follows until the official publication is available:

```bibtex
@software{renalprog2024,
  title = {renalprog: A Deep Learning Framework for Kidney Cancer Progression Modeling},
  author = {[Guillermo Prol-Castelo, Elina Syrri, Nikolaos Manginas, Vasileos Manginas, Nikos Katzouris, Davide Cirillo, George Paliouras, Alfonso Valencia]},
  year = {2025},
  url = {https://github.com/gprolcas/renalprog},
  note = {Preprint in preparation}
}
```

**Status**: Preprint in preparation | **Expected publication**: 2026

For the most up-to-date citation information, please visit the [Citation page](https://gprolcastelo.github.io/renalprog/citation/) in our documentation.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](docs/docs/contributing/guidelines.md) for guidelines.

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- TCGA Research Network for KIRC data
- Barcelona Supercomputing Center (BSC-CNS)
- [EVENFLOW EU Project](https://evenflow-project.eu)
- All contributors and collaborators

## Contact

For questions and support, please open an issue on [GitHub](https://github.com/gprolcastelo/renalprog/issues).
