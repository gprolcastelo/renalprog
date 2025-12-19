<div align="center">
  <img src="assets/images/kidneys.png" alt="renalprog logo" width="200"/>
</div>

# renalprog

**A Python package for simulating kidney cancer progression with synthetic data generation and machine learning**

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

<p align="center">
  <sub>Logo: <a href="https://www.flaticon.com/free-icons/kidneys" title="kidneys icons">Kidneys icons created by Smashicons - Flaticon</a></sub>
</p>


!!! warning
    This site is under construction. Some pages may be missing.

---

## Overview

`renalprog` is a comprehensive bioinformatics pipeline for analyzing kidney cancer (KIRC) progression using deep learning and pathway enrichment analysis. The package integrates **Variational Autoencoders (VAEs)** with **differential expression analysis** and **gene set enrichment** to model and predict cancer progression trajectories.

### Scientific Context

Cancer progression is a complex, dynamic process involving multiple molecular alterations across time. Traditional static analyses fail to capture the temporal dynamics of tumor evolution. `renalprog` addresses this challenge by:

1. **Learning latent representations** of gene expression data using deep generative models (VAEs)
2. **Generating synthetic trajectories** between cancer stages in latent space
3. **Identifying enriched biological pathways** along progression trajectories
4. **Classifying cancer stages** using interpretable machine learning

This approach enables researchers to:

- Identify key biological pathways driving cancer progression
- Predict patient outcomes based on molecular profiles
- Generate testable hypotheses about therapeutic targets
- Understand the temporal dynamics of tumor evolution

---

## Key Features

### 🔬 Data Processing
- Automated filtering of low-expression genes
- Robust outlier detection using Mahalanobis distance
- Normalization and batch effect correction
- Integration with TCGA and other genomics datasets

### 🧠 Deep Learning Models
- **Variational Autoencoder (VAE)** for unsupervised representation learning
- **Conditional VAE (CVAE)** for stage-specific modeling
- Support for custom architectures and hyperparameters
- GPU acceleration for large-scale datasets

### 🔄 Trajectory Generation
- Generate synthetic patient trajectories between cancer stages
- Interpolation in latent space with biological constraints
- Multiple trajectory types (early-to-late, stage-specific, custom)
- Quality control and validation metrics

### 📊 Stage Classification
- XGBoost-based classification of early vs. late stage cancer
- SHAP values for feature importance and interpretability
- Cross-validation and performance evaluation
- Gene signature discovery

### 🧬 Enrichment Analysis
- Integration with **DESeq2** for differential expression
- **GSEA** (Gene Set Enrichment Analysis) for pathway analysis
- Support for Reactome, KEGG, and custom pathway databases
- Parallel processing for large-scale analyses

### 📈 Visualization
- Comprehensive plotting functions for all analysis steps
- UMAP/t-SNE visualizations of latent space
- Pathway enrichment heatmaps
- Classification performance metrics
- Interactive plots with Plotly

---

## Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog

# Create environment (includes Python and R)
mamba env create -f environment.yml
mamba activate renalprog

# Install package
pip install -e .
```

See [Installation Guide](getting-started.md) for detailed instructions.

### Basic Usage

```python
from renalprog import config, dataset, modeling

# Load and process data
data = dataset.load_data('data/processed/rnaseq_maha.csv')
processed = dataset.preprocess(data)

# Train VAE
vae = modeling.VAE(input_dim=processed.shape[1])
vae.train(processed, epochs=100)

# Generate trajectories
trajectories = modeling.generate_trajectories(
    vae=vae,
    start_stage='early',
    end_stage='late',
    n_samples=100
)

# Run enrichment analysis
results = modeling.enrichment_analysis(trajectories)
```

See [Quick Start Tutorial](tutorials/quickstart.md) for a complete example.

---

## Pipeline Overview

The `renalprog` pipeline consists of six main steps:

```mermaid
graph LR
    A[Raw Data] --> B[1. Data Processing]
    B --> C[2. VAE Training]
    C --> D[3. Reconstruction Check]
    D --> E[4. Trajectory Generation]
    E --> F[5. Classification]
    F --> G[6. Enrichment Analysis]
    G --> H[Results & Visualization]
```

### Step 1: Data Processing
- Filter low-expression genes
- Remove outliers using Mahalanobis distance
- Normalize expression values
- Prepare clinical metadata

### Step 2: VAE Training
- Train deep generative models on gene expression data
- Learn low-dimensional latent representations
- Validate reconstruction quality

### Step 3: Reconstruction Validation
- Assess VAE reconstruction accuracy
- Visualize latent space structure
- Identify potential issues

### Step 4: Trajectory Generation
- Generate synthetic patient trajectories
- Interpolate between cancer stages
- Export trajectory gene expression

### Step 5: Classification
- Train XGBoost classifier for stage prediction
- Calculate SHAP values for interpretability
- Identify important gene signatures

### Step 6: Enrichment Analysis
- Differential expression analysis with DESeq2
- Pathway enrichment with GSEA
- Identify biological processes along trajectories

---

## Documentation Structure

### 📚 For New Users
Start with:

1. [Installation Guide](getting-started.md) - Set up your environment
2. [Quick Start Tutorial](tutorials/quickstart.md) - Run your first analysis
3. [Complete Pipeline Tutorial](tutorials/complete-pipeline.md) - End-to-end workflow

### 🔬 For Reproducibility
Reproduce published results:

1. [System Requirements](reproducibility/requirements.md) - Hardware and software needs
2. [Data Preparation](reproducibility/data-preparation.md) - Download and prepare data
3. [Running the Pipeline](reproducibility/pipeline.md) - Step-by-step execution
4. [Expected Results](reproducibility/results.md) - Validate your outputs

### 🛠️ For Developers
Extend and customize:

1. [API Reference](api/index.md) - Complete function documentation
2. [Architecture Guide](advanced/architecture.md) - Design principles
3. [Custom Models](advanced/custom-models.md) - Implement new architectures
4. [Contributing Guidelines](contributing/guidelines.md) - Join development

---

## System Requirements

### Minimum Requirements
- **OS**: Linux, macOS, or Windows (with WSL for enrichment analysis)
- **Python**: 3.9 or higher
- **R**: 4.0 or higher (for enrichment analysis)
- **RAM**: 8 GB
- **Storage**: 20 GB free space

### Recommended Requirements
- **RAM**: 16+ GB
- **CPU**: 8+ cores
- **GPU**: CUDA-capable GPU with 6+ GB VRAM (for VAE training)
- **Storage**: 50+ GB on SSD

### Software Dependencies
- PyTorch 2.0+
- scikit-learn 1.0+
- XGBoost 1.5+
- pandas, numpy, scipy
- R packages: DESeq2, gprofiler2

See [System Requirements](reproducibility/requirements.md) for complete details.

---

## Use Cases

### Cancer Research
- Model tumor evolution over time
- Identify driver pathways in progression
- Predict patient outcomes
- Discover therapeutic targets

### Computational Biology
- Learn representations of high-dimensional genomics data
- Generate synthetic data for validation
- Integrate multi-omics datasets
- Perform pathway-level analysis

### Machine Learning Research
- Apply VAEs to biological data
- Develop interpretable deep learning models
- Benchmark generative models
- Study latent space interpolation

---

## Citation

If you use `renalprog` in your research, please cite:

```bibtex
@software{renalprog2025,
  author = {Prol-Castelo, Guillermo and EVENFLOW Project},
  title = {renalprog: Simulating Kidney Cancer Progression with Generative AI},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/gprolcastelo/renalprog}
}
```

See [How to Cite](citation.md) for additional references.

---

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](license.md) file for details.

---

## Support

- **Issues**: [GitHub Issues](https://github.com/gprolcastelo/renalprog/issues)
- **Discussions**: [GitHub Discussions](https://github.com/gprolcastelo/renalprog/discussions)
- **Email**: Contact the EVENFLOW Project team

---

## Acknowledgments

This work is supported by the EVENFLOW Project and builds upon numerous open-source tools and databases. See [Acknowledgments](acknowledgments.md) for complete credits.

---

## Quick Links

- [Installation Guide](getting-started.md)
- [Quick Start Tutorial](tutorials/quickstart.md)
- [API Reference](api/index.md)
- [Reproducibility Guide](reproducibility/index.md)
- [Contributing Guidelines](contributing/guidelines.md)

