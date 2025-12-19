# Tutorials Overview

This section provides comprehensive tutorials for using `renalprog` to analyze kidney cancer progression. Each tutorial is designed to be self-contained yet builds upon previous steps.

## Tutorial Structure

### 🚀 Quick Start
New to `renalprog`? Start here!

- **[Quick Start](quickstart.md)**: 10-minute introduction to core functionality
- **[Using Pretrained Models](pretrained-models.md)**: Fastest way to reproduce paper results
- **[Data Requirements](data-requirements.md)**: Understanding input data formats

### 📋 Complete Pipeline
Step-by-step walkthrough of the entire analysis pipeline:

1. **[Data Processing](step1-data-processing.md)**: Download, filter, and preprocess TCGA data
2. **[VAE Training](step2-vae-training.md)**: Train variational autoencoders
3. **[Reconstruction Validation](step3-reconstruction.md)**: Assess model quality
4. **[Trajectory Generation](step4-trajectories.md)**: Create synthetic progression paths
5. **[Classification](step5-classification.md)**: Stage prediction and biomarker discovery
6. **[Enrichment Analysis](step6-enrichment.md)**: Pathway analysis with GSEA

### 🎨 Visualization
- **[Visualization Guide](visualization.md)**: Create publication-quality figures

## Learning Paths

### For Biologists
If you're primarily interested in **biological insights**:

1. Start with [Quick Start](quickstart.md)
2. Read [Data Requirements](data-requirements.md) to understand the data
3. Follow [Complete Pipeline](complete-pipeline.md) end-to-end
4. Focus on [Enrichment Analysis](step6-enrichment.md) for pathway interpretation
5. Use [Visualization Guide](visualization.md) for publication figures

### For Computational Scientists
If you want to **customize models** or **develop new methods**:

1. Complete [Quick Start](quickstart.md)
2. Study [VAE Training](step2-vae-training.md) in detail
3. Explore [API Reference](../api/index.md) for implementation details
4. Read [Custom Models](../advanced/custom-models.md) for extending functionality
5. Review [Architecture](../advanced/architecture.md) for design principles

### For Reproducing Published Results
If you want to **reproduce the paper**:

**Option 1: Using Pretrained Models (Recommended)**

1. Follow [Using Pretrained Models](pretrained-models.md) tutorial
2. This is the fastest and most accurate way to reproduce results
3. Uses the exact models from the paper

**Option 2: Training from Scratch**

1. Follow [System Requirements](../reproducibility/requirements.md)
2. Complete [Data Preparation](../reproducibility/data-preparation.md)
3. Execute [Running the Pipeline](../reproducibility/pipeline.md)
4. Validate using [Expected Results](../reproducibility/results.md)

## Tutorial Conventions

### Code Blocks

Python code to execute:
```python
from renalprog import dataset
data = dataset.load_data('path/to/data.csv')
```

Shell commands:
```bash
python scripts/pipeline_steps/1_data_processing.py
```

### Callouts

!!! note "Note"
    Informational notes provide additional context.

!!! tip "Tip"
    Tips offer helpful suggestions and best practices.

!!! warning "Warning"
    Warnings highlight potential issues or common pitfalls.

!!! danger "Danger"
    Critical warnings about data loss or major errors.

!!! example "Example"
    Example outputs or usage patterns.

### File Paths

All file paths are relative to the repository root unless otherwise specified:

```
renalprog/
├── data/
│   ├── raw/           # Downloaded TCGA data
│   ├── interim/       # Intermediate processing outputs
│   └── processed/     # Final processed data
├── models/            # Trained models
├── reports/           # Analysis results
└── scripts/           # Pipeline scripts
```

## Prerequisites

Before starting these tutorials, ensure you have:

1. ✅ Installed `renalprog` ([Installation Guide](../getting-started.md))
2. ✅ Python 3.9+ and R 4.0+ available
3. ✅ At least 16 GB RAM (8 GB minimum)
4. ✅ 50+ GB free disk space
5. ✅ (Optional) CUDA-capable GPU for faster training

## Time Estimates

| Tutorial | Reading Time | Execution Time | Difficulty |
|----------|--------------|----------------|------------|
| Quick Start | 10 min | 15 min | ⭐ Easy |
| Data Processing | 15 min | 30 min | ⭐ Easy |
| VAE Training | 20 min | 2-4 hours* | ⭐⭐ Moderate |
| Reconstruction | 10 min | 10 min | ⭐ Easy |
| Trajectories | 15 min | 30 min | ⭐⭐ Moderate |
| Classification | 20 min | 1 hour | ⭐⭐ Moderate |
| Enrichment | 25 min | 2-6 hours* | ⭐⭐⭐ Advanced |
| Visualization | 15 min | 30 min | ⭐⭐ Moderate |

\* With GPU acceleration and 8+ CPU cores. Times may be significantly longer on older hardware.

## Getting Help

If you encounter issues while following these tutorials:

1. Check the [Troubleshooting Guide](../reproducibility/troubleshooting.md)
2. Review the [API Reference](../api/index.md) for function details
3. Search [GitHub Issues](https://github.com/gprolcastelo/renalprog/issues)
4. Ask in [GitHub Discussions](https://github.com/gprolcastelo/renalprog/discussions)

## Next Steps

Ready to begin? Start with the [Quick Start Tutorial](quickstart.md)!

