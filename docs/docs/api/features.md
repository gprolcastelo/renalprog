# Features API

The `features` module provides functions for feature engineering and quality control of gene expression data.

## Overview

This module includes:

- Low expression gene filtering
- Mahalanobis distance outlier detection
- Gene clustering with tsfresh
- Feature extraction for trajectory analysis


## Core Functions

### filter_low_expression

Filter genes with low or invariant expression across samples.

::: renalprog.features.filter_low_expression

**Example Usage:**

```python
import pandas as pd
from renalprog.features import filter_low_expression

# Load raw expression data (genes × samples)
rnaseq = pd.read_csv("data/raw/KIRC_rnaseq.tsv", sep="\t", index_col=0)
print(f"Original: {rnaseq.shape[0]} genes")

# Filter low expression genes
filtered = filter_low_expression(
    rnaseq,
    mean_threshold=0.5,      # Minimum mean expression
    var_threshold=0.5,       # Minimum variance
    min_sample_fraction=0.2  # Maximum fraction of zero values
)
print(f"Filtered: {filtered.shape[0]} genes")
```

**Filtering Criteria:**

1. **Zero expression threshold**: Remove genes with >20% samples at zero
2. **Mean expression threshold**: Keep genes with mean ≥ 0.5
3. **Variance threshold**: Keep genes with variance ≥ 0.5

### detect_outliers_mahalanobis

Detect and remove outlier samples using robust Mahalanobis distance.

::: renalprog.features.detect_outliers_mahalanobis

**Example Usage:**

```python
import pandas as pd
from renalprog.features import detect_outliers_mahalanobis

# Load gene expression data (genes × samples)
rnaseq = pd.read_csv("data/interim/filtered_expression.csv", index_col=0)

# Detect outliers
cleaned_data, outlier_ids, mahal_distances = detect_outliers_mahalanobis(
    rnaseq,
    alpha=0.05,           # Significance level
    support_fraction=None,  # Auto-determine robust subset
    transpose=True,       # Transpose to samples × genes
    seed=42
)

print(f"Original samples: {rnaseq.shape[1]}")
print(f"Outliers detected: {len(outlier_ids)}")
print(f"Clean samples: {cleaned_data.shape[1]}")
print(f"Outlier IDs: {outlier_ids}")
```

**How It Works:**

1. **Minimum Covariance Determinant (MCD)**: Computes robust covariance estimate
2. **Mahalanobis Distance**: Calculates distance of each sample from the robust center
3. **Chi-square Test**: Identifies outliers exceeding chi-square threshold at significance level α

**Visualization:**

```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

# Plot Mahalanobis distances
n_features = rnaseq.shape[0]
cutoff = chi2.ppf(1 - 0.05, n_features)

plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.hist(mahal_distances, bins=50, edgecolor='black')
plt.axvline(cutoff, color='red', linestyle='--', label=f'Cutoff (α=0.05)')
plt.xlabel('Mahalanobis Distance')
plt.ylabel('Count')
plt.legend()
plt.title('Distribution of Mahalanobis Distances')

plt.subplot(1, 2, 2)
plt.scatter(range(len(mahal_distances)), sorted(mahal_distances))
plt.axhline(cutoff, color='red', linestyle='--')
plt.xlabel('Sample Index (sorted)')
plt.ylabel('Mahalanobis Distance')
plt.title('Sorted Mahalanobis Distances')
plt.tight_layout()
plt.savefig('outlier_detection.png', dpi=300)
```

## Quality Control Workflow

Complete QC pipeline for gene expression data:

```python
from pathlib import Path
import pandas as pd
from renalprog.features import filter_low_expression, detect_outliers_mahalanobis
from renalprog.config import PreprocessingConfig
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def qc_pipeline(rnaseq_path: Path, output_dir: Path):
    """Run complete quality control pipeline."""
    
    # Load configuration
    config = PreprocessingConfig()
    
    # Load data
    logger.info("Loading RNA-seq data...")
    rnaseq = pd.read_csv(rnaseq_path, sep="\t", index_col=0)
    logger.info(f"Initial shape: {rnaseq.shape}")
    
    # Filter low expression genes
    logger.info("Filtering low expression genes...")
    rnaseq_filtered = filter_low_expression(
        rnaseq,
        mean_threshold=config.mean_threshold,
        var_threshold=config.var_threshold,
        min_sample_fraction=config.min_sample_fraction
    )
    logger.info(f"After filtering: {rnaseq_filtered.shape}")
    
    # Detect and remove outliers
    logger.info("Detecting outliers...")
    rnaseq_clean, outliers, distances = detect_outliers_mahalanobis(
        rnaseq_filtered,
        alpha=config.outlier_alpha,
        seed=config.random_state
    )
    logger.info(f"Outliers removed: {len(outliers)}")
    logger.info(f"Final shape: {rnaseq_clean.shape}")
    
    # Save results
    output_dir.mkdir(parents=True, exist_ok=True)
    rnaseq_clean.to_csv(output_dir / "expression_qc.csv")
    
    # Save QC report
    qc_report = {
        'initial_genes': rnaseq.shape[0],
        'initial_samples': rnaseq.shape[1],
        'filtered_genes': rnaseq_filtered.shape[0],
        'genes_removed': rnaseq.shape[0] - rnaseq_filtered.shape[0],
        'outlier_samples': len(outliers),
        'outlier_ids': outliers,
        'final_genes': rnaseq_clean.shape[0],
        'final_samples': rnaseq_clean.shape[1]
    }
    
    pd.DataFrame([qc_report]).to_csv(output_dir / "qc_report.csv", index=False)
    
    return rnaseq_clean, qc_report

# Run pipeline
if __name__ == "__main__":
    rnaseq_clean, report = qc_pipeline(
        rnaseq_path=Path("data/raw/KIRC_rnaseq.tsv"),
        output_dir=Path("data/interim/qc_results")
    )
```

## Advanced Usage

### Custom Filtering Criteria

Define custom filtering thresholds for different datasets:

```python
from renalprog.features import filter_low_expression

# Strict filtering for high-quality datasets
strict_filtered = filter_low_expression(
    rnaseq,
    mean_threshold=1.0,
    var_threshold=1.0,
    min_sample_fraction=0.1
)

# Lenient filtering for smaller datasets
lenient_filtered = filter_low_expression(
    rnaseq,
    mean_threshold=0.1,
    var_threshold=0.1,
    min_sample_fraction=0.3
)
```

### Outlier Detection with Custom Parameters

Adjust sensitivity of outlier detection:

```python
from renalprog.features import detect_outliers_mahalanobis

# Conservative (fewer outliers)
conservative, outliers_con, _ = detect_outliers_mahalanobis(
    rnaseq, alpha=0.01, support_fraction=0.8
)

# Liberal (more outliers)
liberal, outliers_lib, _ = detect_outliers_mahalanobis(
    rnaseq, alpha=0.10, support_fraction=0.6
)

print(f"Conservative: {len(outliers_con)} outliers")
print(f"Liberal: {len(outliers_lib)} outliers")
```

## See Also

- [Dataset API](dataset.md) - Data loading and preparation
- [Configuration API](config.md) - Preprocessing configuration
- [Data Requirements Tutorial](../tutorials/data-requirements.md) - Data preparation guide

