# Dataset API

The `dataset` module provides functions for downloading, loading, and preparing RNA-seq data for analysis.

## Overview

This module handles:

- Downloading preprocessed data from Zenodo (fastest way to get started!)
- Downloading raw TCGA Pan-Cancer Atlas data
- Processing and filtering by cancer type
- Creating train/test splits
- Loading preprocessed data for modeling

!!! tip "Quick Start"
    Use `download_preprocessed_from_zenodo()` to quickly download ready-to-use datasets from Zenodo. This bypasses all preprocessing steps!
    Otherwise, use `download_data()` and `process_downloaded_data()` to process raw TCGA data with custom parameters.

## Core Functions

### download_preprocessed_from_zenodo

Download preprocessed RNAseq and clinical data directly from Zenodo repositories.

::: renalprog.dataset.download_preprocessed_from_zenodo

**Example Usage:**

```python
from renalprog.dataset import download_preprocessed_from_zenodo
from pathlib import Path

# Download preprocessed KIRC data
rnaseq, clinical = download_preprocessed_from_zenodo(
    rnaseq_url='https://zenodo.org/records/17987300/files/Static_KIRC.csv?download=1',
    clinical_url='https://zenodo.org/records/17987300/files/nodes_metadata.csv?download=1',
    output_dir=Path('data/interim/preprocessed_KIRC_data')
)

print(f"RNAseq shape: {rnaseq.shape}")
print(f"Clinical shape: {clinical.shape}")

# Download preprocessed BRCA data
rnaseq_brca, clinical_brca = download_preprocessed_from_zenodo(
    rnaseq_url='https://zenodo.org/records/17986123/files/Static_BRCA.csv?download=1',
    clinical_url='https://zenodo.org/records/17986123/files/nodes_metadata.csv?download=1',
    output_dir=Path('data/interim/preprocessed_BRCA_data')
)
```

**Available Zenodo Datasets:**

| Cancer Type | Zenodo Record | RNAseq File | Clinical File |
|-------------|---------------|-------------|---------------|
| **KIRC** (Kidney) | [17987300](https://zenodo.org/records/17987300) | `Static_KIRC.csv` | `nodes_metadata.csv` |
| **BRCA** (Breast) | [17986123](https://zenodo.org/records/17986123) | `Static_BRCA.csv` | `nodes_metadata.csv` |

!!! tip "Quick Start"
    This is the fastest way to get started! The data is already preprocessed, normalized, and ready for modeling.

!!! info "What's Included"
    - **RNAseq Data**: Log2-transformed, normalized gene expression matrix
    - **Clinical Data**: Patient metadata including stage, survival, demographics
    - **Automatic Validation**: Function checks sample overlap and data integrity

---

### download_data

Download raw TCGA datasets from UCSC Xena (TCGA data).

::: renalprog.dataset.download_data

**Example Usage:**

```python
from renalprog.dataset import download_data
from pathlib import Path

# Download to default location
rnaseq_path, clinical_path, phenotype_path = download_data(
    destination=Path("data/raw"),
    remove_gz=True,
    timeout=300
)

print(f"RNA-seq data: {rnaseq_path}")
print(f"Clinical data: {clinical_path}")
print(f"Phenotype data: {phenotype_path}")
```

### process_downloaded_data

Process downloaded TCGA data for a specific cancer type.

::: renalprog.dataset.process_downloaded_data

**Example Usage:**

```python
from renalprog.dataset import process_downloaded_data
from pathlib import Path

# Process data for KIRC (Kidney Renal Clear Cell Carcinoma)
rnaseq, clinical, phenotype = process_downloaded_data(
    rnaseq_path=Path("data/raw/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"),
    clinical_path=Path("data/raw/Survival_SupplementalTable_S1_20171025_xena_sp"),
    phenotype_path=Path("data/raw/TCGA_phenotype_denseDataOnlyDownload.tsv"),
    cancer_type="KIRC",
    output_dir=Path("data/raw"),
    early_late=False
)

print(f"RNA-seq shape: {rnaseq.shape}")
print(f"Clinical shape: {clinical.shape}")
```

### load_rnaseq_data

Load RNA-seq expression data from a file.

::: renalprog.dataset.load_rnaseq_data

### load_clinical_data

Load clinical and survival data from a file.

::: renalprog.dataset.load_clinical_data

### create_train_test_split

Create stratified train/test splits of the data.

::: renalprog.dataset.create_train_test_split

**Example Usage:**

```python
from renalprog.dataset import load_rnaseq_data, load_clinical_data, create_train_test_split
from pathlib import Path

# Load the data
rnaseq = load_rnaseq_data(Path("data/raw/KIRC_rnaseq.tsv"))
clinical = load_clinical_data(Path("data/raw/KIRC_clinical.tsv"))

# Create train/test split
create_train_test_split(
    rnaseq=rnaseq,
    clinical=clinical,
    test_size=0.2,
    random_state=42,
    output_dir=Path("data/interim/my_experiment")
)
```

### load_train_test_split

Load previously saved train/test split data.

::: renalprog.dataset.load_train_test_split

**Example Usage:**

```python
from renalprog.dataset import load_train_test_split
from pathlib import Path

# Load existing split
train_expr, test_expr, train_clin, test_clin = load_train_test_split(
    Path("data/interim/my_experiment")
)

print(f"Training samples: {len(train_expr)}")
print(f"Test samples: {len(test_expr)}")
```

### map_stages_to_early_late

Map cancer stages to binary early/late categories.

::: renalprog.dataset.map_stages_to_early_late

## Data Format Requirements

### RNA-seq Expression Data

Expected format:
- Rows: Genes (with gene symbols or Ensembl IDs)
- Columns: Samples (patient IDs)
- Values: Log2-transformed TPM or FPKM expression values

```python
# Example structure
#              TCGA-A1-A0SB  TCGA-A1-A0SD  TCGA-A1-A0SE  ...
# GENE_A       5.234         4.891         6.123         ...
# GENE_B       2.456         2.789         2.634         ...
# GENE_C       8.912         9.234         8.756         ...
```

### Clinical Data

Expected columns:
- `sample`: Patient ID matching expression columns
- `OS`: Overall survival status (0=alive, 1=deceased)
- `OS.time`: Overall survival time (days)

Optional columns:
- `age_at_initial_pathologic_diagnosis`: Age at diagnosis
- `gender`: Patient gender
- `tumor_stage`: Tumor stage

```python
# Example structure
#    sample          OS  OS.time  age  gender  stage
# 0  TCGA-A1-A0SB    1   1825     65   MALE    IV
# 1  TCGA-A1-A0SD    0   2190     58   FEMALE  II
```

## Train/Test Splitting

The module provides stratified splitting to maintain class balance:

```python
from renalprog.dataset import (
    load_rnaseq_data, 
    load_clinical_data, 
    create_train_test_split
)
from pathlib import Path

# Load the data
rnaseq = load_rnaseq_data(Path("data/raw/KIRC_rnaseq.tsv"))
clinical = load_clinical_data(Path("data/raw/KIRC_clinical.tsv"))

# Create stratified split preserving early/late survival distribution
create_train_test_split(
    rnaseq=rnaseq,
    clinical=clinical,
    test_size=0.2,  # 20% test set
    random_state=42,  # For reproducibility
    output_dir=Path("data/interim/my_split")
)

# Load the split data
train_expr, test_expr, train_clin, test_clin = load_train_test_split(
    Path("data/interim/my_split")
)

# Check class distribution
import pandas as pd
train_dist = train_clin.value_counts(normalize=True)
test_dist = test_clin.value_counts(normalize=True)

print("Training set distribution:")
print(train_dist)
print("\nTest set distribution:")
print(test_dist)
```

## Data Preprocessing Pipeline

### Option 1: Quick Start with Zenodo

Download preprocessed data directly from Zenodo - no preprocessing needed!

```python
from pathlib import Path
from renalprog.dataset import download_preprocessed_from_zenodo

# Download and you're ready to go!
rnaseq, clinical = download_preprocessed_from_zenodo(
    rnaseq_url='https://zenodo.org/records/17987300/files/Static_KIRC.csv?download=1',
    clinical_url='https://zenodo.org/records/17987300/files/nodes_metadata.csv?download=1',
    output_dir=Path('data/interim/preprocessed_KIRC_data')
)

# Data is already preprocessed and ready for modeling
print(f"Ready to use: {rnaseq.shape[0]:,} genes × {rnaseq.shape[1]:,} samples")
```

---

### Option 2: Custom Preprocessing Pipeline (Advanced)

The standard preprocessing pipeline for raw TCGA data:

1. **Download raw data** from TCGA
2. **Filter by cancer type** (e.g., KIRC)
3. **Filter low expression genes** (see [Features API](features.md))
4. **Remove outliers** using Mahalanobis distance
5. **Create train/test split** with stratification
6. **Normalize** using MinMaxScaler (0-1 range)
7. **Save preprocessed data** for modeling

**Complete Example:**

```python
from pathlib import Path
from renalprog.dataset import (
    download_data, 
    process_downloaded_data,
    load_rnaseq_data,
    load_clinical_data,
    create_train_test_split,
    load_train_test_split
)
from renalprog.features import filter_low_expression, detect_outliers_mahalanobis

# Step 1: Download
rnaseq_path, clinical_path, phenotype_path = download_data(
    destination=Path("data/raw")
)

# Step 2: Process for KIRC
rnaseq, clinical, _ = process_downloaded_data(
    rnaseq_path=rnaseq_path,
    clinical_path=clinical_path,
    phenotype_path=phenotype_path,
    cancer_type="KIRC",
    output_dir=Path("data/raw")
)

# Step 3: Filter low expression
rnaseq_filtered = filter_low_expression(
    rnaseq,
    mean_threshold=0.5,
    var_threshold=0.5
)

# Step 4: Remove outliers
rnaseq_clean, outliers, _ = detect_outliers_mahalanobis(
    rnaseq_filtered,
    alpha=0.05
)

# Step 5: Create train/test split
rnaseq_clean_path = Path("data/interim/rnaseq_clean.csv")
clinical_path = Path("data/raw/KIRC_clinical.tsv")
rnaseq_clean.to_csv(rnaseq_clean_path)

# Load and split
rnaseq_final = load_rnaseq_data(rnaseq_clean_path)
clinical_final = load_clinical_data(clinical_path)

create_train_test_split(
    rnaseq=rnaseq_final,
    clinical=clinical_final,
    test_size=0.2,
    random_state=42,
    output_dir=Path("data/interim/20251218_experiment")
)
```

## See Also

- [Features API](features.md) - Gene filtering and outlier detection
- [Preprocessing Tutorial](../tutorials/data-requirements.md) - Step-by-step data preparation
- [Configuration API](config.md) - Data paths and preprocessing parameters

