# Data Requirements

Understanding the input data formats and requirements for `renalprog`.

## Overview

`renalprog` is designed to work with **gene expression data** and **clinical annotations**. The primary use case is TCGA (The Cancer Genome Atlas) data, but the package can work with any properly formatted RNA-seq or microarray data.

## Input Data Formats

### 1. Gene Expression Data

**Format**: CSV or TSV file with genes as columns and samples as rows (or vice versa).

**Example** (`rnaseq.csv`):
```csv
,Gene1,Gene2,Gene3,...,GeneN
Sample1,12.5,8.3,0.1,...,5.7
Sample2,13.1,7.9,0.0,...,6.2
Sample3,11.8,9.2,0.3,...,5.1
...
```

**Requirements**:

- **Samples**: Minimum 50, recommended 200+
- **Genes**: Minimum 1000, recommended 5000+
- **Values**: 
  - Can be raw counts, TPM, FPKM, or log-transformed
  - Will be preprocessed by the pipeline
  - Must be numeric (no missing values in final data)
- **Index**: Sample IDs as row index OR column names
- **Columns**: Gene IDs or symbols

**Supported Formats**:

- ✅ Samples × Genes (preferred)
- ✅ Genes × Samples (will be transposed)
- ✅ Log-transformed or linear scale
- ✅ Normalized or unnormalized (pipeline handles normalization)

### 2. Clinical Data

**Format**: CSV or TSV file with samples as rows and clinical variables as columns.

**Example** (`clinical.csv`):
```csv
sample_id,ajcc_pathologic_tumor_stage,age,gender,survival_days
Sample1,Stage I,65,M,1825
Sample2,Stage IIIA,72,F,730
Sample3,Stage II,58,M,2190
...
```

**Required Columns**:

- `sample_id`: Must match gene expression sample IDs
- `ajcc_pathologic_tumor_stage`: Cancer stage (or similar staging variable)
  - Examples: "Stage I", "Stage II", "Stage IIIA", "Stage IV"
  - Or: "T1N0M0", "T2N1M0", etc.

**Optional Columns** (highly recommended):

- `age`: Patient age at diagnosis
- `gender` or `sex`: M/F or Male/Female
- `survival_days`: Days of follow-up
- `vital_status`: Alive/Dead or 0/1
- `tumor_grade`: G1, G2, G3, G4
- `histological_type`: Tumor histology

**Stage Format**:

The package can handle various stage formats:

```python
# AJCC staging
"Stage I", "Stage IA", "Stage IB"
"Stage II", "Stage IIA", "Stage IIB"
"Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"
"Stage IV", "Stage IVA", "Stage IVB"

# TNM staging
"T1N0M0", "T2N1M0", "T3N2M1"

# Simplified
"Early", "Late"
"Low", "High"
```

The pipeline automatically groups stages:

- **Early**: Stage I, II
- **Late**: Stage III, IV

### 3. Sample Metadata (Optional)

**Format**: CSV file with additional sample annotations.

**Example** (`phenotype.csv`):
```csv
sample_id,sample_type,tissue_source_site,batch,platform
Sample1,Primary Tumor,TCGA-XX,Batch1,IlluminaHiSeq
Sample2,Primary Tumor,TCGA-YY,Batch1,IlluminaHiSeq
Sample3,Primary Tumor,TCGA-ZZ,Batch2,IlluminaHiSeq
...
```

Useful for:
- Batch effect analysis
- Sample type filtering
- Quality control

## TCGA Data

### Downloading TCGA Data

The package includes utilities to download TCGA data from UCSC Xena Browser:

```python
from renalprog import dataset

# Download KIRC data
rnaseq, clinical, pheno = dataset.download_data(
    destination='data/raw',
    cancer_type='KIRC',
    remove_gz=True
)
```

### TCGA Data Structure

TCGA data includes:

1. **Gene Expression** (20,531 genes initially)
   - Platform: Illumina HiSeq
   - Normalization: RSEM + upper quartile normalized
   - File: `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`

2. **Clinical Data** (survival and staging)
   - Curated clinical annotations
   - File: `Survival_SupplementalTable_S1_20171025_xena_sp`

3. **Phenotype** (sample metadata)
   - Batch information
   - Sample types (primary tumor, normal, metastatic)
   - File: `TCGA_phenotype_denseDataOnlyDownload.tsv`

### TCGA Sample IDs

TCGA uses barcodes like:
```
TCGA-A3-3308-01A-01R-0864-07
│    │  │    │ │ │  │    │
│    │  │    │ │ │  │    └─ Sequencing center
│    │  │    │ │ │  └────── Plate
│    │  │    │ │ └─────────  RNA/DNA type
│    │  │    │ └───────────  Aliquot
│    │  │    └─────────────  Sample type (01=Primary, 11=Normal)
│    │  └──────────────────  Patient ID
│    └─────────────────────  Tissue source site
└──────────────────────────  Project (TCGA)
```

The package handles these automatically.

## Custom (Non-TCGA) Data

### Preparing Your Own Data

If you have your own gene expression data:

#### Step 1: Format Gene Expression

```python
import pandas as pd

# Load your data
# Assume you have a samples × genes matrix
your_data = pd.read_csv('your_rnaseq_data.csv', index_col=0)

# Ensure proper format
assert your_data.shape[0] < your_data.shape[1], "Transpose if needed"
# Should have more genes than samples

# Check for missing values
assert not your_data.isnull().any().any(), "Remove missing values"

# Save in standard format
your_data.to_csv('data/raw/your_data_formatted.csv')
```

#### Step 2: Format Clinical Data

```python
# Create clinical file
clinical = pd.DataFrame({
    'sample_id': your_data.index,
    'ajcc_pathologic_tumor_stage': your_stages,  # e.g., ["Stage I", "Stage III", ...]
    'age': your_ages,
    'gender': your_gender
})

# Map to early/late if needed
stage_map = {
    'Stage I': 'Early',
    'Stage II': 'Early',
    'Stage III': 'Late',
    'Stage IV': 'Late'
}
clinical['stage_binary'] = clinical['ajcc_pathologic_tumor_stage'].map(stage_map)

clinical.to_csv('data/raw/your_clinical_formatted.csv', index=False)
```

#### Step 3: Validate

```python
from renalprog import dataset

# Try loading
rnaseq = dataset.load_data('data/raw/your_data_formatted.csv')
clinical = dataset.load_data('data/raw/your_clinical_formatted.csv')

# Check alignment
assert set(rnaseq.index) == set(clinical['sample_id']), "Sample IDs must match!"

print(f"Samples: {rnaseq.shape[0]}")
print(f"Genes: {rnaseq.shape[1]}")
print(f"Stages: {clinical['ajcc_pathologic_tumor_stage'].value_counts()}")
```

## Data Preprocessing

### Automatic Preprocessing

The pipeline automatically:

1. **Filters low-expression genes**
   - Removes genes with low mean expression
   - Removes genes with low variance
   - Keeps genes expressed in ≥20% of samples

2. **Detects outliers**
   - Uses Mahalanobis distance
   - Removes samples >3 SD from mean
   - Configurable significance level

3. **Normalizes** (optional)
   - Log transformation
   - Z-score normalization
   - Quantile normalization

### Preprocessing Parameters

```python
from renalprog import features

processed, info = features.preprocess_rnaseq(
    data=rnaseq,
    filter_expression=True,
    mean_threshold=0.5,        # Minimum mean expression
    var_threshold=0.5,         # Minimum variance
    min_sample_fraction=0.2,   # Gene must be in ≥20% samples
    detect_outliers=True,
    alpha=0.05,                # Outlier detection significance
    seed=2023
)

# Check what was filtered
print(f"Original: {rnaseq.shape[1]} genes")
print(f"After filtering: {processed.shape[1]} genes")
print(f"Removed: {info['n_outliers']} outlier samples")
```

## Quality Control

### Check Data Quality

```python
import pandas as pd
import numpy as np

# Load data
rnaseq = pd.read_csv('data/processed/rnaseq.csv', index_col=0)

# 1. Check for missing values
missing = rnaseq.isnull().sum().sum()
print(f"Missing values: {missing}")
assert missing == 0, "Data contains missing values!"

# 2. Check for duplicates
duplicates = rnaseq.index.duplicated().sum()
print(f"Duplicate samples: {duplicates}")
assert duplicates == 0, "Duplicate samples found!"

# 3. Check value range
print(f"Min value: {rnaseq.min().min()}")
print(f"Max value: {rnaseq.max().max()}")
print(f"Mean expression: {rnaseq.mean().mean()}")

# 4. Check distribution
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# Overall distribution
rnaseq.values.flatten().hist(bins=50, ax=axes[0])
axes[0].set_title('Expression Value Distribution')
axes[0].set_xlabel('Expression')

# Per-sample means
rnaseq.mean(axis=1).hist(bins=30, ax=axes[1])
axes[1].set_title('Per-Sample Mean Expression')
axes[1].set_xlabel('Mean Expression')

plt.savefig('reports/figures/data_quality.png')
```

### Sample Size Guidelines

| Analysis Type | Minimum Samples | Recommended | Optimal |
|---------------|----------------|-------------|---------|
| VAE Training | 50 | 200 | 500+ |
| Classification | 100 | 300 | 1000+ |
| Trajectory Generation | 20 per stage | 50 per stage | 100+ per stage |
| Enrichment Analysis | N/A | N/A | N/A (uses trajectories) |

**Note**: More samples generally improve results, but you need balanced stages:

```python
# Check stage balance
clinical['ajcc_pathologic_tumor_stage'].value_counts()

# Should have:
# Early stages: ≥25% of samples
# Late stages: ≥25% of samples
```

## Example Datasets

### KIRC (Kidney Renal Clear Cell Carcinoma)

The main dataset used in publication:

- **Samples**: 534 (498 after filtering)
- **Genes**: 20,531 (5,127 after filtering)
- **Stages**: I (216), II (50), III (120), IV (80)
- **Source**: TCGA via UCSC Xena

### Other TCGA Cancers

The pipeline works with other TCGA cancers:

```python
# Download different cancer type
rnaseq, clinical, pheno = dataset.download_data(
    destination='data/raw',
    cancer_type='BRCA'  # or LUAD, COAD, etc.
)
```

Supported: BRCA, LUAD, LUSC, COAD, STAD, LIHC, PRAD, etc.

## Troubleshooting

### "Sample IDs don't match"

```python
# Find mismatches
rnaseq_ids = set(rnaseq.index)
clinical_ids = set(clinical['sample_id'])

only_in_rnaseq = rnaseq_ids - clinical_ids
only_in_clinical = clinical_ids - rnaseq_ids

print(f"Only in RNA-seq: {len(only_in_rnaseq)}")
print(f"Only in clinical: {len(only_in_clinical)}")

# Keep only matching samples
common_ids = rnaseq_ids & clinical_ids
rnaseq = rnaseq.loc[common_ids]
clinical = clinical[clinical['sample_id'].isin(common_ids)]
```

### "Too few genes after filtering"

```python
# Relax filtering thresholds
processed, info = features.preprocess_rnaseq(
    data=rnaseq,
    mean_threshold=0.1,     # Lower threshold
    var_threshold=0.1,      # Lower threshold
    min_sample_fraction=0.1 # Require expression in only 10% samples
)
```

### "Data looks weird"

```python
# Check if data needs transformation
print("Min:", rnaseq.min().min())
print("Max:", rnaseq.max().max())

# If values are large (>100), likely not log-transformed
if rnaseq.max().max() > 100:
    rnaseq_log = np.log2(rnaseq + 1)
    rnaseq = rnaseq_log
```

## Next Steps

- [Quick Start Tutorial](../tutorials/quickstart.md)
- [Data Processing Pipeline](../tutorials/step1-data-processing.md)
- [API Reference](../api/dataset.md)

