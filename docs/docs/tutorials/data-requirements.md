# Data Requirements

Understanding the input data formats and requirements for `renalprog`.

!!! note
    Even though the package has been used with bulk RNA-seq data and taking cancer stages as labels, the VAE-based pipeline is data and label agnostic. You can use this pipeline with other types of tabular data.
    For supervised tasks, you can use a CVAE instead of the VAE. For this, please check the [API reference](../api/models.md).

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

!!! warning
    The datasets used by default in this package are in log2-transformed RSEM normalized counts.
    Before training, the pipeline automatically re-scales the data with the MinMaxScaler from scikit-learn to the [0, 1] range.
    You should not train the models with raw counts.

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

You can specify a different column name for staging when loading data, but you need to pass the name to the functions.
By default, `dataset.load_clinical_data()()` looks for `ajcc_pathologic_tumor_stage`, but can be changed with the `stage_column`parameter.
See the [API reference](../api/dataset.md) for details.

Optional Columns:

- `age`: Patient age at diagnosis
- `gender` or `sex`: M/F or Male/Female

These columns may be used for synthetic trajectory generation. See the [tutorial on synthetic trajectory generation](../tutorials/step4-trajectories.md) for details.



The pipeline automatically groups stages:

- **Early**: Stage I, II
- **Late**: Stage III, IV

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
   - Normalization: log2(RSEM+1) 
   - File: `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`

2. **Clinical Data**
   - Curated clinical annotations, including survival data.
   - File: `Survival_SupplementalTable_S1_20171025_xena_sp`

3. **Phenotype** (sample metadata)
   - Sample types (primary tumor, normal, metastatic)
   - File: `TCGA_phenotype_denseDataOnlyDownload.tsv`

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
    'race': your_races,
    'gender': your_gender
})

# Map to early/late if needed
stage_map = {
    'Stage I': 'early',
    'Stage II': 'early',
    'Stage III': 'late',
    'Stage IV': 'late'
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




## Example Datasets

### KIRC (Kidney Renal Clear Cell Carcinoma)

The main dataset used in publication:

- **Samples**: 533 (530 with stage)
- **Genes**: 20,531 (8,516 after preprocessing for feature filtering)
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

For example: BRCA, LUAD, LUSC, COAD, STAD, LIHC, PRAD, etc.


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

