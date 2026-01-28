# Step 1: Data Processing Pipeline

This guide explains how to download and preprocess TCGA RNA-seq data for the RenalProg pipeline.

## Overview

The data processing pipeline performs the following steps:

1. **Option 1: Download Preprocessed Data** - Download preprocessed data from Zenodo (recommended for reproducibility)
2. **Option 2: Process Raw Data** - Download and process raw TCGA data from scratch

The script provides flexible options for:
- Downloading preprocessed RNAseq data from Zenodo repositories
- Downloading raw TCGA data and processing it locally
- Processing multiple cancer types (KIRC, BRCA)
- Automatic clinical data filtering to match preprocessed samples

## Prerequisites

Before running the data processing pipeline, ensure you have:

- **Python environment**: With required packages installed (pandas, numpy, scipy, etc.)
- **Internet connection**: For downloading data from Zenodo or TCGA
- **Sufficient disk space**: ~2-5 GB for raw data, ~500 MB for preprocessed data

## Usage

### Option 1: Download Preprocessed Data from Zenodo (Recommended)

This is the fastest and most reproducible option. It downloads preprocessed RNAseq data from publicly available Zenodo repositories.

#### Basic Usage - KIRC

```bash
python scripts/pipeline_steps/1_data_processing.py --zenodo_preprocessed
```

This will:
- Download preprocessed KIRC RNAseq data from Zenodo
- Use existing clinical/phenotype files (assumes they're already downloaded)
- Process and filter clinical data to match the RNAseq patients
- Save everything to the appropriate directories

#### With Raw Data Download - KIRC

```bash
python scripts/pipeline_steps/1_data_processing.py --zenodo_preprocessed --download_raw
```

This will:
- Download preprocessed KIRC RNAseq data from Zenodo
- Download raw clinical and phenotype data from TCGA
- Process and filter clinical data to match the RNAseq patients

#### For BRCA Cancer Type

```bash
python scripts/pipeline_steps/1_data_processing.py --zenodo_preprocessed --cancer_type BRCA --download_raw
```

### Option 2: Process Raw TCGA Data from Scratch

This option downloads raw TCGA data and performs all preprocessing steps locally.

#### Process KIRC (assumes raw data already downloaded)

```bash
python scripts/pipeline_steps/1_data_processing.py
```

#### Download and Process KIRC

```bash
python scripts/pipeline_steps/1_data_processing.py --download_raw
```

#### Process BRCA

```bash
python scripts/pipeline_steps/1_data_processing.py --cancer_type BRCA --download_raw
```

## Command-Line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--zenodo_preprocessed` | Flag | False | Download preprocessed data from Zenodo instead of processing raw data |
| `--cancer_type` | String | 'KIRC' | Cancer type to process. Choices: 'KIRC', 'BRCA' |
| `--download_raw` | Flag | False | Download raw TCGA data (otherwise assumes data is already downloaded) |

## Data Sources

### Zenodo Preprocessed Data

#### KIRC
- **Repository**: https://zenodo.org/records/17987300
- **File**: Static_KIRC.csv
- **Direct Download**: https://zenodo.org/records/17987300/files/Static_KIRC.csv?download=1

#### BRCA
- **Repository**: https://zenodo.org/records/17986123
- **File**: Static_BRCA.csv
- **Direct Download**: https://zenodo.org/records/17986123/files/Static_BRCA.csv?download=1

### Raw TCGA Data

When using `--download_raw`, the script downloads:
- RNA-seq expression data from UCSC Xena
- Clinical survival data
- Phenotype data

## Processing Steps

### When Using `--zenodo_preprocessed`

1. **Download preprocessed RNAseq data** from Zenodo
2. **Download clinical/phenotype data** (if `--download_raw` is specified)
3. **Process clinical data**:
   - Load raw clinical and phenotype files
   - Filter for specified cancer type
   - Convert stage information (early/late classification)
   - Filter to match only patients in the preprocessed RNAseq data
4. **Copy control data from repository**:
   - Control samples are included in the repository at `controls/{CANCER_TYPE}/`
   - Copies `rnaseq_controls.csv` and `clinical_controls.csv` to preprocessed data directory
   - No need to download or process controls separately
5. **Save outputs**:
   - `preprocessed_rnaseq.csv`: Preprocessed gene expression data
   - `clinical_data.csv`: Filtered clinical metadata
   - `stages.csv`: Stage information for matched patients
   - `controls/{CANCER_TYPE}_control_rnaseq.csv`: Control RNAseq data (copied from repo)
   - `controls/{CANCER_TYPE}_control_clinical.csv`: Control clinical data (copied from repo)

### When Processing Raw Data

1. **Download raw data** (if `--download_raw` is specified)
2. **Process downloaded data**:
   - Filter for specified cancer type
   - Merge clinical and phenotype information
   - Extract control samples (healthy tissue)
   - Convert stages to early/late classification
3. **Preprocess RNAseq data**:
   - Filter low-expression genes
   - Detect and remove outlier samples using Mahalanobis distance
   - Remove duplicate genes
4. **Process control samples**:
   - Filter control data to match preprocessed genes
5. **Save outputs**:
   - `preprocessed_rnaseq.csv`: Preprocessed gene expression data
   - `preprocessing_info.csv/json`: Metadata about preprocessing steps
   - `stages.csv`: Stage information
   - `controls/preprocessed_control_rnaseq.csv`: Preprocessed control samples

## Output Structure

### Using `--zenodo_preprocessed`

```
data/interim/preprocessed_{CANCER_TYPE}_data/
├── Static_{CANCER_TYPE}.csv                           # Original downloaded file from Zenodo
├── preprocessed_rnaseq.csv                            # Standard format (genes × patients)
├── clinical_data.csv                                  # Filtered clinical data
├── stages.csv                                         # Stage information
└── controls/
    ├── {CANCER_TYPE}_control_rnaseq.csv              # Control RNAseq (copied from repo)
    └── {CANCER_TYPE}_control_clinical.csv            # Control clinical (copied from repo)
```

Also creates (intermediate):
```
data/interim/{CANCER_TYPE}_data/
└── [Processed clinical/phenotype files before filtering]
```

**Note**: Control data is included in the repository at:
```
controls/{CANCER_TYPE}/
├── rnaseq_controls.csv
└── clinical_controls.csv
```

### Processing Raw Data

```
data/raw/
├── EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena
├── Survival_SupplementalTable_S1_20171025_xena_sp
└── TCGA_phenotype_denseDataOnlyDownload.tsv

data/interim/{CANCER_TYPE}_data/
├── {CANCER_TYPE}_rnaseq.csv
├── {CANCER_TYPE}_clinical.csv
└── controls/
    ├── {CANCER_TYPE}_control_rnaseq.csv
    └── {CANCER_TYPE}_control_clinical.csv

data/interim/preprocessed_{CANCER_TYPE}_data/
├── preprocessed_rnaseq.csv
├── preprocessing_info.csv
├── preprocessing_info.json
├── stages.csv
└── controls/
    └── preprocessed_control_rnaseq.csv
```

## Preprocessing Details

When processing raw data, the following steps are applied:

### Gene Expression Filtering
- **Low expression filtering**: Removes genes with mean expression < 0.5 and variance < 0.5
- **Sample fraction**: Genes must be expressed in at least 20% of samples
- **Statistical threshold**: Alpha = 0.05 for statistical tests

### Outlier Detection
- **Method**: Mahalanobis distance with Minimum Covariance Determinant (MCD)
- **Purpose**: Identifies and removes samples with atypical expression patterns
- **Seed**: 2023 (for reproducibility)

### Duplicate Handling
- **Strategy**: Keep first occurrence, remove subsequent duplicates
- **Applied to**: Both gene names and sample IDs

### Stage Classification
- **Early stage**: Stage I and Stage II
- **Late stage**: Stage III and Stage IV

## Complete Example Workflow

### Reproducible Workflow (Zenodo)

For maximum reproducibility using published preprocessed data:

```bash
# Download preprocessed RNAseq and raw clinical data
python scripts/pipeline_steps/1_data_processing.py \
  --zenodo_preprocessed \
  --download_raw \
  --cancer_type KIRC
```

### Custom Workflow (Raw Processing)

For full control over preprocessing parameters:

```bash
# Download and process everything from scratch
python scripts/pipeline_steps/1_data_processing.py \
  --download_raw \
  --cancer_type KIRC
```

## Verification

After running the script, verify the outputs:

```bash
# Check preprocessed data dimensions
python -c "
import pandas as pd
data = pd.read_csv('data/interim/preprocessed_KIRC_data/preprocessed_rnaseq.csv', index_col=0)
clinical = pd.read_csv('data/interim/preprocessed_KIRC_data/clinical_data.csv', index_col=0)
print(f'RNAseq: {data.shape[0]} genes × {data.shape[1]} patients')
print(f'Clinical: {clinical.shape[0]} patients × {clinical.shape[1]} features')
"
```

## Troubleshooting

### Download Failures

If downloads from Zenodo fail:
- Check your internet connection
- Verify the Zenodo URLs are accessible
- Try using a VPN if regional restrictions apply

### Clinical Data Mismatch

If patient IDs don't match between RNAseq and clinical data:
- Ensure you're using the same cancer type for both datasets
- Verify that raw data files are not corrupted
- Check that the preprocessing completed successfully

### Memory Issues

For large datasets:
- Process on a machine with sufficient RAM (recommend 16+ GB)
- Consider processing cancer types separately
- Monitor memory usage during preprocessing

### Missing Raw Data Files

If raw data files are not found:
- Use the `--download_raw` flag to download them automatically
- Or manually download from UCSC Xena and place in `data/raw/`

## Next Steps

After completing the data processing:

1. **Verify outputs**: Check that all expected files were created
2. **Inspect data quality**: Review preprocessing statistics
3. **Proceed to Step 2**: Train VAE models on the preprocessed data

```bash
python scripts/pipeline_steps/2_models.py
```

## Additional Resources

- [TCGA Data Portal](https://portal.gdc.cancer.gov/)
- [UCSC Xena Browser](https://xenabrowser.net/)
- [Zenodo Repository Browser](https://zenodo.org/)
- [Feature Selection Guide](../advanced/ENRICHMENT_ANALYSIS.md)
