# Dynamic Enrichment Analysis

This document describes the dynamic enrichment analysis pipeline for renalprog.

## Overview

The enrichment analysis pipeline performs Gene Set Enrichment Analysis (GSEA) on synthetic cancer progression trajectories. This allows us to identify biological pathways that are enriched at different stages of progression.

## Pipeline Steps

### 1. DESeq2 Differential Expression Analysis

For each synthetic trajectory timepoint:
1. Load trajectory gene expression data (reverse log-transform from RSEM)
2. Load healthy control samples (reverse log-transform from RSEM)
3. Run PyDESeq2 differential expression analysis comparing trajectory vs controls
4. Extract log2 fold-change and adjusted p-values for each gene
5. Rank genes by log2 fold-change
6. Save ranked gene list (`.rnk` file) for GSEA

**Note**: The pipeline uses PyDESeq2 for proper differential expression analysis, not simple fold-change calculations. This ensures statistical rigor and proper handling of count data variance.

### 2. GSEA Analysis

For each ranked gene list:
1. Run GSEA using preranked mode
2. Test against pathway database (ReactomePathways.gmt)
3. Calculate enrichment scores and FDR q-values
4. Generate positive and negative enrichment reports

### 3. Results Combination

Combine all GSEA results into a single dataset:
- One row per (patient, timepoint, pathway)
- Includes enrichment score (ES), normalized ES (NES), and FDR q-value
- Missing pathways filled with NaN values

## Installation Requirements

### Python Dependencies

The enrichment pipeline requires PyDESeq2 for differential expression analysis:

```bash
pip install pydeseq2
```

PyDESeq2 is a Python implementation of the DESeq2 method for differential expression analysis of count data.

**Citation:**
```
Muzellec, B., Telenczuk, M., & Cabeli, V. (2022).
PyDESeq2: a python package for bulk RNA-seq differential expression analysis.
bioRxiv, 2022-12.
```

### GSEA CLI Tool

1. Download GSEA from: https://www.gsea-msigdb.org/gsea/downloads.jsp
2. Extract to project root (creates `GSEA_4.3.2/` directory)
3. Ensure `gsea-cli.sh` (Unix) or `gsea-cli.bat` (Windows) is executable

**Citation:**
```
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005).
Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
Proceedings of the National Academy of Sciences, 102(43), 15545-15550.
```

### Pathway Database

The ReactomePathways.gmt file is included in `data/external/ReactomePathways.gmt`.

**Citation:**
```
Jassal, B., Matthews, L., Viteri, G., Gong, C., Lorente, P., Fabregat, A., ... & D'Eustachio, P. (2020).
The reactome pathway knowledgebase.
Nucleic acids research, 48(D1), D498-D503.
```

## Usage

### Basic Usage

```python
from renalprog.enrichment import EnrichmentPipeline

# Initialize pipeline
pipeline = EnrichmentPipeline(
    trajectory_dir='data/interim/20251216_synthetic_data/kirc/early_to_late',
    output_dir='data/processed/20251217_enrichment',
    cancer_type='kirc',
    n_threads=8,
    n_threads_per_deseq=8,  # Threads per DESeq2 analysis to prevent CPU over-subscription
    gsea_path='./GSEA_4.3.2/gsea-cli.sh',
    pathways_file='data/external/ReactomePathways.gmt'
)

# Run complete pipeline
results = pipeline.run()

# Results are saved to: data/processed/20251217_enrichment/trajectory_enrichment.csv
```

### Parameter Tuning for Parallel Processing

The pipeline supports parallel processing with careful control over CPU usage:

- `n_threads`: Number of trajectory files to process in parallel (default: 4)
- `n_threads_per_deseq`: Number of threads for each PyDESeq2 analysis (default: 8)

**Important**: Total CPU usage = `n_threads × n_threads_per_deseq`. Choose values based on your system:

- **Local PC (8 cores)**: `n_threads=2, n_threads_per_deseq=4` → 8 threads total
- **Workstation (32 cores)**: `n_threads=4, n_threads_per_deseq=8` → 32 threads total  
- **HPC cluster (112 cores)**: `n_threads=14, n_threads_per_deseq=8` → 112 threads total

This prevents CPU over-subscription which can cause NODE_FAIL errors on SLURM clusters.

### Command Line Usage

```bash
# Run complete pipeline
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --trajectory_dir data/interim/20251216_synthetic_data/kirc/early_to_late \
    --output_dir data/processed/20251217_enrichment \
    --cancer_type kirc \
    --n_threads 8 \
    --gsea_path ./GSEA_4.3.2/gsea-cli.sh \
    --pathways_file data/external/ReactomePathways.gmt
```

### Advanced Usage

```python
# Run steps separately
from renalprog.enrichment import (
    process_trajectories_for_deseq,
    run_gsea_parallel,
    combine_gsea_results
)

# Step 1: DESeq2 processing (reverse log-transform + differential expression)
process_trajectories_for_deseq(
    trajectory_dir='data/interim/20251216_synthetic_data/kirc/early_to_late',
    output_dir='data/processed/20251217_enrichment',
    cancer_type='kirc',
    n_threads=8
)

# Step 2: Run GSEA
run_gsea_parallel(
    deseq_dir='data/processed/20251217_enrichment/deseq',
    n_threads=8
)

# Step 3: Combine results
results = combine_gsea_results(
    deseq_dir='data/processed/20251217_enrichment/deseq',
    pathways_file='data/external/ReactomePathways.gmt'
)
```

**Technical Details:**

The DESeq2 processing involves:
1. **Reverse log-transformation**: Input data is log-transformed RSEM values, which are converted back to counts
2. **PyDESeq2 analysis**: Proper variance modeling and statistical testing
3. **Rank file generation**: Genes ranked by log2 fold-change for GSEA preranked mode

### Resume from Checkpoint

If the pipeline fails or is interrupted, you can resume:

```bash
# Skip DESeq if already completed
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --trajectory_dir data/interim/20251216_synthetic_data/kirc/early_to_late \
    --output_dir data/processed/20251217_enrichment \
    --skip_deseq

# Skip GSEA if already completed
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --trajectory_dir data/interim/20251216_synthetic_data/kirc/early_to_late \
    --output_dir data/processed/20251217_enrichment \
    --skip_deseq \
    --skip_gsea
```

## Configuration

### Thread Count

The `n_threads` parameter controls parallelization:
- **Recommended**: Number of CPU cores - 1
- **Minimum**: 1 (sequential processing, very slow)
- **Maximum**: Number of CPU cores
- **Default**: 4

Example for different systems:
```python
# Modern PC (8+ cores)
pipeline = EnrichmentPipeline(..., n_threads=6)

# Older PC (4 cores)
pipeline = EnrichmentPipeline(..., n_threads=2)

# High-performance cluster
pipeline = EnrichmentPipeline(..., n_threads=32)
```

### GSEA Parameters

Default GSEA parameters in `generate_gsea_command()`:
- `collapse`: false (use all genes)
- `nperm`: 1000 permutations
- `set_max`: 500 (maximum pathway size)
- `set_min`: 15 (minimum pathway size)

To modify, edit `renalprog/modeling/enrichment.py`:
```python
def generate_gsea_command(...):
    cmd = (
        f'"{gsea_path}" GSEAPreranked '
        f'-gmx "{gmt_file}" '
        f'-rnk "{rnk_file}" '
        f'-out "{output_dir}" '
        f'-collapse false '
        f'-nperm 2000 '  # Increase permutations
        f'-set_max 1000 '  # Allow larger pathways
        f'-set_min 10'  # Allow smaller pathways
    )
```

## Output Format

### Directory Structure

```
output_dir/
├── deseq/                          # DESeq results
│   ├── early_to_late/             # Transition type
│   │   ├── patient1/              # Patient trajectory
│   │   │   ├── patient1_tp0_foldchange.rnk
│   │   │   ├── patient1_tp1_foldchange.rnk
│   │   │   ├── ...
│   │   │   ├── gsea_tp0/          # GSEA output for timepoint 0
│   │   │   │   ├── gsea_report_for_na_pos_*.tsv
│   │   │   │   └── gsea_report_for_na_neg_*.tsv
│   │   │   └── gsea_tp1/
│   │   └── patient2/
│   └── gsea_commands_*.cmd        # GSEA command files
├── gsea/                           # GSEA working directory
└── trajectory_enrichment.csv       # Final combined results
```

### Final Results Format

`trajectory_enrichment.csv` columns:
- `Patient`: Patient identifier (e.g., "TCGA-3Z-A93Z-01_to_TCGA-A3-A8OW-01")
- `Idx`: Timepoint index (0 to n_samples-1)
- `Transition`: Transition type (e.g., "early_to_late")
- `NAME`: Pathway name (from ReactomePathways.gmt)
- `ES`: Enrichment score
- `NES`: Normalized enrichment score
- `FDR q-val`: False discovery rate q-value

Example:
```csv
Patient,Idx,Transition,NAME,ES,NES,FDR q-val
TCGA-3Z-A93Z-01_to_TCGA-A3-A8OW-01,0,early_to_late,Cell Cycle,0.65,2.13,0.001
TCGA-3Z-A93Z-01_to_TCGA-A3-A8OW-01,0,early_to_late,DNA Repair,0.52,1.87,0.012
TCGA-3Z-A93Z-01_to_TCGA-A3-A8OW-01,1,early_to_late,Cell Cycle,0.71,2.31,0.000
```

## Performance

### Computational Requirements

- **Memory**: ~4 GB per thread
- **Disk Space**: ~10 GB for KIRC dataset (500+ trajectories)
- **Time**: ~2-6 hours for KIRC dataset (depends on threads)

### Benchmarks

Testing on different hardware configurations:

| System | Cores | RAM | Trajectories | Time | Threads Used |
|--------|-------|-----|--------------|------|--------------|
| Modern Desktop (2023) | 16 | 32 GB | 500 | 2 hrs | 12 |
| Laptop (2020) | 8 | 16 GB | 500 | 4 hrs | 6 |
| Workstation (2018) | 4 | 8 GB | 500 | 6 hrs | 3 |
| HPC Cluster | 64 | 256 GB | 500 | 0.5 hrs | 48 |

### Optimization Tips

1. **Use SSD**: GSEA creates many temporary files
2. **Adjust threads**: Set to CPU cores - 1 for best performance
3. **Clean up**: Use `--cleanup` flag to save disk space
4. **Resume**: Use `--skip_deseq` or `--skip_gsea` to resume interrupted runs

## Troubleshooting

### CPU Over-Subscription / NODE_FAIL on HPC

**Error**: `NODE_FAIL` on SLURM, or system reports >500% CPU usage

**Cause**: Too many parallel processes using too many threads each

**Solution**: The pipeline uses two levels of parallelization:
1. `n_threads`: Number of trajectory files processed in parallel
2. `n_threads_per_deseq`: Number of threads for each PyDESeq2 analysis

Total CPU usage = `n_threads × n_threads_per_deseq`

**Fix**: Adjust these parameters based on available cores:

```bash
# For 112-core HPC node
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --n_threads 14 \
    --n_threads_per_deseq 8  # Total: 112 threads

# For 32-core workstation
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --n_threads 4 \
    --n_threads_per_deseq 8  # Total: 32 threads

# For 8-core PC
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --n_threads 2 \
    --n_threads_per_deseq 4  # Total: 8 threads
```

The pipeline sets environment variables (OMP_NUM_THREADS, MKL_NUM_THREADS, etc.) to limit threading per process.

### GSEA Not Found

**Error**: `GSEA CLI not found at ./GSEA_4.3.2/gsea-cli.sh`

**Solution**:
1. Download GSEA from https://www.gsea-msigdb.org/gsea/downloads.jsp
2. Extract to project root
3. Or specify custom path with `--gsea_path`

### Memory Errors

**Error**: `OutOfMemoryError` or system freezes

**Solution**:
1. Reduce `n_threads`
2. Process trajectories in batches
3. Increase system swap space

### GSEA Command Failures

**Error**: GSEA commands fail with non-zero exit code

**Solution**:
1. Check GSEA installation
2. Verify pathway file format (GMT)
3. Check file permissions
4. Review GSEA log files in output directories

### Missing Pathways

**Warning**: Some pathways have all NaN values

**Explanation**: Normal - not all pathways are significant in every sample

### Windows-Specific Issues

**Error**: Cannot run `gsea-cli.sh` on Windows

**Solution**:
1. Install Git Bash or WSL (Windows Subsystem for Linux)
2. Use `gsea-cli.bat` instead if available
3. Or run in WSL environment

## References

1. **GSEA Method**:
   - Subramanian et al. (2005). "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." PNAS 102(43):15545-15550.
   - https://www.gsea-msigdb.org/

2. **PyDESeq2**:
   - Muzellec, B., Telenczuk, M., & Cabeli, V. (2022). "PyDESeq2: a python package for bulk RNA-seq differential expression analysis." bioRxiv, 2022-12.
   - https://github.com/owkin/PyDESeq2

3. **DESeq2 Original Method**:
   - Love, M.I., Huber, W., Anders, S. (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology 15(12):550.

4. **Reactome Pathways**:
   - Jassal et al. (2020). "The reactome pathway knowledgebase." Nucleic Acids Research 48(D1):D498-D503.
   - https://reactome.org/

5. **Original Implementation**:
   - Prol-Castelo, G. (2024). My_BRCA repository
   - Files: `src_deseq_and_gsea_NCSR/py_deseq.py`, `trajectory_analysis.py`

## See Also

- [Classification Pipeline](./CLASSIFICATION.md)
- [Trajectory Generation](./TRAJECTORIES.md)
- [R Analysis Scripts](./R_ANALYSIS.md)

