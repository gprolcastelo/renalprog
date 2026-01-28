# Dynamic Enrichment Analysis

This document describes the dynamic enrichment analysis pipeline for renalprog.

## Overview

The enrichment analysis pipeline performs Gene Set Enrichment Analysis (GSEA) on synthetic cancer progression trajectories. This allows us to identify biological pathways that are enriched at different stages of progression.

## Pipeline Steps

### 1. pyDESeq2 Differential Expression Analysis

For each synthetic trajectory timepoint:

1. Load trajectory gene expression data (reverse log-transform from RSEM)
2. Load healthy control samples (reverse log-transform from RSEM)
3. Run pyDESeq2 differential expression analysis comparing trajectory vs controls
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

**Technical Details:**

The DESeq2 processing involves:

1. **Reverse log-transformation**: Input data is log-transformed RSEM values, which are converted back to RSEM.
2. **PyDESeq2 analysis**: Proper variance modeling and statistical testing.
3. **Rank file generation**: Genes ranked by log fold-change for GSEA preranked mode.



### GSEA Parameters

Default GSEA parameters in `generate_gsea_command()`:
- `collapse`: false
- `nperm`: 1000 permutations
- `set_max`: 500 (maximum pathway size)
- `set_min`: 15 (minimum pathway size)

## Output Format

### Directory Structure

```
output_dir/
├── test_to_test/                                          # Synthetic trajectories
│   ├── early_to_late/                                     # Transition type
│   │   ├── patient1_to_patient2/                          # Patient trajectory
│   │   │   ├── patient1_to_patient_0.rnk                  # Ranked gene list for timepoint 0
│   │   │   ├── patient1_to_patient_1.rnk
│   │   │   ├── ...
│   │   │   ├── reports/                                   # GSEA output for all patients in directory
│   │   │   │   ├── patient1_to_patient_0.GseaPreranked.*  # GSEA output files
│   │   │   │   ├── gsea_report_for_na_pos_*.tsv           # Positive enrichment report
│   │   │   │   └── gsea_report_for_na_neg_*.tsv           # Negative enrichment report
│   │   ├──patient3_to_patient4/
│   │   ├── ...
│   └── gsea_commands_*.cmd                                # GSEA command files
├── full_gsea_reports_kirc.csv                             # Final combined results
└── heatmap_kirc_significantNES.csv                        # Significant NES heatmap data
```

### Final Results Format

`full_gsea_reports_kirc.csv` columns:

- `Patient`: Patient identifier (e.g., "TCGA-CZ-5989-01_to_TCGA-B0-5108-01)
- `Idx`: Timepoint index (0 to n_samples-1)
- `Transition`: Transition type (e.g., "early_to_late")
- `NAME`: Pathway name (from ReactomePathways.gmt)
- `ES`: Enrichment score
- `NES`: Normalized enrichment score
- `FDR q-val`: False discovery rate q-value

Example:
```csv
Patient,Idx,Transition,NAME,ES,NES,FDR q-val
TCGA-CZ-5989-01_to_TCGA-B0-5108-01,0,early_to_late,Cell Cycle,0.65,2.13,0.001
TCGA-CZ-5989-01_to_TCGA-B0-5108-01,0,early_to_late,DNA Repair,0.52,1.87,0.012
TCGA-CZ-5989-01_to_TCGA-B0-5108-01,1,early_to_late,Cell Cycle,0.71,2.31,0.000
```


### GSEA Not Found

**Error**: `GSEA CLI not found at ./GSEA_4.3.2/gsea-cli.sh`

**Solution**:

1. Download GSEA from https://www.gsea-msigdb.org/gsea/downloads.jsp
2. Extract to project root
3. Or specify custom path with `--gsea_path`

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

## See Also

- [Trajectory Generation](../tutorials/step4-trajectories.md)
- [Classification Pipeline](../tutorials/step5-classification.md)
- [R Analysis Scripts](R_ANALYSIS.md)

