# R Analysis Scripts Usage Guide

This directory contains production-grade R scripts for enrichment and differential expression analysis. Both scripts follow the same professional standards including comprehensive logging, command-line argument parsing, and extensive documentation.

## Quick Start

### 1. Install Dependencies

```bash
Rscript scripts/r_analysis/install_r_packages.R
```

This will install all required CRAN and Bioconductor packages:
- **CRAN**: gprofiler2, ggplot2, optparse, dplyr, effsize
- **Bioconductor**: limma, edgeR

### 2. Run Differential Expression Analysis

```bash
# With default parameters
Rscript scripts/r_analysis/differential_expression.R

# With custom parameters
Rscript scripts/r_analysis/differential_expression.R \
  --data-path data/my_rnaseq.csv \
  --clinical-path data/my_clinical.csv \
  --alpha 0.05 \
  --description my_analysis
```

### 3. Run Gene Enrichment Analysis

```bash
# With default parameters
Rscript scripts/r_analysis/gene_enrichment.R

# With custom parameters
Rscript scripts/r_analysis/gene_enrichment.R \
  --input data/my_genes.csv \
  --sources GO,KEGG,REAC \
  --top-n 15 \
  --alpha 0.05
```

## Detailed Documentation

### Differential Expression Analysis

**Script:** `differential_expression.R`

**Purpose:** Compares gene expression between two groups (e.g., early vs. late stage tumors) using:
1. **Limma** - Empirical Bayes moderated t-tests (parametric)
2. **Wilcoxon** - Non-parametric rank-sum tests with Cohen's d effect sizes

**Parameters:**
- `--data-path` - RNA-Seq expression matrix (genes × samples)
- `--clinical-path` - Clinical metadata with stage information
- `--output-dir` - Output directory base path
- `--alpha` - Significance threshold (default: 0.01)
- `--group1` - First group name (default: "early")
- `--group2` - Second group name (default: "late")
- `--description` - Output subdirectory suffix

**Output:** Two parallel directories (limma and wilcox) containing:
- `sorted_table.csv` / `wilcox_results.csv` - All genes with statistics
- `significant_genes.csv` - Significant genes only
- `combined_results_all_genes.csv` - All genes from both methods
- `combined_results_both_significant.csv` - Consensus genes (significant in both)
- `summary_statistics.csv` - Summary counts

**Key Features:**
- Automatic clinical stage mapping (Stage I-IV → early/late)
- Auto-transposition detection for data orientation
- Error handling with informative logging
- Consensus identification (genes significant in both methods)
- Comprehensive timestamped logging

### Gene Enrichment Analysis

**Script:** `gene_enrichment.R`

**Purpose:** Functional enrichment analysis using g:Profiler API

**Parameters:**
- `--input` - Gene list (CSV file)
- `--sources` - Databases to query (GO, KEGG, REAC, WP, etc.)
- `--top-n` - Number of top terms per database
- `--output` - Output directory
- `--alpha` - FDR threshold
- `--organism` - Organism identifier

**Output:**
- `gostplot.pdf` / `gostplot.png` - Visualization
- `gostres.csv` / `gostres.tsv` - Complete results
- `plot_data_gprofiler.csv` - Plot-ready data

## Examples

### Example 1: Full Pipeline Run

```bash
# Install dependencies (one-time)
Rscript scripts/r_analysis/install_r_packages.R

# Run differential expression with all defaults
Rscript scripts/r_analysis/differential_expression.R

# Run enrichment on the significant genes
Rscript scripts/r_analysis/gene_enrichment.R \
  --input data/processed/20250117_de_analysis_limma/preprocessed/significant_genes.csv \
  --sources GO,KEGG,REAC \
  --output results/enrichment_analysis
```

### Example 2: Custom Analysis

```bash
# Differential expression with custom parameters
Rscript scripts/r_analysis/differential_expression.R \
  --data-path data/my_rnaseq.csv \
  --clinical-path data/my_clinical.csv \
  --alpha 0.05 \
  --group1 untreated \
  --group2 treated \
  --description batch_corrected

# Enrichment with strict parameters
Rscript scripts/r_analysis/gene_enrichment.R \
  --input data/important_genes.csv \
  --sources GO \
  --top-n 10 \
  --fdr 0.01 \
  --output results/strict_enrichment
```

### Example 3: Mouse Data Analysis

```bash
Rscript scripts/r_analysis/gene_enrichment.R \
  --input data/mouse_genes.csv \
  --sources GO,KEGG \
  --organism mmusculus \
  --top-n 20 \
  --output results/mouse_enrichment
```

## Output File Descriptions

### Differential Expression Files

**sorted_table.csv / wilcox_results.csv**
- One row per gene
- Columns include:
  - `logFC` (Limma) - Log2 fold change
  - `AveExpr` (Limma) - Average expression
  - `adj.P.Val` / `adj.p.value` - Adjusted p-value
  - `effect_size` (Wilcox) - Cohen's d
  - `cohens_d_lower`, `cohens_d_upper` (Wilcox) - 95% CI

**combined_results_both_significant.csv**
- Genes significant in BOTH methods (most robust)
- Useful for downstream validation and prioritization

**summary_statistics.csv**
- Quick overview of significant gene counts
- Per method and consensus

### Gene Enrichment Files

**gostplot.pdf / gostplot.png**
- Bar plot of top enriched terms per database
- P-value on x-axis, term names on y-axis

**gostres.csv / gostres.tsv**
- Complete g:Profiler results
- Includes term IDs, p-values, enrichment statistics

**plot_data_gprofiler.csv**
- Processed data ready for visualization
- Source, term_name, p_value, log_p_value columns

## Troubleshooting

### Package Installation Issues

If `install_r_packages.R` fails:

```bash
# Install Bioconductor manager first
Rscript -e "install.packages('BiocManager')"

# Then install specific packages
Rscript -e "BiocManager::install('limma')"
Rscript -e "BiocManager::install('edgeR')"
```

### Data Format Issues

**Expected data formats:**
- Expression data: CSV with genes as rows, samples as columns
- Clinical data: CSV with row names as sample IDs
- Gene lists: CSV with one gene per line

**Stage mapping:**
- Limma/Wilcox expect: "Stage I", "Stage II", "Stage III", "Stage IV"
- Or pre-mapped: "early", "late"

### Getting Help

View detailed help for any script:

```bash
# The help function is built-in but not directly accessible from CLI
# However, you can examine the documentation in the script:
head -500 scripts/r_analysis/differential_expression.R
```

## Support

For issues, questions, or improvements, refer to the script documentation:
- Help function in each script (comprehensive parameter docs)
- Comments throughout code
- Error messages with logging

---

**Last Updated:** 2025-01-17
**Version:** 1.0.0-alpha
**Status:** Experimental
