# R Analysis Pipeline

This guide provides comprehensive documentation for the R analysis scripts in the renalprog package, covering differential expression analysis and gene enrichment analysis.

## Overview

The R analysis pipeline provides two main capabilities:

1. **Differential Expression Analysis** - Identify genes that are differentially expressed between disease stages
2. **Gene Enrichment Analysis** - Identify biological pathways and processes enriched in gene lists

Both scripts are include a command-line interface with argument parsing, comprehensive logging, and error handling.

## Prerequisites

### R Installation

**Required:** R version 4.0 or higher

Check your R version:
```bash
R --version
```

### Installation Methods

This section summarizes different methods to install the required R packages.
Check the [detailed installation instructions](../tutorials/installation.md).

#### Method 1: Conda/Mamba (Recommended)

If using the renalprog conda environment:

```bash
# Activate environment
mamba activate renalprog

# Install R packages via conda
mamba install -c conda-forge -c bioconda \
    r-base \
    r-gprofiler2 \
    r-ggplot2 \
    r-optparse \
    r-dplyr \
    r-effsize \
    bioconductor-limma \
    bioconductor-edger
```

Advantages:

- All dependencies in one environment
- No admin privileges needed
- Reproducible across systems
- Automatic dependency resolution

#### Method 2: Automated Script

Use the provided installation script:

```bash
# From repository root
Rscript scripts/r_analysis/install_r_packages.R
```

This installs:

- **CRAN packages**: gprofiler2, ggplot2, optparse, dplyr, effsize
- **Bioconductor packages**: limma, edgeR

#### Method 3: Manual Installation

In an R console:

```r
# CRAN packages
install.packages(c("gprofiler2", "ggplot2", "optparse", "dplyr", "effsize"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR"))
```

---

## Script 1: Differential Expression Analysis

### Purpose

Identifies genes that are differentially expressed between two groups (e.g., early vs. late stage cancer) using two complementary statistical methods:

1. **Limma** - Linear Models for Microarray Data (parametric, empirical Bayes)
2. **Wilcoxon** - Non-parametric rank-sum test with Cohen's _d_ effect sizes

### Quick Start

```bash
# Basic usage with default parameters
Rscript scripts/r_analysis/differential_expression.R

# Custom analysis
Rscript scripts/r_analysis/differential_expression.R \
    --data-path data/interim/preprocessed_KIRC_data/preprocessed_rnaseq.csv \
    --clinical-path data/interim/preprocessed_KIRC_data/clinical_data.csv \
    --output-dir results/differential_expression \
    --alpha 0.01 \
    --description KIRC_analysis

# View help
Rscript scripts/r_analysis/differential_expression.R --help
```

### Command-Line Arguments

| Argument | Short | Type | Default | Description |
|----------|-------|------|---------|-------------|
| `--data-path` | `-d` | string | `preprocessed_rnaseq.csv` | RNA-seq expression matrix file |
| `--clinical-path` | `-c` | string | `clinical_data.csv` | Clinical metadata file |
| `--output-dir` | `-o` | string | `results/differential_expression` | Base output directory |
| `--alpha` | `-a` | numeric | 0.01 | Significance threshold (FDR/adjusted p-value) |
| `--group1` | `-g1` | string | "early" | First comparison group name |
| `--group2` | `-g2` | string | "late" | Second comparison group name |
| `--description` | `-x` | string | (timestamp) | Output subdirectory suffix |
| `--help` | `-h` | flag | - | Display help message |

### Input Data Format

#### RNA-seq Expression Matrix

Can be in either orientation (auto-detected):

- **Genes × Samples** (preferred)
- **Samples × Genes** (auto-transposed)

**CSV format:**
```
,Sample1,Sample2,Sample3,...
Gene1,5.23,4.89,6.12,...
Gene2,2.45,2.78,2.63,...
Gene3,8.91,9.23,8.75,...
```

**Requirements:**

- First column: Gene names/IDs
- Remaining columns: Sample expression values
- Values should be normalized (log2-transformed recommended)

#### Clinical Metadata

**CSV format:**
```
sample,stage,other_columns,...
Sample1,early,...
Sample2,late,...
Sample3,early,...
```

**Requirements:**

- Must have a `stage` column containing group labels
- Accepts: "early", "late", "Early", "Late", "Stage I", "Stage II", etc.
- Automatic mapping: Stage I/II → early, Stage III/IV → late

### Output Files

The script creates a timestamped directory structure:

```
results/differential_expression/YYYYMMDD_HHMMSS_description/
├── limma/
│   ├── sorted_table.csv              # All genes with limma statistics
│   ├── significant_genes.csv         # Significant genes (FDR < alpha)
│   └── log.txt                       # Limma analysis log
├── wilcox/
│   ├── wilcox_results.csv            # All genes with Wilcoxon statistics
│   ├── significant_genes.csv         # Significant genes (FDR < alpha)
│   └── log.txt                       # Wilcoxon analysis log
├── combined_results_all_genes.csv            # All genes from both methods
├── combined_results_both_significant.csv     # Consensus genes (both methods)
├── summary_statistics.csv                    # Summary of results
└── differential_expression_analysis.log      # Complete analysis log
```

### Output File Descriptions

#### Limma Results (`limma/sorted_table.csv`)

| Column | Description |
|--------|-------------|
| `gene` | Gene identifier |
| `logFC` | Log2 fold change (group2 vs group1) |
| `AveExpr` | Average expression across all samples |
| `t` | t-statistic |
| `P.Value` | Raw p-value |
| `adj.P.Val` | FDR-adjusted p-value (Benjamini-Hochberg) |
| `B` | Log-odds of differential expression |

#### Wilcoxon Results (`wilcox/wilcox_results.csv`)

| Column | Description |
|--------|-------------|
| `gene` | Gene identifier |
| `p_value` | Raw p-value from Wilcoxon test |
| `fdr` | FDR-adjusted p-value (Benjamini-Hochberg) |
| `cohens_d` | Cohen's d effect size |
| `mean_group1` | Mean expression in group1 |
| `mean_group2` | Mean expression in group2 |
| `mean_diff` | Difference in means (group2 - group1) |

#### Combined Results (`combined_results_both_significant.csv`)

Consensus genes significant in **both** Limma and Wilcoxon tests.

| Column | Description |
|--------|-------------|
| `gene` | Gene identifier |
| `limma_logFC` | Log2 fold change from Limma |
| `limma_adj.P.Val` | FDR from Limma |
| `wilcox_fdr` | FDR from Wilcoxon |
| `wilcox_cohens_d` | Effect size from Wilcoxon |
| `avg_rank` | Average rank across both methods |

---

## Script 2: Gene Enrichment Analysis

### Purpose

Performs functional enrichment analysis to identify biological pathways, processes, and functions that are over-represented in a gene list using the g:Profiler web service.

### Quick Start

```bash
# Basic usage with default gene list
Rscript scripts/r_analysis/gene_enrichment.R

# Custom gene list and parameters
Rscript scripts/r_analysis/gene_enrichment.R \
    --input data/external/my_genes.csv \
    --output reports/enrichment_analysis \
    --sources GO,KEGG,REAC,WP \
    --top-n 20 \
    --alpha 0.05

# View help
Rscript scripts/r_analysis/gene_enrichment.R --help
```

### Command-Line Arguments

| Argument | Short | Type | Default | Description                          |
|----------|-------|------|---------|--------------------------------------|
| `--input` | `-i` | string | `genes.csv` | Input gene list (single-column CSV)  |
| `--output` | `-o` | string | (auto-generated) | Output directory                     |
| `--sources` | `-s` | string | "GO,REAC,KEGG,WP" | Databases to query (comma-separated) |
| `--top-n` | `-n` | integer | 20 | Number of top terms per database     |
| `--organism` | - | string | "hsapiens" | Organism identifier                  |
| `--alpha` | `-a` | numeric | 0.05 | FDR significance threshold           |
| `--width` | `-w` | numeric | 25 | Plot width (inches)                  |
| `--height` | - | numeric | 15 | Plot height (inches)                 |
| `--dpi` | - | integer | 600 | Plot resolution                      |

### Data Sources

| Source | Description |
|--------|-------------|
| `GO` | Gene Ontology (BP, MF, CC) |
| `KEGG` | KEGG Pathways |
| `REAC` | Reactome Pathways | 
| `WP` | WikiPathways |

### Input Data Format

**Gene list CSV:**
```
gene
TP53
VEGFA
HIF1A
VHL
PTEN
...
```

**Requirements:**
- CSV file with a column named "gene" (or first column if no header)
- Gene symbols (e.g., "TP53") or Ensembl IDs
- One gene per row
- No duplicates

**Alternative formats accepted:**
```
# Single column, no header
TP53
VEGFA
HIF1A
```

### Output Files

```
reports/enrichment_analysis/YYYYMMDD_HHMMSS/
├── gostplot.pdf                  # Publication-quality plot (vector)
├── gostplot.png                  # High-resolution raster plot
├── gostres.csv                   # Full enrichment results (CSV)
├── gostres.tsv                   # Full enrichment results (TSV)
├── plot_data_gprofiler.csv       # Plot-ready data (CSV)
├── plot_data_gprofiler.tsv       # Plot-ready data (TSV)
└── gene_enrichment.log           # Analysis log
```

### Output File Descriptions

#### Enrichment Results (`gostres.csv`)

| Column | Description |
|--------|-------------|
| `term_id` | Database term identifier (e.g., GO:0007049) |
| `term_name` | Human-readable term name |
| `p_value` | Raw p-value |
| `adjusted_p_value` | FDR-adjusted p-value |
| `source` | Database source (GO, KEGG, etc.) |
| `term_size` | Total genes in term |
| `query_size` | Genes queried |
| `intersection_size` | Genes found in term |
| `intersection` | List of genes in term |
| `precision` | Intersection / query size |
| `recall` | Intersection / term size |

#### Plot Data (`plot_data_gprofiler.csv`)

Pre-processed data ready for custom plotting:
- Filtered to significant terms (FDR < threshold)
- Sorted by significance
- Top N terms per database

---

## Complete Pipeline Workflow

### Step-by-Step Example

```bash
# Step 1: Install R packages (one time)
Rscript scripts/r_analysis/install_r_packages.R

# Step 2: Run differential expression analysis
Rscript scripts/r_analysis/differential_expression.R \
    --data-path data/interim/preprocessed_KIRC_data/preprocessed_rnaseq.csv \
    --clinical-path data/interim/preprocessed_KIRC_data/clinical_data.csv \
    --alpha 0.01 \
    --description KIRC_early_vs_late

# Step 3: Extract significant genes for enrichment
# Option A: Use limma results
head -n 101 results/differential_expression/*/limma/significant_genes.csv > my_genes.csv

# Option B: Use consensus genes (both methods)
head -n 101 results/differential_expression/*/combined_results_both_significant.csv > my_genes.csv

# Step 4: Run enrichment analysis
Rscript scripts/r_analysis/gene_enrichment.R \
    --input my_genes.csv \
    --sources GO,KEGG,REAC \
    --top-n 15 \
    --alpha 0.05
```

### Integration with Python Pipeline

The R scripts integrate seamlessly with the Python renalprog pipeline if R is installed in the same conda/mamba environment.

---

## Troubleshooting

### Package Installation Issues

#### Problem: BiocManager not found

```r
# Solution: Install BiocManager first
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR"))
```


### Data Format Issues

#### Problem: "No stage column found"

**Solution:** Ensure clinical data has a `stage` column:
```r
# Check column names
clinical <- read.csv("clinical_data.csv")
colnames(clinical)
```



### Output Issues

#### Problem: Plots not generated

**Solutions:**
```r
# Check ggplot2 installation
install.packages("ggplot2")

# Reduce plot dimensions
Rscript scripts/r_analysis/gene_enrichment.R --width 15 --height 10

# Check write permissions
dir.create("output_directory", recursive = TRUE)
```

---

## Citation

If you use these R scripts in your research, please cite:

**For differential expression (Limma):**
> Ritchie, M. E., Phipson, B., Wu, D. I., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015).
> "limma powers differential expression analyses for RNA-sequencing and microarray studies."
> Nucleic Acids Research, 43(7), e47.


**For enrichment analysis (g:Profiler):**
> Kolberg, L., Raudvere, U., Kuzmin, I., Vilo, J., & Peterson, H. (2020).
> "gprofiler2--an R package for gene list functional enrichment analysis and namespace conversion toolset g:Profiler."
> F1000Research, 9, ELIXIR-709.

---

## Additional Resources

### Documentation
- [Limma User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)
- [g:Profiler](https://biit.cs.ut.ee/gprofiler/)
- [edgeR User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

### Related Tutorials
- [Full Pipeline](../tutorials/complete-pipeline.md)
- [Dynamic Enrichment Analysis Tutorial](../tutorials/step6-enrichment.md)
- [Classification Tutorial](../tutorials/step5-classification.md)

### External Links
- [Bioconductor](https://bioconductor.org/)
- [CRAN](https://cran.r-project.org/)