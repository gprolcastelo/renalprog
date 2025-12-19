# R Analysis Scripts

This directory contains R scripts for gene enrichment analysis as part of the renalprog pipeline.

## Prerequisites

- R version 4.0 or higher
- Internet connection (for installing packages and querying g:Profiler API)

## Setup

### Method 1: Conda/Mamba (Recommended for renalprog users)

If you're using the conda/mamba environment for renalprog:

```bash
# Activate your renalprog environment
mamba activate renalprog

# Install R and required packages
mamba install -c conda-forge r-base r-gprofiler2 r-ggplot2 r-optparse
```

Or use the provided environment file:

```bash
# Create environment with all dependencies (Python + R)
mamba env create -f environment.yml
mamba activate renalprog
```

**Advantages:**
- All dependencies managed in one environment
- Better reproducibility
- No admin/sudo privileges required
- Cross-platform compatibility

### Method 2: R's install.packages()

If you have R installed system-wide or prefer CRAN packages:

```bash
# From the repository root
Rscript scripts/r_analysis/install_r_packages.R
```

Or on Windows PowerShell:

```powershell
Rscript scripts\r_analysis\install_r_packages.R
```

### Method 3: Manual Installation

If the automatic installation fails, you can install packages manually in R:

```r
# Start R
R

# Install CRAN packages
install.packages(c("gprofiler2", "ggplot2", "optparse"))
```

## Required R Packages

| Package | Version | Source | Purpose |
|---------|---------|--------|---------|
| `gprofiler2` | Latest | CRAN | Gene enrichment analysis via g:Profiler API |
| `ggplot2` | Latest | CRAN | Data visualization |
| `optparse` | Latest | CRAN | Command-line argument parsing |

## Scripts

### `gene_enrichment.R`

Performs functional enrichment analysis on gene lists using the g:Profiler web service.

**Usage:**

```bash
# Basic usage with default parameters
Rscript scripts/r_analysis/gene_enrichment.R

# Custom input file and parameters
Rscript scripts/r_analysis/gene_enrichment.R \
  --input data/external/my_genes.csv \
  --output reports/figures/my_enrichment \
  --sources "GO,KEGG,REAC" \
  --top_n 20 \
  --fdr 0.05

# View help
Rscript scripts/r_analysis/gene_enrichment.R --help
```

**Parameters:**

- `--input, -i`: Path to CSV file with gene names (default: `data/external/important_genes_shap.csv`)
- `--output, -o`: Output directory for results (default: auto-generated with date)
- `--sources, -s`: Comma-separated databases to query (default: `"GO,REAC,KEGG,WP"`)
  - Available: GO, KEGG, REAC, WP, TF, MIRNA, HPA, CORUM, HP
- `--top_n, -n`: Number of top terms per database (default: 20)
- `--organism`: Organism identifier (default: `"hsapiens"`)
- `--fdr`: FDR significance threshold (default: 0.05)
- `--width`: Plot width in inches (default: 25)
- `--height`: Plot height in inches (default: 15)
- `--dpi`: Plot resolution (default: 600)

**Output Files:**

The script creates a directory with the following files:

- `gostplot.pdf` - Publication-quality plot
- `gostplot.png` - High-resolution PNG
- `gostres.csv` - Full enrichment results (CSV)
- `gostres.tsv` - Full enrichment results (TSV)
- `plot_data_gprofiler.csv` - Processed data for plotting (CSV)
- `plot_data_gprofiler.tsv` - Processed data for plotting (TSV)

**Example Output:**

```
[2025-12-17 10:30:15] INFO: === Gene Enrichment Analysis Pipeline ===
[2025-12-17 10:30:15] INFO: Reading gene data from: data/external/important_genes_shap.csv
[2025-12-17 10:30:15] SUCCESS: Loaded 150 genes
[2025-12-17 10:30:16] INFO: Starting g:Profiler enrichment analysis...
[2025-12-17 10:30:16] INFO: Using databases: GO, KEGG, REAC, WP
[2025-12-17 10:30:16] INFO: Querying 150 genes against g:Profiler
[2025-12-17 10:30:22] SUCCESS: Found 245 significant terms
[2025-12-17 10:30:22] INFO: Processing top 20 terms for each of 4 sources
[2025-12-17 10:30:22] INFO:   - GO: 20 terms
[2025-12-17 10:30:22] INFO:   - KEGG: 18 terms
[2025-12-17 10:30:22] INFO:   - REAC: 20 terms
[2025-12-17 10:30:22] INFO:   - WP: 15 terms
[2025-12-17 10:30:22] SUCCESS: g:Profiler pipeline completed successfully
[2025-12-17 10:30:23] SUCCESS: Plot created successfully
[2025-12-17 10:30:25] SUCCESS: Saved plot: gostplot.pdf
[2025-12-17 10:30:27] SUCCESS: Saved plot: gostplot.png
[2025-12-17 10:30:27] SUCCESS: Saved results: gostres.csv
[2025-12-17 10:30:27] SUCCESS: === Analysis completed successfully ===
```

## Integration with Python Pipeline

The R enrichment analysis is designed to integrate with the Python-based renalprog pipeline:

1. **Step 1-5**: Run Python pipeline to generate classification models and important genes
2. **Step 6**: Run R enrichment analysis on important genes
3. Results are saved in `reports/figures/` for downstream analysis

## Troubleshooting

### Package Installation Issues

**Problem**: Installation fails for `gprofiler2`

**Solution**: Try installing from GitHub:
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("cran/gprofiler2")
```

### API Connection Issues

**Problem**: g:Profiler query fails with network error

**Solutions**:
- Check internet connection
- Verify firewall isn't blocking R
- Try again later (API might be temporarily down)
- Use a VPN if accessing from restricted network

### Memory Issues

**Problem**: R runs out of memory with large gene lists

**Solutions**:
- Reduce `--top_n` parameter
- Limit number of databases queried with `--sources`
- Increase R memory limit: `R --max-mem-size=8G` (adjust as needed)

### Plot Generation Issues

**Problem**: Plots are not generated or are malformed

**Solutions**:
- Ensure `ggplot2` is properly installed: `install.packages("ggplot2")`
- Try reducing plot dimensions: `--width 15 --height 10`
- Check that output directory has write permissions

## Citation

If you use the enrichment analysis in your research, please cite:

**g:Profiler:**
> Liis Kolberg, Uku Raudvere, Ivan Kuzmin, Jaak Vilo, Hedi Peterson (2020).
> gprofiler2 -- an R package for gene list functional enrichment analysis and
> namespace conversion toolset g:Profiler. F1000Research, 9:ELIXIR-709.

## Additional Resources

- [g:Profiler web interface](https://biit.cs.ut.ee/gprofiler/)
- [gprofiler2 R package documentation](https://cran.r-project.org/web/packages/gprofiler2/)
- [g:Profiler API documentation](https://biit.cs.ut.ee/gprofiler/page/docs)

