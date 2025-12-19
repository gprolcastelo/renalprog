# R Setup Guide for renalprog

This guide provides detailed instructions for setting up and using the R components of the renalprog package.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Installation](#installation)
3. [Verification](#verification)
4. [Usage](#usage)
5. [Troubleshooting](#troubleshooting)
6. [Advanced Configuration](#advanced-configuration)

## Prerequisites

### Required Software

- **R**: Version 4.0 or higher
  - Download from: https://cran.r-project.org/
  - Verify installation: `R --version`

- **Rscript**: Command-line interface to R (included with R)
  - Verify: `Rscript --version`

### System Requirements

- **Operating System**: Windows, macOS, or Linux
- **Internet Connection**: Required for:
  - Installing R packages from CRAN
  - Querying g:Profiler API during enrichment analysis
- **Disk Space**: ~500 MB for R packages and dependencies
- **RAM**: At least 4 GB recommended

## Installation

### Option 1: Conda/Mamba Environment (Recommended)

**Best for**: Users of the renalprog package who want everything in one reproducible environment.

#### Quick Start with environment.yml

```bash
# Clone repository
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog

# Create environment with Python, R, and all dependencies
mamba env create -f environment.yml
mamba activate renalprog

# Install the package
pip install -e .
```

#### Manual Conda/Mamba Setup

```bash
# Activate your renalprog environment
mamba activate renalprog

# Install R and required packages
mamba install -c conda-forge r-base r-gprofiler2 r-ggplot2 r-optparse
```

**Advantages:**
- ✅ All dependencies in one environment
- ✅ No admin/sudo privileges required
- ✅ Better reproducibility across systems
- ✅ Automatic dependency resolution
- ✅ Cross-platform compatibility
- ✅ Easy to share and recreate environments

**Verify Installation:**
```bash
# Check R is available
Rscript --version

# Check packages are installed
Rscript -e "library(gprofiler2); library(ggplot2); library(optparse)"
```

### Option 2: Automated Installation Script

**Best for**: Users with system-wide R installation.

Using the provided installation script:

**On Linux/macOS:**
```bash
Rscript scripts/r_analysis/install_r_packages.R
```

**On Windows PowerShell:**
```powershell
Rscript scripts\r_analysis\install_r_packages.R
```

**On Windows Command Prompt:**
```cmd
Rscript scripts\r_analysis\install_r_packages.R
```

**Using Makefile (Linux/macOS):**
```bash
make install-r
```

### Option 3: Manual Installation

**Best for**: Custom setups or troubleshooting.

If the automated script fails, install packages manually:

1. **Start R:**
   ```bash
   R
   ```

2. **Install packages:**
   ```r
   # Set CRAN mirror
   options(repos = c(CRAN = "https://cloud.r-project.org/"))
   
   # Install required packages
   install.packages("gprofiler2")
   install.packages("ggplot2")
   install.packages("optparse")
   
   # Verify installations
   library(gprofiler2)
   library(ggplot2)
   library(optparse)
   ```

3. **Exit R:**
   ```r
   quit(save = "no")
   ```

### Option 4: Using renv (For Reproducibility)

**Best for**: R-centric workflows requiring precise package versions.

For a reproducible environment:

```r
# Install renv if not already installed
install.packages("renv")

# Initialize renv in the project (from R console in project directory)
renv::init()

# Install packages
renv::install("gprofiler2")
renv::install("ggplot2")
renv::install("optparse")

# Create a snapshot
renv::snapshot()
```

## Verification

### Check Package Installation

Run this verification script:

```bash
Rscript -e "
packages <- c('gprofiler2', 'ggplot2', 'optparse')
installed <- sapply(packages, require, character.only = TRUE, quietly = TRUE)
if (all(installed)) {
  cat('✓ All R packages are installed correctly\n')
  cat('Package versions:\n')
  for (pkg in packages) {
    cat(sprintf('  %s: %s\n', pkg, packageVersion(pkg)))
  }
} else {
  cat('✗ Missing packages:\n')
  cat(paste('  -', packages[!installed], collapse = '\n'), '\n')
  quit(status = 1)
}
"
```

### Expected Output

```
✓ All R packages are installed correctly
Package versions:
  gprofiler2: 0.2.3
  ggplot2: 3.4.4
  optparse: 1.7.3
```

## Usage

### Basic Usage

Run gene enrichment analysis with default settings:

```bash
Rscript scripts/r_analysis/gene_enrichment.R
```

### Custom Parameters

```bash
Rscript scripts/r_analysis/gene_enrichment.R \
  --input data/external/important_genes_shap.csv \
  --output reports/figures/my_enrichment_results \
  --sources "GO,KEGG,REAC" \
  --top_n 20 \
  --fdr 0.05
```

### Using Makefile

```bash
# Run enrichment analysis
make enrichment

# Run full pipeline including enrichment
make pipeline-full
```

### Integration with Python Pipeline

```bash
# Step 1-5: Python pipeline
python scripts/pipeline_steps/1_data_processing.py
python scripts/pipeline_steps/2_models.py
python scripts/pipeline_steps/3_check_reconstruction.py
python scripts/pipeline_steps/4_trajectories.py
python scripts/pipeline_steps/5_classification.py

# Step 6: R enrichment analysis
Rscript scripts/r_analysis/gene_enrichment.R \
  --input data/external/important_genes_shap.csv \
  --sources "GO,KEGG,REAC,WP"
```

## Troubleshooting

### Common Issues

#### Issue 1: "Rscript: command not found"

**Cause**: R is not installed or not in system PATH

**Solutions**:
1. Install R from https://cran.r-project.org/
2. Add R to PATH:
   - **Windows**: Add `C:\Program Files\R\R-4.x.x\bin` to system PATH
   - **Linux/macOS**: Usually automatic, or add to `.bashrc`/`.zshrc`:
     ```bash
     export PATH="/usr/local/bin/R:$PATH"
     ```

#### Issue 2: Package Installation Fails

**Cause**: Network issues, permissions, or missing system libraries

**Solutions**:

1. **Check internet connection**
2. **Run R as administrator/sudo** (if permission error):
   ```bash
   sudo R
   # Then install packages
   ```

3. **Install system dependencies** (Linux):
   ```bash
   # Ubuntu/Debian
   sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
   
   # CentOS/RHEL
   sudo yum install libcurl-devel openssl-devel libxml2-devel
   ```

4. **Use different CRAN mirror**:
   ```r
   options(repos = c(CRAN = "https://cran.rstudio.com/"))
   install.packages("gprofiler2")
   ```

#### Issue 3: g:Profiler API Connection Fails

**Cause**: Network restrictions, firewall, or API downtime

**Solutions**:
1. Check internet connection
2. Try from different network
3. Check g:Profiler status: https://biit.cs.ut.ee/gprofiler/
4. Use VPN if accessing from restricted network
5. Increase timeout:
   ```r
   options(timeout = 300)  # 5 minutes
   ```

#### Issue 4: Memory Issues with Large Gene Lists

**Cause**: Insufficient RAM

**Solutions**:
1. Increase R memory limit:
   ```bash
   R --max-mem-size=8G
   ```

2. Reduce parameters:
   ```bash
   Rscript scripts/r_analysis/gene_enrichment.R --top_n 10
   ```

3. Limit databases queried:
   ```bash
   Rscript scripts/r_analysis/gene_enrichment.R --sources "GO"
   ```

#### Issue 5: Plot Generation Fails

**Cause**: Missing graphics libraries or corrupted ggplot2

**Solutions**:

1. **Reinstall ggplot2**:
   ```r
   remove.packages("ggplot2")
   install.packages("ggplot2")
   ```

2. **Install system graphics libraries** (Linux):
   ```bash
   sudo apt-get install libcairo2-dev libxt-dev
   ```

3. **Check output directory permissions**:
   ```bash
   chmod 755 reports/figures
   ```

## Advanced Configuration

### Using Custom CRAN Mirror

Create `~/.Rprofile` with:
```r
options(repos = c(CRAN = "https://cran.rstudio.com/"))
```

### Parallel Package Installation

For faster installation:
```r
install.packages("gprofiler2", Ncpus = 4)
```

### Environment Variables

Set environment variables for R:

**Linux/macOS** (`~/.bashrc` or `~/.zshrc`):
```bash
export R_LIBS_USER="~/R/library"
export R_MAX_NUM_DLLS=500
```

**Windows** (System Environment Variables):
```
R_LIBS_USER=C:\Users\YourName\R\library
R_MAX_NUM_DLLS=500
```

### Using RStudio

You can also run the enrichment script from RStudio:

1. Open RStudio
2. Set working directory to project root
3. Open `scripts/r_analysis/gene_enrichment.R`
4. Modify parameters at the top if needed
5. Run the script

### Batch Processing Multiple Gene Lists

Create a batch script:

```bash
#!/bin/bash
# batch_enrichment.sh

for gene_file in data/external/gene_lists/*.csv; do
    output_name=$(basename "$gene_file" .csv)
    Rscript scripts/r_analysis/gene_enrichment.R \
        --input "$gene_file" \
        --output "reports/figures/enrichment_${output_name}" \
        --sources "GO,KEGG,REAC"
done
```

## Additional Resources

- **R Project**: https://www.r-project.org/
- **CRAN**: https://cran.r-project.org/
- **g:Profiler**: https://biit.cs.ut.ee/gprofiler/
- **gprofiler2 Documentation**: https://cran.r-project.org/web/packages/gprofiler2/
- **ggplot2 Documentation**: https://ggplot2.tidyverse.org/

## Getting Help

If you encounter issues not covered in this guide:

1. Check the error message carefully
2. Search for the error on:
   - Stack Overflow: https://stackoverflow.com/questions/tagged/r
   - R-help mailing list archives
3. Check package documentation: `?gprofiler2::gost`
4. Open an issue on the GitHub repository

## Version Information

This guide was written for:
- R >= 4.0
- gprofiler2 >= 0.2.0
- ggplot2 >= 3.0.0
- optparse >= 1.7.0

Check your versions:
```r
sessionInfo()
```

