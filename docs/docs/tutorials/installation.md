# Installation

`renalprog` is a heavy package. Installation requires about 4.5 GB of disk space and a stable internet connection.

### Prerequisites

- Python 3.9 or higher
- R 4.0+ (for enrichment analysis) - **Can be installed via conda/mamba** (recommended)
- Conda or Mamba (recommended for environment management)
- CUDA-capable GPU (optional, for faster VAE training)

### Recommended: Using Mamba/Conda + uv

#### Quick Setup with environment.yml (Easiest)

```bash
# Clone repository
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog

# Create environment from file (includes Python, R, and all dependencies)
mamba env create -f environment.yml
mamba activate renalprog

# Install the package in editable mode
pip install -e .
```

#### Manual Setup

```bash
# Clone repository
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog

# Create environment with Python 3.9 AND R
mamba create -n renalprog "python==3.9" "r-base>=4.0"
mamba activate renalprog

# Install R packages via conda (recommended for reproducibility)
mamba install -c conda-forge r-gprofiler2 r-ggplot2 r-optparse

# Install uv for faster Python package management
pip install uv

# Install Python package
uv pip install -e .

# Install testing dependencies
uv pip install pytest pytest-cov
```

**Note**: Installing R via conda/mamba ensures all dependencies are managed in the same environment, improving reproducibility.

### Alternative: Using pip + venv

```bash
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -e .
```

### With Development Dependencies

```bash
# Using uv
uv pip install -e ".[dev]"

# Or using pip
pip install -e ".[dev]"
```

### R Dependencies for Enrichment Analysis

The package includes R scripts for gene enrichment analysis. You can install R dependencies in several ways:

#### Option 1: Via Conda/Mamba (Recommended)

If using conda/mamba environment:

```bash
# Activate your environment
mamba activate renalprog

# Install R and packages (if not done during initial setup)
mamba install -c conda-forge r-base r-gprofiler2 r-ggplot2 r-optparse
```

#### Option 2: Via R's install.packages()

If you already have R installed system-wide or prefer CRAN:

```bash
# Install R packages using the provided script
Rscript scripts/r_analysis/install_r_packages.R
```

On Windows PowerShell:
```powershell
Rscript scripts\r_analysis\install_r_packages.R
```

**Required R packages:**
- `r-gprofiler2` / `gprofiler2` - Gene enrichment via g:Profiler API
- `r-ggplot2` / `ggplot2` - Visualization
- `r-optparse` / `optparse` - Command-line parsing

### GSEA Installation (Required for Step 6)

For dynamic enrichment analysis, you need to install GSEA:

1. Download from: https://www.gsea-msigdb.org/gsea/index.jsp
2. Extract to project root (creates `GSEA_4.3.2/` directory)
3. See [GSEA Installation Guide](../advanced/GSEA_INSTALLATION.md) for detailed instructions