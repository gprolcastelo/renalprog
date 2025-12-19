# R Quick Reference for renalprog

## Environment Setup (Recommended)

```bash
# Quick start with conda/mamba
mamba env create -f environment.yml
mamba activate renalprog
pip install -e .

# Or use setup script
bash setup_environment.sh          # Linux/macOS
.\setup_environment.ps1            # Windows PowerShell
```

## Installation

```bash
# Install R packages
Rscript scripts/r_analysis/install_r_packages.R

# Or with Makefile
make install-r

# Install everything (Python + R)
make install-all
```

## Usage

```bash
# Run with defaults
Rscript scripts/r_analysis/gene_enrichment.R

# Run with custom parameters
Rscript scripts/r_analysis/gene_enrichment.R \
  --input data/external/important_genes_shap.csv \
  --sources "GO,KEGG,REAC" \
  --top_n 20

# Using Makefile
make enrichment
```

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--input, -i` | `data/external/important_genes_shap.csv` | Input gene list CSV |
| `--output, -o` | Auto-generated | Output directory |
| `--sources, -s` | `"GO,REAC,KEGG,WP"` | Databases to query |
| `--top_n, -n` | `20` | Top terms per database |
| `--organism` | `"hsapiens"` | Organism (hsapiens, mmusculus, etc.) |
| `--fdr` | `0.05` | FDR threshold |
| `--width` | `25` | Plot width (inches) |
| `--height` | `15` | Plot height (inches) |
| `--dpi` | `600` | Plot resolution |

## Available Databases

- `GO` - Gene Ontology
- `KEGG` - KEGG pathways
- `REAC` - Reactome
- `WP` - WikiPathways
- `TF` - TRANSFAC
- `MIRNA` - miRTarBase
- `HPA` - Human Protein Atlas
- `CORUM` - CORUM complexes
- `HP` - Human Phenotype Ontology

## Output Files

- `gostplot.pdf` - Vector plot
- `gostplot.png` - Raster plot
- `gostres.csv` - Full results (CSV)
- `gostres.tsv` - Full results (TSV)
- `plot_data_gprofiler.csv` - Plot data (CSV)
- `plot_data_gprofiler.tsv` - Plot data (TSV)

## Full Pipeline

```bash
# Run complete pipeline with R analysis
make pipeline-full

# Or step-by-step
python scripts/pipeline_steps/1_data_processing.py
python scripts/pipeline_steps/2_models.py
python scripts/pipeline_steps/3_check_reconstruction.py
python scripts/pipeline_steps/4_trajectories.py
python scripts/pipeline_steps/5_classification.py
Rscript scripts/r_analysis/gene_enrichment.R
```

## Required R Packages

- gprofiler2
- ggplot2
- optparse

## Documentation

- `scripts/r_analysis/README.md` - Detailed usage
- `docs/R_SETUP_GUIDE.md` - Setup & troubleshooting
- `R_SETUP_SUMMARY.md` - Complete summary

