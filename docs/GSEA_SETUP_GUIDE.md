# GSEA Installation and Setup Guide

This guide explains how to install and configure GSEA (Gene Set Enrichment Analysis) for use with the renalprog pipeline.

## What is GSEA?

GSEA (Gene Set Enrichment Analysis) is a computational method that determines whether a priori defined sets of genes show statistically significant, concordant differences between two biological states. It's widely used in genomics research to interpret gene expression data.

**Citation:**
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences*, 102(43), 15545-15550.

## Installation

### 1. Download GSEA

1. Visit the GSEA downloads page: https://www.gsea-msigdb.org/gsea/downloads.jsp
2. You will need to register for a free account
3. Download the **GSEA Desktop Application** (Java-based, cross-platform)
4. Choose the appropriate version for your system:
   - **Recommended**: GSEA 4.3.2 or later
   - Available for: Windows, macOS, Linux

### 2. Extract GSEA

Extract the downloaded archive to your project root directory:

```bash
# Linux/macOS
tar -xzf GSEA_Linux_4.3.2.tar.gz

# Windows (use 7-Zip or WinRAR)
# Extract GSEA_Win_4.3.2.zip
```

Your directory structure should look like:
```
renalprog/
├── GSEA_4.3.2/
│   ├── gsea-cli.sh (Linux/macOS)
│   ├── gsea-cli.bat (Windows)
│   ├── gsea.jar
│   └── lib/
├── data/
├── scripts/
└── ...
```

### 3. Install Java

GSEA requires Java 11 or later. Check if Java is installed:

```bash
java -version
```

If Java is not installed:

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get update
sudo apt-get install openjdk-11-jdk
```

**macOS:**
```bash
brew install openjdk@11
```

**Windows:**
Download and install from: https://adoptium.net/

### 4. Make GSEA Executable (Linux/macOS)

```bash
chmod +x GSEA_4.3.2/gsea-cli.sh
```

### 5. Test GSEA Installation

Test that GSEA works:

```bash
# Linux/macOS
./GSEA_4.3.2/gsea-cli.sh --help

# Windows
bash GSEA_4.3.2/gsea-cli.sh --help
```

You should see GSEA help output.

## Pathway Databases

### ReactomePathways.gmt

The pipeline uses Reactome pathway database by default.

**Download:**
1. The file is included in `data/external/ReactomePathways.gmt`
2. If missing, download from: https://reactome.org/download-data

**Citation:**
Gillespie, M., et al. (2022). The reactome pathway knowledgebase 2022. *Nucleic acids research*, 50(D1), D687-D692.

### Alternative: KEGG Pathways

KEGG pathways are also supported (included as `c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt`).

**Citation:**
Kanehisa, M., & Goto, S. (2000). KEGG: kyoto encyclopedia of genes and genomes. *Nucleic acids research*, 28(1), 27-30.

## Configuration

### Default Paths

The pipeline uses these default paths:
- GSEA: `./GSEA_4.3.2/gsea-cli.sh`
- Pathways: `data/external/ReactomePathways.gmt`

### Custom Paths

You can specify custom paths when running the pipeline:

```bash
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --trajectory_dir data/interim/.../trajectories \
    --output_dir data/processed/enrichment \
    --gsea_path /path/to/your/gsea-cli.sh \
    --pathways_file /path/to/your/pathways.gmt
```

## Troubleshooting

### Issue: "GSEA CLI not found"

**Solution:** Ensure GSEA is extracted to the correct location and the path is correct:
```bash
ls -l GSEA_4.3.2/gsea-cli.sh
```

### Issue: "Permission denied"

**Solution (Linux/macOS):**
```bash
chmod +x GSEA_4.3.2/gsea-cli.sh
```

### Issue: "Java not found"

**Solution:** Install Java 11 or later (see Installation section).

### Issue: "bash: command not found" (Windows)

**Solution:** Install Git for Windows (includes Git Bash) or use WSL (Windows Subsystem for Linux).

### Issue: GSEA runs but no results

**Solution:** 
1. Check the GSEA output logs in the output directory
2. Ensure your ranked gene list (.rnk file) is properly formatted
3. Verify pathway file (.gmt) is valid
4. Use the debug script: `python scripts/debug_gsea_output.py --deseq_dir <path>`

## Performance Tuning

### Parallel Execution

The pipeline runs GSEA in parallel. Adjust the number of threads:

```bash
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --n_threads 8  # Use 8 parallel threads
```

**Recommendations:**
- Local computer: 4-8 threads
- Server: 16-32 threads
- Adjust based on available CPU cores and memory

### Memory Requirements

GSEA can be memory-intensive. For large analyses:
- Minimum: 4 GB RAM
- Recommended: 8-16 GB RAM
- Large datasets: 32+ GB RAM

To increase Java heap size, edit `gsea-cli.sh` and modify the `-Xmx` parameter:
```bash
# Change from default
java -Xmx4g ...  # 4 GB

# To larger memory
java -Xmx16g ... # 16 GB
```

## HPC/Cluster Usage

For HPC systems with SLURM or other job schedulers, the pipeline can generate job arrays:

```bash
# Generate GSEA commands without running
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --skip_gsea \
    ...

# Commands are saved in: output_dir/deseq/gsea_commands_*.cmd
```

Then submit as job array (SLURM example):
```bash
sbatch --array=1-N scripts/run_gsea_array.sh
```

## Further Reading

- [GSEA User Guide](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html)
- [MSigDB Database](https://www.gsea-msigdb.org/gsea/msigdb/)
- [Reactome Documentation](https://reactome.org/documentation)

## Support

For GSEA-specific issues, consult:
- [GSEA Forum](https://groups.google.com/g/gsea-help)
- [GSEA Documentation](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html)

For renalprog pipeline issues:
- Check the main README.md
- Review the pipeline documentation in `docs/`

