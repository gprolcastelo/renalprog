# GSEA Installation Guide

This guide provides instructions for installing and configuring GSEA (Gene Set Enrichment Analysis) for use with renalprog.

## What is GSEA?

GSEA is a computational method that determines whether a priori defined sets of genes show statistically significant, concordant differences between two biological states.

**Citation:**
```
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005).
Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
Proceedings of the National Academy of Sciences, 102(43), 15545-15550.
DOI: 10.1073/pnas.0506580102
```

## Installation

### Step 1: Download GSEA

1. Visit the GSEA downloads page:
   - https://www.gsea-msigdb.org/gsea/downloads.jsp

2. Register for an account (free for academic use)

3. Download the latest version of GSEA:
   - **For Windows/Mac/Linux**: Download "GSEA_X.X.X.zip" (X.X.X = version number)
   - Recommended version: GSEA 4.3.2 or later

### Step 2: Extract GSEA

1. Extract the downloaded ZIP file to your renalprog project root:
   ```
   renalprog/
   ├── GSEA_4.3.2/          # ← Extract here
   │   ├── gsea-cli.sh      # Unix/Linux/Mac
   │   ├── gsea-cli.bat     # Windows
   │   └── ...
   ├── data/
   ├── renalprog/
   └── ...
   ```

2. The extraction should create a folder named `GSEA_X.X.X/` (e.g., `GSEA_4.3.2/`)

### Step 3: Make Executable (Unix/Linux/Mac)

```bash
cd GSEA_4.3.2
chmod +x gsea-cli.sh
```

### Step 4: Test Installation

#### On Unix/Linux/Mac:
```bash
./GSEA_4.3.2/gsea-cli.sh --help
```

#### On Windows (Git Bash or WSL):
```bash
bash GSEA_4.3.2/gsea-cli.sh --help
```

#### On Windows (Command Prompt):
```cmd
GSEA_4.3.2\gsea-cli.bat --help
```

You should see the GSEA help message if installation was successful.

## System Requirements

### Java

GSEA requires Java 11 or later.

**Check Java version:**
```bash
java -version
```

**If Java is not installed:**

- **Windows**: Download from https://adoptium.net/
- **Mac**: `brew install openjdk@11`
- **Ubuntu/Debian**: `sudo apt-get install openjdk-11-jdk`
- **CentOS/RHEL**: `sudo yum install java-11-openjdk`

### Memory

GSEA can be memory-intensive. Recommended:
- **Minimum**: 4 GB RAM
- **Recommended**: 8 GB RAM per parallel thread
- **Large datasets**: 16+ GB RAM

To increase GSEA memory limit, edit the GSEA configuration file or set Java heap size.

## Pathway Databases

### ReactomePathways.gmt

The renalprog package includes ReactomePathways.gmt in `data/external/`.

**Citation:**
```
Jassal, B., Matthews, L., Viteri, G., Gong, C., Lorente, P., Fabregat, A., ... & D'Eustachio, P. (2020).
The reactome pathway knowledgebase.
Nucleic acids research, 48(D1), D498-D503.
DOI: 10.1093/nar/gkz1031
```

### Additional GMT Files

You can download additional gene set databases from MSigDB:
- https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

Available collections:
- **H**: Hallmark gene sets (50 gene sets)
- **C1**: Positional gene sets (326 gene sets)
- **C2**: Curated gene sets (>6000 gene sets)
  - CGP: Chemical and genetic perturbations
  - CP: Canonical pathways (includes KEGG, Reactome, BioCarta)
- **C3**: Regulatory target gene sets
- **C4**: Computational gene sets
- **C5**: Ontology gene sets (GO)
- **C6**: Oncogenic signature gene sets
- **C7**: Immunologic signature gene sets
- **C8**: Cell type signature gene sets

Place downloaded GMT files in `data/external/` and reference them with the `--pathways_file` argument.

## Configuration

### Default Configuration

Renalprog uses these GSEA parameters by default:

```python
# In renalprog/modeling/enrichment.py
gsea_params = {
    'collapse': 'false',      # Don't collapse probe sets
    'nperm': 1000,            # Number of permutations
    'set_max': 500,           # Maximum gene set size
    'set_min': 15,            # Minimum gene set size
    'scoring_scheme': 'weighted',  # Scoring scheme
    'mode': 'Max_probe',      # Normalization mode
}
```

### Custom Configuration

To modify GSEA parameters, edit `renalprog/modeling/enrichment.py`:

```python
def generate_gsea_command(
    gsea_path: Path,
    rnk_file: Path,
    gmt_file: Path,
    output_dir: Path
) -> str:
    cmd = (
        f'"{gsea_path}" GSEAPreranked '
        f'-gmx "{gmt_file}" '
        f'-rnk "{rnk_file}" '
        f'-out "{output_dir}" '
        f'-collapse false '
        f'-nperm 2000 '        # ← Increase permutations
        f'-set_max 1000 '      # ← Allow larger gene sets
        f'-set_min 10 '        # ← Allow smaller gene sets
        f'-scoring_scheme weighted '
    )
    return cmd
```

## Troubleshooting

### "GSEA CLI not found"

**Problem**: Script cannot find GSEA installation

**Solution**:
1. Verify GSEA is extracted to project root
2. Check the folder name matches `GSEA_X.X.X/`
3. Specify custom path: `--gsea_path /path/to/gsea-cli.sh`

### "Java not found" or "UnsupportedClassVersionError"

**Problem**: Java is not installed or version is too old

**Solution**:
1. Install Java 11 or later
2. Update PATH environment variable
3. Set JAVA_HOME environment variable

### "OutOfMemoryError"

**Problem**: GSEA runs out of memory

**Solution 1**: Increase Java heap size
```bash
# Edit gsea-cli.sh and modify:
java -Xmx8g -jar gsea.jar  # 8 GB heap
```

**Solution 2**: Reduce parallel threads
```bash
# Use fewer threads
python scripts/pipeline_steps/6_enrichment_analysis.py ... --n_threads 2
```

**Solution 3**: Process in batches
```python
# Process subset of trajectories
from pathlib import Path

trajectory_files = list(Path('data/interim/trajectories').glob('*.csv'))

# Process in chunks
chunk_size = 50
for i in range(0, len(trajectory_files), chunk_size):
    chunk = trajectory_files[i:i+chunk_size]
    # Process chunk...
```

### "Permission denied" (Unix/Linux/Mac)

**Problem**: GSEA script is not executable

**Solution**:
```bash
chmod +x GSEA_4.3.2/gsea-cli.sh
```

### Windows-Specific Issues

**Problem**: Cannot run `gsea-cli.sh` on Windows

**Solution 1**: Use Git Bash
```bash
# Install Git for Windows (includes Git Bash)
# Then run:
bash GSEA_4.3.2/gsea-cli.sh --help
```

**Solution 2**: Use WSL (Windows Subsystem for Linux)
```bash
# Install WSL
wsl --install

# Then in WSL terminal:
./GSEA_4.3.2/gsea-cli.sh --help
```

**Solution 3**: Use batch file (if available)
```cmd
GSEA_4.3.2\gsea-cli.bat --help
```

### "No such file or directory" when running GSEA commands

**Problem**: Path contains spaces or special characters

**Solution**: Use quotes in paths
```python
# In enrichment.py, paths are already quoted:
cmd = f'"{gsea_path}" GSEAPreranked -gmx "{gmt_file}" ...'
```

## Alternative: Using GSEApy (Python Implementation)

If you cannot install GSEA CLI, you can modify renalprog to use GSEApy instead:

```bash
pip install gseapy
```

**Note**: Results may differ slightly from official GSEA CLI. The CLI version is considered the gold standard.

To switch to GSEApy, modify `renalprog/modeling/enrichment.py`:

```python
import gseapy as gp

def run_gsea_python(rnk_file, gmt_file, output_dir):
    """Alternative GSEA using Python."""
    pre_res = gp.prerank(
        rnk=rnk_file,
        gene_sets=gmt_file,
        outdir=output_dir,
        permutation_num=1000,
        min_size=15,
        max_size=500,
    )
    return pre_res
```

## License and Citation

### GSEA License

GSEA software is distributed under a custom license:
- **Free for academic and non-profit use**
- **Commercial use requires a license**
- See: https://www.gsea-msigdb.org/gsea/login.jsp

### Required Citations

If you use GSEA in your research, please cite:

**Primary Citation:**
```
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005).
Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
Proceedings of the National Academy of Sciences, 102(43), 15545-15550.
```

**If using MSigDB gene sets:**
```
Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., & Tamayo, P. (2015).
The molecular signatures database hallmark gene set collection.
Cell systems, 1(6), 417-425.
```

**If using Reactome pathways:**
```
Jassal, B., Matthews, L., Viteri, G., Gong, C., Lorente, P., Fabregat, A., ... & D'Eustachio, P. (2020).
The reactome pathway knowledgebase.
Nucleic acids research, 48(D1), D498-D503.
```

## Additional Resources

- **GSEA User Guide**: https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
- **MSigDB User Guide**: https://www.gsea-msigdb.org/gsea/msigdb/
- **GSEA Forum**: https://groups.google.com/g/gsea-help
- **Reactome**: https://reactome.org/
- **GSEApy Documentation**: https://gseapy.readthedocs.io/

## Contact

For renalprog-specific issues, please open an issue on GitHub.

For GSEA software issues, contact the GSEA team through their forum or email.

