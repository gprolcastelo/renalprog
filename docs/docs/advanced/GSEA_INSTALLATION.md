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

1. Visit the GSEA downloads page: <https://www.gsea-msigdb.org/gsea/downloads.jsp>

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

- **Windows**: check Windows build of Java [here](https://learn.microsoft.com/en-us/java/openjdk/install?tabs=exe%2Chomebrew%2Cubuntu).
- **Ubuntu/Debian**: `sudo apt-get install openjdk-11-jdk`

## Pathway Databases

### ReactomePathways.gmt

The renalprog package includes ReactomePathways.gmt in `data/external/`.

### Additional GMT Files

You can download additional gene set databases from MSigDB: <https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp>


Place downloaded GMT files in `data/external/` and reference them with the `--pathways_file` argument.

## Configuration

### Default Configuration

`renalprog` uses these GSEA parameters by default:

* --collapse false Don't collapse probe sets
* --nperm 1000            # Number of permutations
* --set_max 500           # Maximum gene set size
* --set_min 15            # Minimum gene set size




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

**Solution 3**: create a batch file to run GSEA
```cmd
GSEA_4.3.2\gsea-cli.bat --help
```

## License and Citation

### GSEA License

GSEA software is distributed under a custom license:
- **Free for academic and non-profit use**
- **Commercial use requires a license**
- See: <https://www.gsea-msigdb.org/gsea/login.jsp>

###  Citations

For the current implementation of GSEA in `renalprog`, or the package more generally, please see the [citation page](../citation.md). 

Besides, if you use GSEA in your research, please cite:

**GSEA:**
```
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005). 
Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. 
Proceedings of the National Academy of Sciences, 102(43), 15545-15550.
```

If using Reactome pathways:

```
Milacic, M., Beavers, D., Conley, P., Gong, C., Gillespie, M., Griss, J., ... & D’Eustachio, P. (2024). 
The reactome pathway knowledgebase 2024. 
Nucleic acids research, 52(D1), D672-D678.
```

If using other gene sets, please cite accordingly.



## Additional Resources

- **GSEA User Guide**: <https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html>
- **Reactome**: <https://reactome.org/>

## Contact

For renalprog-specific issues, please open an issue on GitHub.

For GSEA software issues, contact the GSEA team through their forum or email.

