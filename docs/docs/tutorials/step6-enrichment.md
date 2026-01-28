# Step 6: Enrichment Analysis Pipeline

This guide explains how to run the enrichment analysis pipeline for trajectory data using DESeq2 and GSEA (Gene Set Enrichment Analysis).

## Overview

The enrichment pipeline performs three main steps:

1. **DESeq2 Analysis**: Performs differential expression analysis on trajectory data
2. **GSEA Analysis**: Runs Gene Set Enrichment Analysis on the DESeq2 results
3. **Trajectory Formatting**: Creates the final trajectory dataset with enrichment results

The pipeline generates "greasy files" (batch job files) that need to be executed manually, allowing you to choose the best execution method for your computing environment.

## Prerequisites

Before running the enrichment pipeline, ensure you have:

- **GSEA installed**: The pipeline expects GSEA CLI to be available at `GSEA_4.3.2/gsea-cli.sh`
  - See [GSEA Installation Guide](../advanced/GSEA_INSTALLATION.md) for setup instructions
- **Python environment**: With required packages (DESeq2, pandas, etc.)
- **Trajectory data**: CSV files containing trajectory information (see [previous step](step5-classification.md)).
- **Preprocessed data**: RNASeq and clinical metadata files.
- **Control data**: Control RNASeq and clinical metadata.
  - **Note**: Control data is included in the repository at `controls/{CANCER_TYPE}/` with files:
    - `rnaseq_controls.csv`
    - `clinical_controls.csv`
  - You can use repo controls by setting `USE_REPO_CONTROLS=true`
  - Alternatively, controls from `data/interim/preprocessed_{CANCER_TYPE}_data/controls/` will be used by default
- **Pathways file**: GMT format pathway file (default: `data/external/ReactomePathways.gmt`).

## Usage

### Basic Usage (with defaults)

For KIRC cancer type with default paths:

```bash
bash scripts/enrichment/pipeline.sh kirc
```

### Custom Configuration

You can override default paths by setting environment variables:

```bash
# Example: Custom trajectory directory
CMD_DIR=/path/to/custom/trajectories bash scripts/enrichment/pipeline.sh kirc

# Example: Use controls from repository root
USE_REPO_CONTROLS=true bash scripts/enrichment/pipeline.sh kirc

# Example: Multiple custom parameters
CMD_DIR=/path/to/trajectories \
DATA_DIR=/path/to/rnaseq.csv \
METADATA_DIR=/path/to/clinical.csv \
PATHWAYS_FILE=/path/to/custom_pathways.gmt \
bash scripts/enrichment/pipeline.sh kirc
```

## Configuration Parameters

### Required Parameter

- **cancer_type**: Cancer type identifier (e.g., `kirc`, `brca`, `luad`)
  - Passed as the first argument to the script

### Optional Parameters (Environment Variables)

| Parameter | Description | Default (KIRC) |
|-----------|-------------|----------------|
| `CMD_DIR` | Directory containing trajectory CSV files | `data/processed/synthetic_data/kirc/recnet/early_to_late/test_to_test/` |
| `ERR_OUT_DIR` | Directory for error/output logs | `data/processed/synthetic_data/kirc/recnet/early_to_late/err_out_kirc_test` |
| `ST_FILE` | Source-target file path | `data/processed/patient_trajectories_KIRC/random_connections_to_5_next_test.csv` |
| `DATA_DIR` | RNASeq data file path | `data/interim/preprocessed_KIRC_data/KIRC_rnaseq.csv` |
| `METADATA_DIR` | Clinical metadata file path | `data/interim/preprocessed_KIRC_data/KIRC_clinical.csv` |
| `CONTROL_DATA` | Control RNASeq data file path | `data/interim/preprocessed_KIRC_data/controls/KIRC_control_rnaseq.csv` |
| `CONTROL_METADATA` | Control clinical metadata file path | `data/interim/preprocessed_KIRC_data/controls/KIRC_control_clinical.csv` |
| `META_NODES` | Nodes metadata file path | `data/processed/patient_trajectories_KIRC/nodes_metadata.csv` |
| `STAGE_TRANSITION` | Stage transition identifier | `early_to_late` |
| `PATHWAYS_FILE` | Pathways GMT file | `data/external/ReactomePathways.gmt` |
| `USE_REPO_CONTROLS` | Use controls from repo root | `false` (uses preprocessed controls by default) |

### Default Paths for Other Cancer Types

For cancer types other than KIRC, the script uses different default paths. You can override these by setting the environment variables explicitly.

## Available GSEA Parameters
GSEA parameters can be passed as arguments to the enrichment functions, making it easy to customize GSEA behavior without modifying configuration files.

All functions that use `build_gsea_command` accept the following optional parameters:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pathway_file` | str | `'data/external/ReactomePathways.gmt'` | Path to the GMT file containing gene sets |
| `mode` | str | `'Max_probe'` | GSEA collapse mode for handling multiple probes per gene |
| `norm` | str | `'meandiv'` | Normalization mode for gene set enrichment scores |
| `nperm` | int | `1000` | Number of permutations for significance testing |
| `rnd_seed` | str | `'timestamp'` | Random seed for reproducibility |
| `scoring_scheme` | str | `'weighted'` | Scoring scheme for enrichment calculation |
| `set_max` | int | `500` | Maximum size of gene sets to include |
| `set_min` | int | `15` | Minimum size of gene sets to include |

## Functions Supporting GSEA Parameters

The following functions now accept GSEA parameters as kwargs:

1. **`build_gsea_command()`** - Core function that builds the GSEA command
2. **`get_rnk_single_patient()`** - Generate GSEA command for a single patient
3. **`fun_single_patient_and_gsea()`** - Perform analysis for a single patient
4. **`fun_synth_single_patient_and_gsea()`** - Perform analysis for synthetic patients

## Pipeline Execution Steps

### Step 1: Generate DESeq2 Greasy File

The script will automatically generate a greasy file containing all DESeq2 commands:

```
scripts/enrichment/greasy_deseq_file_<cancer_type>.txt
```

**You must execute this file manually.** The script will display instructions when this file is generated.

#### Execution Options for DESeq2 Greasy File:

Even though the commands can be run sequentially, this would be very slow for large datasets.
It is recommended to parallelize the execution, and highly recommended to use HPC.

Here are several options to run the DESeq2 greasy file:

**Option 1: Sequential Execution**
```bash
bash scripts/enrichment/greasy_deseq_file_kirc.txt
```

**Option 2: GNU Parallel (Recommended)**
```bash
# Run with 8 parallel jobs
parallel -j 8 < scripts/enrichment/greasy_deseq_file_kirc.txt

# Run with all available cores
parallel < scripts/enrichment/greasy_deseq_file_kirc.txt
```

**Option 3: SLURM Job Scheduler for HPC**
```bash
#!/bin/bash
#SBATCH --job-name=deseq_kirc
#SBATCH --output=logs/deseq_%j.out
#SBATCH --error=logs/deseq_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4

# Load required modules
module load parallel

# Run greasy file with parallel
parallel -j 10 < scripts/enrichment/greasy_deseq_file_kirc.txt
```

### Step 2: Generate GSEA Greasy File

After DESeq2 completes, the script will generate a GSEA greasy file:

```
<CMD_DIR>/greasy_<cancer_type>.sh
```

**You must execute this file manually.** The script will display instructions when this file is generated.

#### Execution Options for GSEA Greasy File:

As mentioned before, sequential execution is possible but slow. Parallel execution is recommended.

**Option 1: Sequential Execution**
```bash
bash data/processed/.../greasy_kirc.sh
```

**Option 2: GNU Parallel (Recommended)**
```bash
# Run with 8 parallel jobs
parallel -j 8 < data/processed/.../greasy_kirc.sh

# Monitor progress with a progress bar
parallel --progress -j 8 < data/processed/.../greasy_kirc.sh
```

**Option 3: SLURM with GNU Parallel**
```bash
#!/bin/bash
#SBATCH --job-name=gsea_kirc
#SBATCH --output=logs/gsea_%j.out
#SBATCH --error=logs/gsea_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1

# Load required modules
module load parallel
module load java  # Required for GSEA

# Export GSEA function
gsea() {
    "$(pwd)/GSEA_4.3.2/gsea-cli.sh" "$@"
}
export -f gsea

# Run greasy file with parallel
parallel -j 8 < data/processed/.../greasy_kirc.sh
```

### Step 3: Final Trajectory Formatting

After GSEA completes, you need to run the trajectory formatting step manually. 
The script will display the exact command to execute with all the correct parameters.

**You must execute this command manually.** The script will display instructions.

#### Execution:

Run the command displayed by the script:

```bash
python scripts/enrichment/trajectory_formatting.py \
  --path_synth <CMD_DIR> \
  --pathways_file <PATHWAYS_FILE> \
  --save_dir <parent_dir> \
  --cancer_type <cancer_type>
```

Replace the placeholders with the actual values shown in the script output. 
This creates the final trajectory dataset with enrichment results.

## Complete Example Workflow

Here's a complete example for KIRC with custom paths:

```bash
# Step 1: Set up environment variables
export CMD_DIR="/data/trajectories/kirc/test"
export DATA_DIR="/data/preprocessed/kirc/rnaseq.csv"
export METADATA_DIR="/data/preprocessed/kirc/clinical.csv"
export CONTROL_DATA="/data/controls/kirc/rnaseq_control.csv"
export CONTROL_METADATA="/data/controls/kirc/clinical_control.csv"
export PATHWAYS_FILE="/data/pathways/ReactomePathways.gmt"

# Step 2: Run the enrichment pipeline script
bash scripts/enrichment/pipeline.sh kirc

# The script will pause after generating the DESeq2 greasy file
# Step 3: Execute the DESeq2 greasy file
parallel -j 8 < scripts/enrichment/greasy_deseq_file_kirc.txt

# The script will then pause after generating the GSEA greasy file
# Step 4: Execute the GSEA greasy file
parallel -j 8 < data/processed/.../greasy_kirc.sh

# Step 5: Execute the trajectory formatting step
# Use the exact command shown by the script output
python scripts/enrichment/trajectory_formatting.py \
  --path_synth /data/trajectories/kirc/test \
  --pathways_file /data/pathways/ReactomePathways.gmt \
  --save_dir /data/trajectories/kirc \
  --cancer_type kirc

# Results will be saved to the parent directory of CMD_DIR
```

## Background Execution

For long-running pipelines, you can run the script in the background:

```bash
nohup bash scripts/enrichment/pipeline.sh kirc > enrichment.log 2>&1 &
```

Monitor progress:
```bash
tail -f enrichment.log
```

## Output Files

### Generated Greasy Files:
- `scripts/enrichment/greasy_deseq_file_<cancer_type>.txt`: DESeq2 commands
- `<CMD_DIR>/greasy_<cancer_type>.sh`: GSEA commands

### Final Output:
- `<parent_dir>/<cancer_type>_trajectory_enrichment.csv`: Final trajectory dataset with enrichment results
- DESeq2 results in `<CMD_DIR>/` directory
- GSEA results in `<CMD_DIR>/` subdirectories


## Troubleshooting

### GSEA Not Found
If you see errors about GSEA not being found:
- Verify GSEA is installed at `GSEA_4.3.2/gsea-cli.sh`
- Check execution permissions: `chmod +x GSEA_4.3.2/gsea-cli.sh`
- See [GSEA Installation Guide](../advanced/GSEA_INSTALLATION.md)

### No CSV Files Found
If the script reports no CSV files in CMD_DIR:
- Verify CMD_DIR contains trajectory CSV files
- Check file permissions
- Ensure the trajectory generation step completed successfully

## Additional Resources

- [GSEA Documentation](https://www.gsea-msigdb.org/gsea/doc/)
- [pyDESeq2 Documentation](https://pydeseq2.readthedocs.io/en/latest/index.html)
- [GNU Parallel Tutorial](https://www.gnu.org/software/parallel/parallel_tutorial.html)
- [Enrichment Analysis Guide](../advanced/ENRICHMENT_ANALYSIS.md)

## Next Steps

After completing the enrichment pipeline:
- Analyze the final trajectory enrichment results
- Visualize pathway changes along trajectories
- Generate pathway heatmaps (see [Pathway Heatmap Guide](../advanced/PATHWAY_HEATMAP.md))
- Perform downstream analysis and interpretation
