# Enrichment Analysis

Pathway enrichment analysis using PyDESeq2 and GSEA.

## Overview

The enrichment module provides tools for:

1. **Differential Expression Analysis**: Using PyDESeq2 to identify significantly changed genes
2. **Gene Set Enrichment Analysis (GSEA)**: Pathway-level analysis using GSEA CLI
3. **Pathway Visualization**: Heatmap generation for pathway enrichment across trajectories

## Main Classes

::: renalprog.enrichment.EnrichmentPipeline
    options:
      show_source: true
      members:
        - __init__
        - run
        - _run_deseq_processing
        - _run_gsea_commands
        - _combine_gsea_results

## Functions

### Differential Expression

::: renalprog.enrichment.run_deseq2_analysis
    options:
      show_source: true

::: renalprog.enrichment.create_rnk_file
    options:
      show_source: true

### GSEA Analysis

::: renalprog.enrichment.run_gsea_command
    options:
      show_source: true

::: renalprog.enrichment.parse_gsea_results
    options:
      show_source: true

### Visualization

::: renalprog.enrichment.generate_pathway_heatmap
    options:
      show_source: true

::: renalprog.enrichment.plot_heatmap_regulation
    options:
      show_source: true

## Usage Examples

### Running Complete Enrichment Pipeline

```python
from renalprog.enrichment import EnrichmentPipeline
from pathlib import Path

# Initialize pipeline
pipeline = EnrichmentPipeline(
    cancer_type="KIRC",
    trajectory_dir=Path("data/processed/trajectories"),
    output_dir=Path("data/processed/enrichment"),
    gsea_path=Path("./GSEA_4.3.2/gsea-cli.sh"),
    pathways_file=Path("data/external/ReactomePathways.gmt"),
    n_threads=8,
    n_threads_per_deseq=8
)

# Run full pipeline
results = pipeline.run()
```

### Running DESeq2 Analysis

```python
from renalprog.enrichment import run_deseq2_analysis
import pandas as pd

# Load trajectory data
trajectory_data = pd.read_csv("trajectory_001.csv", index_col=0)
control_data = pd.read_csv("control.csv", index_col=0)

# Run DESeq2
results_df = run_deseq2_analysis(
    trajectory_samples=trajectory_data,
    control_samples=control_data,
    n_threads=8
)

# Results contain: log2FoldChange, pvalue, padj, etc.
print(results_df.head())
```

### Creating RNK Files for GSEA

```python
from renalprog.enrichment import create_rnk_file

# Create ranked gene list from DESeq2 results
rnk_file = create_rnk_file(
    deseq_results=results_df,
    output_path="analysis/genes.rnk"
)
```

### Running GSEA

```python
from renalprog.enrichment import run_gsea_command

# Run GSEA on ranked gene list
gsea_output = run_gsea_command(
    rnk_file="analysis/genes.rnk",
    gmt_file="data/external/ReactomePathways.gmt",
    output_dir="analysis/gsea_results",
    label="trajectory_001",
    gsea_path="./GSEA_4.3.2/gsea-cli.sh"
)
```

### Generating Pathway Heatmaps

```python
from renalprog.enrichment import generate_pathway_heatmap

# Generate heatmaps from enrichment results
heatmap_data, figures = generate_pathway_heatmap(
    enrichment_file="data/processed/enrichment/trajectory_enrichment.csv",
    output_dir="data/processed/enrichment",
    fdr_threshold=0.05,
    n_timepoints=50
)

# figures contains:
# - "top_50_changing": Most variable pathways
# - "top_50_upregulated": Most upregulated pathways
# - "top_50_downregulated": Most downregulated pathways
# - "high_level": Reactome high-level pathways
# - "literature": Literature-curated pathways
```

## Configuration

### EnrichmentPipeline Parameters

- `cancer_type`: Cancer type identifier (e.g., "KIRC", "BRCA")
- `trajectory_dir`: Directory containing trajectory CSV files
- `output_dir`: Directory for output files
- `gsea_path`: Path to GSEA CLI executable
- `pathways_file`: Path to GMT file with pathway definitions
- `n_threads`: Number of parallel threads for processing
- `n_threads_per_deseq`: Number of threads per DESeq2 job
- `memory_per_job_gb`: Memory limit per DESeq2 job (default: 12 GB)
- `total_memory_gb`: Total available memory (default: 224 GB)

### GSEA Parameters

- `nperm`: Number of permutations (default: 1000)
- `set_min`: Minimum gene set size (default: 15)
- `set_max`: Maximum gene set size (default: 500)
- `scoring_scheme`: GSEA scoring method (default: "weighted")
- `norm`: Normalization method (default: "meandiv")

## Pathway Collections

### High-Level Reactome Pathways

29 top-level biological processes:
- Autophagy
- Cell Cycle
- DNA Repair
- Immune System
- Metabolism
- Signal Transduction
- And more...

### Literature-Curated Pathways

33 pathways from literature review:
- VHL/HIF pathway
- PI3K/AKT/MTOR pathway
- Warburg effect
- TCA cycle
- And more...

## Output Files

### DESeq2 Results

- `{trajectory_id}_deseq_results.csv`: Complete DESeq2 output
- `{trajectory_id}.rnk`: Ranked gene list for GSEA

### GSEA Results

- `{trajectory_id}/gsea_report_for_na_pos_{timestamp}.tsv`: Positive enrichment
- `{trajectory_id}/gsea_report_for_na_neg_{timestamp}.tsv`: Negative enrichment
- `{trajectory_id}/ranked_gene_list_{timestamp}.tsv`: Ranked genes with scores

### Combined Results

- `trajectory_enrichment.csv`: All GSEA results combined
- `pathway_heatmap_*.png/pdf/svg`: Pathway heatmap visualizations

## Performance Considerations

### Memory Management

DESeq2 analysis is memory-intensive. The pipeline automatically:
- Limits concurrent jobs based on available memory
- Allocates memory per job (default: 12 GB)
- Monitors memory usage

For large datasets:
```python
pipeline = EnrichmentPipeline(
    ...,
    n_threads=8,  # Reduce parallelism
    n_threads_per_deseq=4,  # Reduce threads per job
    memory_per_job_gb=16  # Increase memory per job
)
```

### CPU Utilization

- DESeq2 jobs run in parallel (up to `n_threads`)
- Each job uses `n_threads_per_deseq` threads
- Total CPU usage ≈ `n_threads` × `n_threads_per_deseq`

Recommended settings:
- Small dataset (<100 trajectories): `n_threads=8`, `n_threads_per_deseq=8`
- Large dataset (>100 trajectories): `n_threads=4`, `n_threads_per_deseq=12`

## Troubleshooting

### Out of Memory Errors

```python
# Reduce parallel jobs
pipeline = EnrichmentPipeline(..., n_threads=4)

# Increase memory per job
pipeline = EnrichmentPipeline(..., memory_per_job_gb=20)
```

### GSEA Not Found

Ensure GSEA is installed and path is correct:
```bash
# Test GSEA
./GSEA_4.3.2/gsea-cli.sh --help
```

### No Results Generated

Check logs for errors:
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## See Also

- [Enrichment Analysis Tutorial](../tutorials/step6-enrichment.md)
- [GSEA Installation Guide](../GSEA_INSTALLATION.md)
- [PyDESeq2 Documentation](https://pydeseq2.readthedocs.io/)

