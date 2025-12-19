# Enrichment Analysis Quick Reference

Quick reference guide for trajectory-based enrichment analysis.

## Prerequisites

```bash
# Install PyDESeq2
pip install pydeseq2

# Download GSEA CLI tool
# From: https://www.gsea-msigdb.org/gsea/downloads.jsp
# Extract to project root
```

## Quick Start

### 1. Run Complete Pipeline

```bash
python scripts/pipeline_steps/6_enrichment_analysis.py \
    --trajectory_dir data/processed/20251216_synthetic_data/kirc/early_to_late \
    --output_dir data/processed/20251219_enrichment \
    --cancer_type kirc \
    --n_threads 8 \
    --gsea_path ./GSEA_4.3.2/gsea-cli.sh \
    --pathways_file data/external/ReactomePathways.gmt
```

### 2. Generate Pathway Heatmaps

```bash
python scripts/pipeline_steps/6b_generate_pathway_heatmap.py \
    --enrichment_file data/processed/20251219_enrichment/trajectory_enrichment.csv \
    --output_dir data/processed/20251219_enrichment \
    --fdr_threshold 0.05
```

## Python API

### Basic Usage

```python
from renalprog.enrichment import EnrichmentPipeline

pipeline = EnrichmentPipeline(
    trajectory_dir='data/processed/trajectories',
    output_dir='data/processed/enrichment',
    cancer_type='kirc',
    n_threads=8
)

results = pipeline.run()
```

### Generate Heatmaps

```python
from renalprog.enrichment import generate_pathway_heatmap

heatmap_data, figures = generate_pathway_heatmap(
    enrichment_file='data/processed/enrichment/trajectory_enrichment.csv',
    output_dir='data/processed/enrichment',
    fdr_threshold=0.05
)
```

## Pipeline Steps

1. **DESeq2 Analysis**: Differential expression using PyDESeq2
   - Reverse log-transform RSEM data
   - Compare trajectory vs healthy controls
   - Generate ranked gene lists

2. **GSEA Execution**: Pathway enrichment analysis
   - Run GSEAPreranked in parallel
   - Test against Reactome pathways
   - Calculate enrichment scores and FDR

3. **Results Combination**: Aggregate all results
   - Combine all timepoints and trajectories
   - Filter by FDR threshold
   - Save to CSV

4. **Visualization**: Generate heatmaps
   - Top varying pathways
   - Top upregulated/downregulated
   - Selected high-level and literature pathways

## Output Files

```
output_dir/
├── deseq/                              # DESeq2 results
│   └── early_to_late/
│       └── patient_id/
│           ├── *.rnk                   # Ranked gene lists
│           └── gsea_*/                 # GSEA results
├── trajectory_enrichment.csv           # Combined results
└── pathway_heatmaps/                   # Heatmap visualizations
    ├── top_varying_pathways.pdf
    ├── top_upregulated_pathways.pdf
    ├── top_downregulated_pathways.pdf
    ├── high_level_pathways.pdf
    └── literature_pathways.pdf
```

## Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `n_threads` | Parallel processing threads | 4 |
| `fdr_threshold` | FDR q-value cutoff | 0.05 |
| `gsea_nperm` | GSEA permutations | 1000 |
| `set_min` | Min pathway size | 15 |
| `set_max` | Max pathway size | 500 |

## Common Issues

### Memory errors
- Reduce `n_threads`
- Process in smaller batches

### GSEA not found
- Check `gsea_path` points to correct location
- Ensure execute permissions on gsea-cli.sh

### No significant pathways
- Lower `fdr_threshold`
- Check input data quality
- Verify pathway database matches gene IDs

## See Also

- [Full Documentation](ENRICHMENT_ANALYSIS.md)
- [API Reference](api/enrichment.md)
- [Trajectory Generation](TRAJECTORIES.md)

