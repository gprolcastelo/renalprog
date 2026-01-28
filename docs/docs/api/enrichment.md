# Enrichment Analysis

Pathway enrichment analysis using PyDESeq2 and GSEA.

## Overview

The `renalprog.enrichment` module provides functions for differential expression analysis and Gene Set Enrichment Analysis (GSEA) to identify biological pathways enriched in cancer progression trajectories.

## Pipeline Usage

For running the complete enrichment pipeline, see the scripts in `scripts/enrichment/`:

- **`py_deseq.py`** - Process trajectories of synthetic patient trajectories (or real ones, if you could ever get these).
- **`py_deseq_real.py`** - Process real patient data.
- **`trajectory_formatting.py`** - Combine GSEA results across multiple patients or trajectories.

### Running the Pipeline

**For synthetic trajectories:**
```bash
python scripts/enrichment/py_deseq.py \
  --cancer_type kirc \
  --traj_dir data/interim/trajectories \
  --source_target_file data/interim/source_target.csv \
  --patient_stage_file data/interim/patient_stages.csv \
  --stage_transition early_to_late \
  --file data/interim/trajectories/early_to_late/patient_001.csv \
  --nperm 1000 \
  --pathway_file data/external/ReactomePathways.gmt
```

**For real patients:**
```bash
python scripts/enrichment/py_deseq_real.py \
  --cancer_type kirc \
  --data_dir data/interim/preprocessed/rnaseq.csv \
  --metadata_dir data/interim/preprocessed/clinical.csv \
  --nperm 1000 \
  --pathway_file data/external/ReactomePathways.gmt
```

### Tutorial

For a complete step-by-step guide, see the [Step 6: Enrichment Analysis Tutorial](../tutorials/step6-enrichment.md).

## API Reference

### Core Functions

#### build_gsea_command

::: renalprog.enrichment.build_gsea_command

**Example:**
```python
from renalprog.enrichment import build_gsea_command

# Basic usage with defaults
cmd = build_gsea_command(
    save_path_rnk_fun="data/interim/analysis/patient_001",
    filename_rnk_fun="patient_001.rnk"
)

# Custom GSEA parameters
cmd = build_gsea_command(
    save_path_rnk_fun="data/interim/analysis/patient_001",
    filename_rnk_fun="patient_001.rnk",
    pathway_file="data/external/c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt",
    nperm=5000,
    rnd_seed="42",
    set_max=800,
    set_min=10
)
print(cmd)
```

---

#### get_rnk_single_patient

::: renalprog.enrichment.get_rnk_single_patient

**Example:**
```python
from renalprog.enrichment import get_rnk_single_patient
import pandas as pd

# Prepare fold change data (from DESeq2 analysis)
fold_change_df = pd.DataFrame({
    'log2FoldChange': [2.5, -1.8, 0.5, 3.2, -2.1]
}, index=['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5'])

# Generate GSEA command for single patient
gsea_cmd = get_rnk_single_patient(
    path_above="data/interim/analysis",
    genes_here=['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5'],
    samples_real_l=fold_change_df,
    pat_i="TCGA-A3-3306",
    index_pat=None,
    foldchange=True,
    nperm=1000
)
```

---

#### fun_apply_deseq

::: renalprog.enrichment.fun_apply_deseq

**Example:**
```python
from renalprog.enrichment import fun_apply_deseq
import pandas as pd
import numpy as np

# Prepare RNA-seq data (log2(RSEM+1) format)
rnaseq_data = pd.DataFrame(
    np.random.uniform(0, 10, (100, 5)),  # 100 genes, 5 samples
    index=[f'GENE{i}' for i in range(100)],
    columns=['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5']
)

# Prepare clinical metadata
clinical_data = pd.DataFrame({
    'ajcc_pathologic_tumor_stage': ['Stage I', 'Stage I', 'Stage III', 'Stage III', 'Stage III']
}, index=['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'])

# Run DESeq2 analysis
results_df = fun_apply_deseq(
    rnaseq_in=rnaseq_data,
    clinical_in=clinical_data
)

# Results contain: log2FoldChange, pvalue, padj, etc.
print(results_df[['log2FoldChange', 'pvalue', 'padj']].head())
```

---

#### fun_single_patient_and_gsea

::: renalprog.enrichment.fun_single_patient_and_gsea

**Example:**
```python
from renalprog.enrichment import fun_single_patient_and_gsea
import pandas as pd

# Load your data
rnaseq_df = pd.read_csv('data/interim/rnaseq.csv', index_col=0)
clinical_controls_df = pd.read_csv('data/processed/controls/clinical_control.csv', index_col=0)
rnaseq_controls_df = pd.read_csv('data/processed/controls/rna_control.csv', index_col=0)
genes = rnaseq_df.columns.tolist()

# Analyze single patient with custom GSEA parameters
gsea_cmd = fun_single_patient_and_gsea(
    patient_here='TCGA-A3-3306',
    patient_stage='Stage II',
    rnaseq_data=rnaseq_df,
    path_above_in='data/interim/analysis',
    clinical_controls=clinical_controls_df,
    rnaseq_controls=rnaseq_controls_df,
    gene_list=genes,
    foldchange=True,
    # GSEA parameters
    pathway_file='data/external/ReactomePathways.gmt',
    nperm=5000,
    rnd_seed='42'
)

# The function generates .rnk files and returns the GSEA command to run
print(gsea_cmd)
```

---

#### fun_synth_single_patient_and_gsea

::: renalprog.enrichment.fun_synth_single_patient_and_gsea

**Example:**
```python
from renalprog.enrichment import fun_synth_single_patient_and_gsea
import pandas as pd

# Load synthetic trajectory data
synth_traj = pd.read_csv('data/interim/trajectories/early_to_late/traj_001.csv', index_col=0)
clinical_controls_df = pd.read_csv('data/processed/controls/clinical_control.csv', index_col=0)
rnaseq_controls_df = pd.read_csv('data/processed/controls/rna_control.csv', index_col=0)
genes = synth_traj.index.tolist()

# Analyze synthetic patient at interpolation point 25
synth_patient_data = pd.DataFrame(synth_traj.iloc[:, 25])  # Column 25
gsea_cmd = fun_synth_single_patient_and_gsea(
    patient_here='TCGA-A3-3306_to_TCGA-A3-3307_25',
    stage_trans_i='early_to_late',
    rnaseq_data=synth_patient_data,
    path_above_in='data/interim/enrichment',
    clinical_controls=clinical_controls_df,
    rnaseq_controls=rnaseq_controls_df,
    gene_list=genes,
    index_pat=25,
    foldchange=True,
    # GSEA parameters
    pathway_file='data/external/ReactomePathways.gmt',
    nperm=1000,
    rnd_seed='timestamp'
)
```

---

### Pathway Visualization

#### generate_pathway_heatmap

::: renalprog.enrichment.generate_pathway_heatmap

**Example:**
```python
from renalprog.enrichment import generate_pathway_heatmap
import pandas as pd

# Load enrichment results (combined GSEA output)
enrichment_df = pd.read_csv(
    'data/processed/enrichment/trajectory_enrichment.csv',
    index_col=0
)

# Generate pathway heatmaps
heatmap_data, figures = generate_pathway_heatmap(
    enrichment_df=enrichment_df,
    output_dir='data/processed/enrichment/heatmaps',
    fdr_threshold=0.05,
    colorbar=True,
    legend=False,
    yticks_fontsize=10,
    show=False
)

# Save figures
for name, fig in figures.items():
    fig.savefig(f'data/processed/enrichment/heatmaps/{name}.png', dpi=300, bbox_inches='tight')

# Inspect top pathways
print(heatmap_data['top_50_changing'].head())
```

---

## GSEA Parameters

All enrichment functions support customizable GSEA parameters. See the [Tutorial Step 6: Enrichment Tutorial](../tutorials/step6-enrichment.md) for detailed information.

### Common Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pathway_file` | `'data/external/ReactomePathways.gmt'` | Path to gene set database |
| `mode` | `'Max_probe'` | Collapse mode for multiple probes |
| `norm` | `'meandiv'` | Normalization method |
| `nperm` | `1000` | Number of permutations |
| `rnd_seed` | `'timestamp'` | Random seed |
| `scoring_scheme` | `'weighted'` | Scoring scheme |
| `set_max` | `500` | Maximum gene set size |
| `set_min` | `15` | Minimum gene set size |



## See Also

- [GSEA Installation Guide](../advanced/GSEA_INSTALLATION.md)
- [Tutorial Step 6: Enrichment Tutorial](../tutorials/step6-enrichment.md)
- [Dynamic Enrichment Analysis Overview](../advanced/ENRICHMENT_ANALYSIS.md)
