# Step 4: Patient Trajectory Generation

This guide explains how to construct patient trajectories and generate synthetic temporal gene expression data representing disease progression.

## Overview

The trajectory generation pipeline performs the following steps:

1. **Load Data and Models**: Load preprocessed data and trained VAE/Reconstruction Network models
2. **Calculate Distances**: Compute Wasserstein distances between all patient pairs
3. **Link Patients**: Create progression trajectories by linking patients from early to late stages
4. **Build Networks**: Construct trajectory networks showing disease progression paths
5. **Generate Positive Trajectories**: Create synthetic data for disease progression (early→late)
6. **Generate Negative Controls**: Create synthetic data for same-stage transitions (no progression)
7. **Train/Test Split**: Generate separate trajectories for training and testing sets

The pipeline produces:
- **Positive trajectories**: Synthetic gene expression data representing disease progression
- **Negative trajectories**: Control data representing no progression (same stage)
- **Trajectory metadata**: Information about patient connections and progression paths

## Prerequisites

Before running the trajectory generation pipeline, ensure you have:

- **Trained models**: VAE and Reconstruction Network from Step 2
- **Preprocessed data**: From Step 1 data processing
- **Train/test split**: Created in Step 2
- **Python environment**: With PyTorch, pandas, scipy installed
- **Sufficient disk space**: ~1-2 GB for synthetic trajectories

## Usage

### Basic Usage

```bash
python scripts/pipeline_steps/4_trajectories.py
```

This will:
- Load preprocessed KIRC data
- Load trained models
- Calculate patient distances
- Generate progression trajectories
- Create synthetic temporal data
- Save all outputs

## Configuration Parameters

Edit the script to customize parameters:

### Data Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cancer_type` | `"KIRC"` | Cancer type identifier |
| `early_late` | `True` | Use early/late stage classification (vs. I/II/III/IV) |

### Trajectory Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `start_with_first_stage` | `True` | Trajectories must start from early stage |
| `link_next` | `5` | Number of potential next patients to link to |
| `distance_metric` | `"wasserstein"` | Distance metric for comparing patients |
| `n_timepoints` | `50` | Number of timepoints in synthetic trajectories |
| `interpolation_method` | `"linear"` | Interpolation method in latent space |

### Hardware

| Parameter | Default | Description |
|-----------|---------|-------------|
| `force_cpu` | `True` | Force CPU usage even if GPU is available |

## Processing Steps

### Step 1: Load and Preprocess Data

Loads data and metadata:
- Gene expression data (genes × patients)
- Clinical metadata (stage, histological type, race, gender)
- Processes stage information (early: I/II, late: III/IV)

### Step 2: Generate Positive Trajectories

Calculates all possible patient-to-patient transitions:

1. **Distance Calculation**:
   - Computes Wasserstein distance between gene expression distributions
   - Only considers early→late transitions (disease progression)
   - Stores distances for all valid pairs

2. **Patient Linking**:
   - For each patient, finds closest `link_next` patients in next stage
   - Creates connections representing likely progression paths
   - Ensures trajectories start from early stage
   - If race and/or gender are provided as metadata, only links patients with matching attributes, i.e., connections are only made between patients of the same race and gender.

3. **Network Building**:
   - Constructs directed graph of patient transitions
   - Identifies complete trajectories (early→late paths)
   - Removes isolated nodes

### Step 3: Generate Negative Control Trajectories

Creates control trajectories with same-stage transitions:
- Calculates distances within same stage (early→early, late→late)
- Links patients within same stage
- Generates control synthetic data (no progression)

### Step 4: Train/Test Split Trajectories

If train/test split exists from Step 2:
- Creates separate trajectory networks for train and test patients
- Ensures no data leakage between sets
- Generates train-to-train and test-to-test synthetic data

### Step 5: Generate Synthetic Temporal Data

For each trajectory:

1. **Latent Space Interpolation**:
   - Encodes start and end patients to latent space using VAE
   - Performs linear interpolation between latent points
   - Creates `n_timepoints` intermediate representations

2. **Decoding**:
   - Decodes interpolated latent points back to gene expression
   - Applies Reconstruction Network for refinement
   - Normalizes using MinMaxScaler (fitted on train data only)

3. **Save Data**:
   - Each trajectory saved as CSV file
   - Filename format: `{patient1}_to_{patient2}_to_{patientN}.csv`
   - Rows: timepoints, Columns: genes

## Output Structure

```
data/processed/patient_trajectories_KIRC/
├── nodes_metadata.csv                         # Patient clinical metadata
├── all_distances.csv                          # All computed distances
├── random_connections_to_5_next.csv           # Selected patient links
├── trajectories.csv                           # Complete trajectory definitions
├── trajectory_mapping.csv                     # Patient to trajectory mapping
├── random_connections_to_5_next_train.csv     # Train set links
├── random_connections_to_5_next_test.csv      # Test set links
├── negatives_random_connections_to_5_next.csv # Negative control links
└── negatives_trajectories.csv                 # Negative trajectory definitions

data/processed/YYYYMMDD_synthetic_data/kirc/recnet/
├── early_to_late/                             # Positive trajectories
│   ├── train_to_train/
│   │   ├── TCGA-XX-1234_to_TCGA-YY-5678.csv
│   │   └── ...
│   ├── test_to_test/
│   │   ├── TCGA-AA-9999_to_TCGA-BB-0000.csv
│   │   └── ...
│   └── all/
│       └── {patient1}_to_{patient2}.csv       # All trajectories
└── negatives/                                 # Negative control trajectories
    ├── early_to_early/
    └── late_to_late/
```

## Trajectory File Format

Each trajectory CSV file contains:

**Structure**:
- **Rows**: Timepoints (0 to n_timepoints-1)
- **Columns**: Genes (same as preprocessed data)
- **Values**: Normalized gene expression values [0, 1]

**Example**:
```csv
,GENE1,GENE2,GENE3,...
0,0.234,0.567,0.123,...
1,0.245,0.578,0.134,...
2,0.256,0.589,0.145,...
...
49,0.567,0.890,0.456,...
```

## Interpolation Methods

### Linear Interpolation (Default)

Performs straight-line interpolation in latent space:

```
z(t) = (1-t) × z_start + t × z_end
```

Where:
- `t ∈ [0, 1]`: Normalized time
- `z_start`: Latent representation of starting patient
- `z_end`: Latent representation of ending patient

## Understanding the Output

### Trajectory Metadata Files

**nodes_metadata.csv**:
- Patient IDs with clinical information
- Used for filtering and trajectory construction

**all_distances.csv**:
- All computed patient-to-patient distances
- Columns: patient1, patient2, distance, stage1, stage2

**random_connections_to_5_next.csv**:
- Selected patient links for trajectories
- Each patient linked to up to 5 nearest neighbors in next stage

**trajectories.csv**:
- Complete trajectory definitions
- Each row is one trajectory with patient sequence

**trajectory_mapping.csv**:
- Maps each patient to their trajectory ID
- Used for classification in Step 5

### Synthetic Data Files

Each CSV file represents one trajectory:
- Filename shows patient progression path
- 50 timepoints (rows) by default
- Same genes as input data (columns)

## Advanced Usage

### Adjust Trajectory Density

Control number of connections per patient:

```python
# Link each patient to more/fewer neighbors
link_next = 10  # More trajectories
link_next = 3   # Fewer trajectories
```

### Modify Temporal Resolution

Change number of timepoints:

```python
# More detailed trajectories
n_timepoints = 100
```

### Different Stage Classifications

Use 4-stage system instead of early/late:

```python
early_late = False  # Use Stage I, II, III, IV
```

## Biological Interpretation

### What Do Trajectories Represent?

Trajectories represent potential disease progression paths:
- **Starting point**: Patient with early-stage cancer
- **Intermediate points**: Hypothetical gene expression states during progression
- **Ending point**: Patient with late-stage cancer

### Positive vs. Negative Trajectories

**Positive (early→late)**:
- Represent actual disease progression
- Show gene expression changes during cancer advancement
- Used for training progression classifiers

**Negative (same stage)**:
- Control trajectories with no progression
- Help distinguish true progression from random variation
- Important for model specificity

## Next Steps

After generating trajectories:

1. **Inspect outputs**: Check trajectory files and metadata
2. **Visualize trajectories**: Plot gene expression changes over time
3. **Proceed to Step 5**: Train classification models on trajectory data

```bash
python scripts/pipeline_steps/5_classification.py
```
