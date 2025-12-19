# Renalprog Package Migration Plan

**Project**: Migrate My_BRCA repository to renalprog with cookiecutter data science structure
**Target Package Name**: renalprog
**Description**: A Python package for simulating kidney cancer progression with synthetic data generation and machine learning.
**License**: Apache 2.0
**Deadline**: December 12, 2025
**Current Date**: December 17, 2025 (5 days after deadline - extending timeline)
**Status**: Days 1-3 complete, Day 6 complete, Day 7 partially complete

## Progress Summary

- ✅ **Day 1**: Core structure, config, utils (COMPLETE)
- ✅ **Day 2**: VAE models and training (COMPLETE)
- ✅ **Day 3**: Trajectory generation (COMPLETE - Dec 16)
- ✅ **Day 6**: Dynamic enrichment analysis (COMPLETE - Dec 17)
- ✅ **Day 7**: Classification pipeline (COMPLETE - Dec 16)
- ⏭️ **Next**: Days 8-10 remaining tasks

## Timeline

- **Days 1-3 (Nov 27-29)**: Core structure + configuration
- **Days 4-7 (Nov 30-Dec 3)**: Pipeline steps 1-6 migration
- **Days 8-10 (Dec 4-6)**: Pipeline steps 7-9 orchestration
- **Days 11-12 (Dec 7-8)**: Testing implementation
- **Days 13-14 (Dec 9-10)**: Documentation + validation
- **Day 15 (Dec 11)**: Buffer for issues

## Repository Structure (Target)

```
renalprog/
├── LICENSE                      <- Apache 2.0 license
├── Makefile                     <- Convenience commands
├── README.md                    <- Top-level README
├── pyproject.toml              <- Package metadata and tool configuration
├── setup.cfg                    <- flake8 configuration
├── requirements.txt             <- Dependencies
├── data/
│   ├── external/               <- Third party data (gene lists, pathways)
│   ├── interim/                <- Intermediate transformed data
│   ├── processed/              <- Final canonical datasets
│   └── raw/                    <- Original immutable data (TCGA-KIRC)
├── docs/                       <- mkdocs documentation
├── models/                     <- Trained models and predictions
├── notebooks/                  <- Jupyter notebooks for exploration
├── references/                 <- Data dictionaries and manuals
├── reports/                    <- Generated analysis
│   └── figures/               <- Generated graphics
├── scripts/                    <- Standalone scripts
│   ├── r_analysis/            <- R scripts for DESeq2, clusterProfiler
│   └── pipeline_steps/        <- Individual pipeline step scripts
├── tests/                      <- Test suite
│   ├── test_dataset.py
│   ├── test_features.py
│   ├── test_modeling.py
│   └── fixtures/              <- Small test data
└── renalprog/                  <- Source code package
    ├── __init__.py
    ├── config.py               <- Configuration and paths
    ├── dataset.py              <- Data loading and splitting
    ├── features.py             <- Feature engineering and preprocessing
    ├── plots.py                <- Visualization functions
    ├── modeling/
    │   ├── __init__.py
    │   ├── train.py           <- Model training (VAE, classifiers)
    │   └── predict.py         <- Inference and trajectory generation
    └── utils/
        ├── __init__.py
        └── common.py          <- Shared utility functions
```

## Pipeline Steps Mapping

### Step 0: Data Download (Manual)
**Current**: Manual download from Xena Browser
**Target**: Document in README.md with instructions
**Files**: None (user manual step)

### Step 1: Preprocessing
**Current**: `notebooks/1_Preprocessing.ipynb`, `maha_outliers.py`
**Target**: `renalprog/features.py` - functions `filter_low_expression()`, `detect_outliers_mahalanobis()`
**Data**: `data/raw/` → `data/interim/YYYYMMDD_preprocessed_KIRC/`

### Step 1.1: Train/Test Split
**Current**: `src/train_test_split.py`
**Target**: `renalprog/dataset.py` - function `create_train_test_split()`
**Data**: Saves to `data/interim/YYYYMMDD_train_test_split/`

### Step 2: Train VAE
**Current**: `src/python_VAE.py`, `src/models/my_model.py`, `src/models/train_model.py`
**Target**: `renalprog/modeling/train.py` - classes `VAE`, `CVAE`, function `train_vae()`
**Data**: Saves model to `models/YYYYMMDD_VAE_KIRC/`

### Step 3: Adjust Reconstruction
**Current**: `src/adjust_reconstruction.py`
**Target**: `renalprog/modeling/train.py` - class `NetworkReconstruction`, function `train_postprocessing_network()`
**Data**: Saves model to `models/YYYYMMDD_VAE_KIRC/postprocessing.pth`

### Step 4: Assess Reconstruction
**Current**: `src_deseq_and_gsea_NCSR/sdanalyses.py`
**Target**: `renalprog/modeling/predict.py` - function `evaluate_reconstruction()`
**Data**: Saves metrics to `data/processed/YYYYMMDD_reconstruction_metrics/`

### Step 5: Create Patient Connections
**Current**: `notebooks/4_1_trajectories.ipynb` (connection logic using deepgraph/networkx)
**Target**: `renalprog/modeling/predict.py` - function `create_patient_connections()`
**Data**: Saves to `data/interim/YYYYMMDD_patient_connections.csv`

### Step 6: Trajectory Generation
**Current**: 
- Forward: `src_deseq_and_gsea_NCSR/synthetic_data_generation.py`
- Controls: `src_deseq_and_gsea_NCSR/create_controls.py`, `src_deseq_and_gsea_NCSR/generate_control_trajectories.py`
**Target**: `renalprog/modeling/predict.py` - functions `generate_trajectories()`, `generate_control_trajectories()`
**Data**: Saves to `data/interim/YYYYMMDD_synthetic_data/kirc/`

### Step 7: Classification
**Current**: `notebooks/kirc_classification_trajectory.ipynb`, `src/xgb_classifier.py`
**Target**: `renalprog/modeling/train.py` - functions `train_stage_classifier()`, `apply_classifier_to_trajectories()`
**Data**: Saves to `models/YYYYMMDD_kirc_classification/` and predictions to `data/processed/`

### Step 8: Apply Classifier to Trajectories
**Current**: `notebooks/kirc_classification_trajectory.ipynb`
**Target**: `renalprog/modeling/predict.py` - function `classify_trajectories()`
**Data**: Saves to `data/processed/YYYYMMDD_trajectory_classifications.csv`

### Step 9: Gene Clustering
**Current**: `src_deseq_and_gsea_NCSR/tsfresh_and_clustering.py`
**Target**: `renalprog/features.py` - function `cluster_genes_tsfresh()`
**Data**: Saves to `data/interim/YYYYMMDD_tsfresh/`

### Step 10: Static Enrichment
**Current**: `notebooks/clusterProfiler.ipynb` (R)
**Target**: `scripts/r_analysis/cluster_profiler.R` + orchestration in `renalprog/modeling/predict.py`
**Data**: Saves to `data/processed/YYYYMMDD_static_enrichment/`

### Step 11: Dynamic Enrichment
**Current**: `src_deseq_and_gsea_NCSR/full_bash.sh`, `py_deseq.py`, `trajectory_analysis.py`, GSEA CLI
**Target**: `renalprog/modeling/predict.py` - function `dynamic_enrichment_analysis()` using gseapy or keep bash
**Data**: Saves to `data/processed/YYYYMMDD_dynamic_enrichment/`

### Step 12: Differential Gene Expression
**Current**: `src_deseq_and_gsea_NCSR/differential_expression.R`
**Target**: `scripts/r_analysis/differential_expression.R` + orchestration
**Data**: Saves to `data/processed/YYYYMMDD_limma/` and `YYYYMMDD_wilcox/`

### Step 13: Generate Paper Figures
**Current**: `paper_figures.ipynb`, `supplementary.ipynb`
**Target**: `renalprog/plots.py` - functions for each figure type
**Notebooks**: Keep as `notebooks/paper_figures.ipynb` that imports from `renalprog.plots`

## Key Files to Migrate

### Source Code
- `src/models/my_model.py` → `renalprog/modeling/train.py` (VAE, CVAE, AE classes)
- `src/models/train_model.py` → `renalprog/modeling/train.py` (training loops, split_data, set_seed)
- `src/python_VAE.py` → `renalprog/modeling/train.py` (VAE training script)
- `src/adjust_reconstruction.py` → `renalprog/modeling/train.py` (postprocessing network)
- `src/train_test_split.py` → `renalprog/dataset.py`
- `src/xgb_classifier.py` → `renalprog/modeling/train.py` (XGBoost classifier)
- `src/data/fun_interpol.py` → `renalprog/modeling/predict.py` (interpolation functions)
- `src/fun_utils.py` → `renalprog/utils/common.py`
- `maha_outliers.py` → `renalprog/features.py`
- `src_deseq_and_gsea_NCSR/synthetic_data_generation.py` → `renalprog/modeling/predict.py`
- `src_deseq_and_gsea_NCSR/sdanalyses.py` → `renalprog/modeling/predict.py`
- `src_deseq_and_gsea_NCSR/tsfresh_and_clustering.py` → `renalprog/features.py`
- `src_deseq_and_gsea_NCSR/py_deseq.py` → `scripts/pipeline_steps/py_deseq.py`
- `src_deseq_and_gsea_NCSR/trajectory_analysis.py` → `scripts/pipeline_steps/trajectory_analysis.py`
- `src_deseq_and_gsea_NCSR/differential_expression.R` → `scripts/r_analysis/differential_expression.R`

### Data Files
- `data/external/50_important_genes_for_KIRC.csv` → keep
- `data/external/80_important_genes_for_KIRC.csv` → keep
- `data/external/final_genes_classification.csv` → keep
- `data/external/ReactomePathways.gmt` → keep
- `data/external/biomart_human_genes_GRCh37.p13.txt` → keep
- `data/external/c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt` → keep

### Configuration Files
- `requirements.txt` → update and keep
- `setup.py` → migrate to `pyproject.toml`
- `LICENSE` → keep (already Apache 2.0)
- `Makefile` → update for new structure

## Implementation Tasks

### Day 1 (Nov 27) - Core Structure ✅ COMPLETE
- [x] Create migration plan document
- [x] Create `renalprog/` package structure with all `__init__.py` files
- [x] Create `renalprog/config.py` with path management and hyperparameters
- [x] Create `renalprog/utils/common.py` with `set_seed()` and utility functions
- [x] Create `pyproject.toml` with package metadata
- [x] Update `setup.cfg` for flake8
- [x] Implement `renalprog/dataset.py` (AHEAD OF SCHEDULE)
- [x] Implement `renalprog/features.py` (AHEAD OF SCHEDULE)
- [x] Implement `renalprog/plots.py` (AHEAD OF SCHEDULE)
- [x] Create comprehensive test suite (AHEAD OF SCHEDULE)
- [x] Create all documentation files (AHEAD OF SCHEDULE)

### Day 2 (Nov 28) - VAE Model Migration ✅ COMPLETE
**Primary Goal**: Migrate VAE models and training implementation

- [x] Migrate VAE class from `src/models/my_model.py` to `renalprog/modeling/train.py`
- [x] Migrate CVAE class from `src/models/my_model.py` to `renalprog/modeling/train.py`
- [x] Migrate AE class from `src/models/my_model.py` to `renalprog/modeling/train.py`
- [x] Migrate VAE_plus_bias (postprocessing network) class
- [x] Migrate loss functions from `src/models/train_model.py`
- [x] Migrate training loop from `src/models/train_model.py`
- [x] Implement `train_vae()` function with full functionality
- [x] Add model checkpointing and loading
- [x] Create tests for VAE models (`tests/test_modeling.py`)
- [x] Fix configuration parameters (BETA, CHECKPOINT_FREQ)
- [x] Complete integration with checkpointing system
- [x] Documentation and code quality review

**Output**: 685 lines of production-ready code + 15 comprehensive tests
**Additional Work Completed**:
- [x] Implement NetworkReconstruction (postprocessing network)
- [x] Implement train_vae_with_postprocessing() full pipeline
- [x] Replace matplotlib with Plotly for all visualizations
- [x] Add multi-format export (HTML, PNG, PDF, SVG) for all plots
- [x] Add kaleido dependency for static image export

**Status**: ✅ **ALL 15/15 TESTS PASSING + POSTPROCESSING + PLOTLY** - Day 2 fully complete!

### Day 3 (Nov 29) - Trajectory Generation ✅ COMPLETE (Dec 16)
**Primary Goal**: Implement trajectory generation and patient connections

- [x] Migrate interpolation functions from `src/data/fun_interpol.py`
- [x] Migrate trajectory generation from `src_deseq_and_gsea_NCSR/synthetic_data_generation.py`
- [x] Implement patient connection logic from `notebooks/4_1_trajectories.ipynb`
- [x] Implement control trajectory generation
- [x] Create tests for trajectory functions
- [x] Test trajectory generation on small dataset

**Status**: ✅ **TRAJECTORY GENERATION COMPLETE** - Implemented in predict.py

### Day 4 (Nov 30) - Prediction and Trajectories
- [ ] Implement `renalprog/modeling/predict.py` with trajectory generation
- [ ] Migrate interpolation functions
- [ ] Implement patient connection logic
- [ ] Implement control trajectory generation

### Day 5 (Dec 1) - Classification
- [ ] Migrate XGBoost classifier to modeling module
- [ ] Implement training on real data (all genes + important genes)
- [ ] Implement trajectory classification
- [ ] Test classifier on small dataset

### Day 6 (Dec 2) - Enrichment Analysis ✅ COMPLETE (Dec 17)
**Primary Goal**: Implement dynamic enrichment analysis pipeline

- [x] Create `scripts/pipeline_steps/6_enrichment_analysis.py` main orchestration script
- [x] Migrate DESeq analysis (`py_deseq.py`) to generate fold-change data
- [x] Implement GSEA command generation and parallel execution
- [x] Migrate trajectory analysis (`trajectory_analysis.py`) to combine GSEA results
- [x] Add multiprocessing support with configurable thread count
- [x] Create visualization functions for enrichment results
- [x] Document GSEA installation and GMT file usage
- [x] Test enrichment pipeline on small dataset

**Output**: 
- 750+ lines in `renalprog/modeling/enrichment.py` (EnrichmentPipeline class)
- 210 lines in `scripts/pipeline_steps/6_enrichment_analysis.py` (orchestration)
- 400+ lines in `tests/test_enrichment.py` (comprehensive test suite)
- Complete documentation in `docs/docs/ENRICHMENT_ANALYSIS.md`
- GSEA installation guide in `docs/docs/GSEA_INSTALLATION.md`

**Key Features**:
- Parallel DESeq fold-change calculation
- Parallel GSEA execution with configurable thread count
- Automatic GSEA command generation
- Results combination across all trajectories
- Missing pathway handling
- Progress bars and logging
- Checkpoint/resume capability
- Cleanup utilities

**External Dependencies**:
- GSEA CLI tool (to be downloaded from https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
- ReactomePathways.gmt (already in data/external/ with citation)

**Status**: ✅ **ENRICHMENT PIPELINE COMPLETE** - Day 6 fully complete!

### Day 7 (Dec 3) - Classification ✅ COMPLETE (Dec 16)
**Primary Goal**: Implement static classification and trajectory classification

- [x] Migrate XGBoost classifier to modeling module
- [x] Implement classification_benchmark() with Optuna optimization
- [x] Implement training on real data (all genes + important genes)
- [x] Implement trajectory classification
- [x] Create complete pipeline script (5_classification.py)
- [x] Add visualization functions for metrics and trajectories
- [x] Test classifier integration

**Output**: 350+ lines in train.py + 590 lines in 5_classification.py
**Additional Work Completed**:
- [x] Multi-seed training for robustness
- [x] Best model selection based on Cohen's Kappa
- [x] Publication-ready Plotly visualizations
- [x] Train-to-train and test-to-test trajectory processing
- [x] Comprehensive documentation

**Status**: ✅ **CLASSIFICATION PIPELINE COMPLETE** - Day 7 fully complete!

### Day 8 (Dec 4) - Visualization
- [ ] Migrate plotting functions to `renalprog/plots.py`
- [ ] Create functions for paper figures
- [ ] Update notebooks to use renalprog.plots
- [ ] Test figure generation

### Day 9 (Dec 5) - Scripts and Documentation
- [ ] Organize standalone scripts in `scripts/` folder
- [ ] Create command-line interfaces with argparse
- [ ] Write docstrings for all public functions
- [ ] Start README.md

### Day 10 (Dec 6) - Testing Framework
- [ ] Create test data fixtures (10 samples, 100 genes)
- [ ] Implement `tests/test_dataset.py`
- [ ] Implement `tests/test_features.py`
- [ ] Implement `tests/test_modeling.py`

### Day 11 (Dec 7) - Integration Testing
- [ ] Create end-to-end test with small dataset
- [ ] Test full pipeline execution
- [ ] Fix bugs and issues
- [ ] Verify outputs

### Day 12 (Dec 8) - Performance Testing
- [ ] Test on CPU without GPU
- [ ] Optimize slow operations
- [ ] Verify <5 minute test execution
- [ ] Profile memory usage

### Day 13 (Dec 9) - Documentation
- [ ] Complete README.md with installation and usage
- [ ] Write CONTRIBUTING.md
- [ ] Create example notebooks
- [ ] Document data download process

### Day 14 (Dec 10) - Final Validation
- [ ] Run full pipeline on test data
- [ ] Verify all outputs are generated
- [ ] Check code quality (flake8)
- [ ] Review and cleanup

### Day 15 (Dec 11) - Buffer
- [ ] Address any remaining issues
- [ ] Final testing
- [ ] Prepare for handoff

## Key Decisions Made

1. **R Integration**: Keep R scripts for DESeq2 and clusterProfiler in `scripts/r_analysis/` called via subprocess
2. **GSEA**: Use gseapy Python library for better integration (fallback to bash if needed)
3. **Model Checkpoints**: Use dated folders but add semantic versioning support in config
4. **Test Data**: Create synthetic minimal dataset programmatically for easy redistribution
5. **Package Name**: renalprog (renal progression)

## Notes

- Current repo is labeled BRCA but actually works on KIRC data
- Pipeline has 13 main steps plus manual data download
- Heavy dependencies: PyTorch, sklearn, R, bioinformatics tools
- Both Python and R code in pipeline
- Large models and data - need careful path management
- Parallel processing used extensively - keep multiprocessing support
