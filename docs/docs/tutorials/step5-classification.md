# Step 5: Classification Pipeline

This guide explains how to train classification models on static data and apply them to synthetic trajectory data for disease progression prediction.

## Overview

The classification pipeline performs the following steps:

1. **Load Data**: Load preprocessed data and train/test splits
2. **Train Classifiers**: Train multiple XGBoost classifiers with different seeds
3. **Select Best Model**: Choose best model based on Cohen's Kappa score
4. **Apply to Trajectories**: Classify synthetic trajectory data
5. **Visualize Results**: Create plots showing progression predictions over time

The pipeline trains classifiers on:
- **Static data**: Real patient gene expression (early vs. late stage)

Then applies classifiers to:
- **Train-to-train trajectories**: Synthetic progressions from training patients
- **Test-to-test trajectories**: Synthetic progressions from test patients (unseen)

## Prerequisites

Before running the classification pipeline, ensure you have:

- **Preprocessed data**: From Step 1 data processing
- **Train/test split**: Created in Step 2 (VAE training)
- **Synthetic trajectories**: Generated in Step 4
- **Python environment**: With XGBoost, pandas, plotly installed
- **Sufficient compute**: Multi-core CPU recommended (uses all available cores)

## Usage

### Basic Usage

```bash
python scripts/pipeline_steps/5_classification.py
```

This will:
- Load preprocessed KIRC data
- Train 10 XGBoost classifiers with different seeds
- Select best model based on Cohen's Kappa
- Apply to train-to-train and test-to-test trajectories
- Generate visualization plots

## Configuration Parameters

Edit the script to customize parameters:

### Data Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cancer_type` | `"KIRC"` | Cancer type identifier |
### Classification Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_seeds` | `10` | Number of classifiers to train with different seeds |
| `n_trials` | `100` | Number of Optuna trials for hyperparameter tuning |
| `n_boosting_rounds` | `100` | Number of XGBoost boosting rounds |
| `num_threads` | `os.cpu_count()-1` | Number of CPU threads to use |

### Trajectory Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_timepoints` | `50` | Number of timepoints in trajectories |

## Processing Steps

### Step 1: Load Data

Loads preprocessed data and train/test splits:
- Gene expression data (samples × genes)
- Clinical metadata (stage labels)
- Train/test patient lists
- Optional: Important genes for feature selection

### Step 2: Train Multiple Classifiers

Trains multiple XGBoost classifiers:

1. **Hyperparameter Optimization**:
   - Uses Optuna for Bayesian optimization
   - Optimizes Cohen's Kappa score
   - Searches over learning rate, depth, regularization, etc.

2. **Training**:
   - Trains on training set only
   - Uses stratified cross-validation
   - Records all metrics (accuracy, precision, recall, F1, Kappa)

3. **Multiple Seeds**:
   - Trains with different random seeds
   - Captures model variance
   - Enables robust model selection

### Step 3: Select Best Model

Selects best performing model:
- Ranks models by Cohen's Kappa score
- Saves best model and metadata
- Logs performance metrics

### Step 4: Apply to Trajectories

Classifies synthetic trajectory data:

1. **Train-to-Train**:
   - Trajectories from training patients
   - Shows classifier behavior on training distribution
   - Expected to predict progression well

2. **Test-to-Test**:
   - Trajectories from held-out test patients
   - **True test of generalization**
   - Most important for model evaluation

3. **Time-Course Classification**:
   - Applies classifier to each timepoint
   - Tracks predicted stage probability over time
   - Visualizes progression dynamics

### Step 5: Visualize Results

Generates multiple visualizations:
- Classification metrics boxplots
- Trajectory predictions over time
- Stage probability heatmaps
- Individual trajectory plots

## Output Files

### Classification Models and Metrics

```
models/YYYYMMDD_classification_KIRC/
├── xgboost_model.json                    # Best trained model
├── best_model_metadata.csv               # Model metadata (seed, kappa, etc.)
├── classification_metrics.csv            # All models' metrics
├── classification_summary.csv            # Summary statistics
└── splits_cv/
    └── fold_*.csv                        # Cross-validation splits

reports/figures/YYYYMMDD_classification_KIRC/
└── boxplot_metrics.html/png/pdf/svg     # Metrics visualization
```

### Trajectory Classifications

```
data/processed/YYYYMMDD_trajectory_classifications/
├── train_to_train/
│   └── predictions_train_to_train.csv    # Training trajectory predictions
├── test_to_test/
│   └── predictions_test_to_test.csv      # Test trajectory predictions
└── figures/
    ├── train_to_train/
    │   ├── heatmap_trajectories.html     # Prediction heatmap
    │   ├── individual_trajectories/      # Individual plots
    │   └── summary_statistics.csv
    └── test_to_test/
        ├── heatmap_trajectories.html
        ├── individual_trajectories/
        └── summary_statistics.csv
```




## Interpreting Results

### Good Classification Performance

Indicators of good performance:
- **Cohen's Kappa > 0.60**: Substantial agreement
- **Smooth progression**: Gradual shift from early→late predictions in trajectories

### Poor Classification Performance

Warning signs:
- **Kappa < 0.40**: Weak classification
- **Sharp jumps in predictions**: Unrealistic progression

### Expected Trajectory Behavior

For early→late trajectories:
1. **Initial timepoints**: High probability of early stage
2. **Middle timepoints**: Gradual transition zone
3. **Final timepoints**: High probability of late stage

## Advanced Usage

### Feature Selection

Use only a set of genes to train the classifier and perform trajectory classification:

```python
USE_IMPORTANT_GENES = True

# Specify path to important genes, e.g.:
important_genes_path = EXTERNAL_DATA_DIR / "genes.csv"
```

### Adjust Hyperparameter Search

Modify optimization ranges in the `classification_benchmark` function:

```python
# Example: Limit tree depth
trial.suggest_int("max_depth", 3, 6)  # Instead of 3, 10

# Example: Adjust learning rate
trial.suggest_float("learning_rate", 0.001, 0.1)  # Instead of 0.001, 0.3
```

### Custom Metrics

Track additional metrics in training loop:

```python
# Add to metrics dictionary
metrics = {
    "Accuracy": accuracy,
    "Custom Metric": custom_metric,
    # ... other metrics
}
```

## Next Steps

After classification:

1. **Analyze results**: Review metrics and trajectory predictions
2. **Proceed to Step 6**: Perform pathway enrichment analysis

```bash
python scripts/enrichment/pipeline.sh
```

## Additional Resources

- [XGBoost Documentation](https://xgboost.readthedocs.io/)
- [Cohen's Kappa](https://en.wikipedia.org/wiki/Cohen%27s_kappa)
- [Optuna Tutorial](https://optuna.readthedocs.io/)
- [Feature Importance with SHAP](https://shap.readthedocs.io/)
- [Classification Metrics Guide](https://scikit-learn.org/stable/modules/model_evaluation.html)
