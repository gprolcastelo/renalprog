"""
KIRC Classification Pipeline - Step 5
=====================================

This script performs static classification of cancer stages using XGBoost,
then applies the best-performing model to synthetic trajectory data.

Pipeline steps:
1. Load preprocessed data and train/test splits
2. Train multiple XGBoost classifiers with different seeds
3. Select best model based on Cohen's Kappa
4. Apply model to synthetic trajectories (train-to-train and test-to-test)
5. Visualize trajectory classifications

Author: Renalprog Team
Date: 2025-12-16
"""

import os

# Set the maximum number of threads to a higher value
os.environ["NUMEXPR_MAX_THREADS"] = "112"
import sys
import logging
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import xgboost as xgb
from tqdm import tqdm
import plotly.express as px
import plotly.graph_objects as go
from renalprog.config import DATA_DIR, MODELS_DIR, REPORTS_DIR
from renalprog.modeling.train import classification_benchmark
from renalprog.utils import set_seed

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))


# ============================================================================
# CONFIGURATION
# ============================================================================

# Set up logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Paths
INTERIM_DATA_DIR = DATA_DIR / "interim"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
EXTERNAL_DATA_DIR = DATA_DIR / "external"

# Parameters
cancer_type = "KIRC"
stage_col_name = "stage"
today = datetime.now().strftime("%Y%m%d")

# Gene selection
USE_IMPORTANT_GENES = (
    False  # Set to True to use only important genes, False to use all genes
)

# Classification parameters
n_seeds = 10
n_trials = 10
n_boosting_rounds = 100
num_threads = max(1, os.cpu_count() - 1)

# Trajectory parameters
n_timepoints = 50


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def load_classification_data(
    preprocessed_dir,
    train_test_split_dir,
    important_genes_path=None,
    use_important_genes=False,
    stage_col_name="ajcc_pathologic_tumor_stage",
):
    """
    Load data for classification.

    Args:
        preprocessed_dir: Directory containing preprocessed data
        train_test_split_dir: Directory containing train/test split files
        important_genes_path: Optional path to important genes file
        use_important_genes: Whether to use only important genes (True) or all genes (False)

    Returns:
        Tuple of (X_data, y_data, train_patients, test_patients, selected_genes, description)
    """
    logger.info("Loading classification data...")

    # Load preprocessed data
    data_path = preprocessed_dir / "preprocessed_rnaseq.csv"
    metadata_path = preprocessed_dir / "clinical_data.csv"

    data = pd.read_csv(data_path, index_col=0).T  # Transpose to samples × genes
    metadata = pd.read_csv(metadata_path, index_col=0)

    logger.info(f"Data shape: {data.shape}")
    logger.info(f"Metadata shape: {metadata.shape}")

    # Load train/test split
    y_train = pd.read_csv(train_test_split_dir / "y_train.csv", index_col=0)
    y_test = pd.read_csv(train_test_split_dir / "y_test.csv", index_col=0)

    train_patients = y_train.index
    test_patients = y_test.index

    logger.info(f"Train patients: {len(train_patients)}")
    logger.info(f"Test patients: {len(test_patients)}")

    # Select genes based on parameter
    if use_important_genes:
        if important_genes_path and important_genes_path.exists():
            logger.info(f"Using important genes from: {important_genes_path}")
            important_genes = pd.read_csv(important_genes_path, index_col=0)
            selected_genes = [g for g in important_genes.index if g in data.columns]
            logger.info(f"Selected {len(selected_genes)} important genes")
            description = "important_genes"
        else:
            logger.warning(f"Important genes file not found at: {important_genes_path}")
            logger.warning("Falling back to using all genes")
            selected_genes = data.columns.tolist()
            description = "all_genes"
    else:
        logger.info("Using all genes")
        selected_genes = data.columns.tolist()
        description = "all_genes"

    logger.info(f"Total genes selected: {len(selected_genes)}")

    # Prepare data for classification (train set only)
    X_data = data.loc[train_patients, selected_genes]
    y_data = metadata.loc[train_patients, stage_col_name]

    logger.info(f"Classification data shape: X={X_data.shape}, y={y_data.shape}")

    # Create beautified stage distribution table
    stage_counts = y_data.value_counts().sort_index()
    total = len(y_data)

    logger.info("Stage distribution:")
    logger.info("  ┌─────────────┬───────┬────────────┐")
    logger.info("  │ Stage       │ Count │ Percentage │")
    logger.info("  ├─────────────┼───────┼────────────┤")
    for stage, count in stage_counts.items():
        percentage = (count / total) * 100
        logger.info(f"  │ {stage:11s} │ {count:5d} │ {percentage:6.2f}%   │")
    logger.info("  ├─────────────┼───────┼────────────┤")
    logger.info(f"  │ Total       │ {total:5d} │ 100.00%   │")
    logger.info("  └─────────────┴───────┴────────────┘")

    return X_data, y_data, train_patients, test_patients, selected_genes, description


def train_multiple_classifiers(X_data, y_data, n_seeds, output_dir):
    """
    Train multiple XGBoost classifiers with different random seeds.

    Args:
        X_data: Feature matrix
        y_data: Target labels
        n_seeds: Number of models to train
        output_dir: Directory to save results

    Returns:
        Tuple of (models, results_df, seeds)
    """
    logger.info("-" * 80)
    logger.info(f"Training {n_seeds} XGBoost classifiers...")
    logger.info("-" * 80)

    # Generate random seeds
    seeds = np.random.randint(0, 1e9, n_seeds)

    # Storage for results
    models = []
    metrics_list = []
    data_splits = []
    optimized_params = []

    # Temporarily reduce logging level to suppress verbose output
    original_level = logging.getLogger("renalprog.modeling.train").level
    logging.getLogger("renalprog.modeling.train").setLevel(logging.WARNING)

    # Train models with clean progress bar
    pbar = tqdm(
        seeds,
        desc="Training classifiers",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
    )

    for seed_i in pbar:
        # Update progress description
        pbar.set_postfix_str(f"Seed: {seed_i}")

        result = classification_benchmark(
            X_data=X_data,
            y_data=y_data,
            classification_type="weighted",
            num_classes=2,
            seed=int(seed_i),
            test_size=0.2,
            n_br=n_boosting_rounds,
            num_threads=num_threads,
            n_trials=n_trials,
        )

        model, metrics, y_test_le, y_pred, data_cv, params = result

        models.append(model)
        metrics_list.append(metrics)
        data_splits.append(data_cv)
        optimized_params.append(params)

    pbar.close()

    # Restore original logging level
    logging.getLogger("renalprog.modeling.train").setLevel(original_level)

    # Save train-test splits
    for split, seed in zip(data_splits, seeds):
        X_train, X_test, y_train_split, y_test_split = split
        y_train_split.to_csv(output_dir / f"y_train_{seed}.csv")
        y_test_split.to_csv(output_dir / f"y_test_{seed}.csv")

    # Save optimized parameters
    df_params = pd.DataFrame(optimized_params, index=seeds)
    df_params.to_csv(output_dir / "optimized_params.csv")

    # Compile results
    results_df = pd.concat(metrics_list)
    results_df.reset_index(inplace=True, drop=True)

    # Drop duplicate rows (keep one entry per seed)
    even_indices = np.arange(0, len(results_df), 2)
    results_df.drop(index=even_indices, inplace=True)
    results_df.reset_index(inplace=True, drop=True)

    results_df.to_csv(output_dir / "classification_results.csv")

    # Extract column name to avoid backslash in f-string
    cohens_kappa_col = "Cohen's Kappa"
    logger.info(f"Mean Cohen's Kappa: {results_df[cohens_kappa_col].mean():.4f}")

    return models, results_df, seeds


def plot_metrics(results_df, figures_dir):
    """
    Create boxplot of classification metrics.

    Args:
        results_df: DataFrame with classification metrics
        figures_dir: Directory to save figures
    """
    logger.info("Creating metrics visualization...")

    # Melt dataframe for plotting
    df_melted = results_df.melt(var_name="Metric", value_name="Score")

    # Create boxplot
    fig = px.box(
        df_melted,
        x="Metric",
        y="Score",
        color="Metric",
        notched=False,
        color_discrete_sequence=px.colors.qualitative.Safe,
    )

    # Layout customization
    fig.update_layout(
        width=600,
        height=450,
        font=dict(family="Arial, sans-serif", size=12),
        xaxis=dict(
            title="Score",
            tickfont=dict(size=13),
            linecolor="black",
            showgrid=False,
        ),
        yaxis=dict(
            title=None,
            tickfont=dict(size=13),
            range=[0, 1.05],
            linecolor="black",
            gridcolor="lightgrey",
            gridwidth=0.5,
        ),
        showlegend=False,
        plot_bgcolor="white",
        paper_bgcolor="white",
    )

    fig.update_traces(
        width=0.6,
        marker=dict(size=4, opacity=0.5),
    )

    # Save figure
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig.write_html(figures_dir / "boxplot_metrics.html")
    fig.write_image(figures_dir / "boxplot_metrics.png", scale=2)
    fig.write_image(figures_dir / "boxplot_metrics.pdf")
    fig.write_image(figures_dir / "boxplot_metrics.svg")

    logger.info(f"Metrics plots saved to: {figures_dir}")


def save_best_model(models, results_df, seeds, model_dir):
    """
    Save the best performing model.

    Args:
        models: List of trained models
        results_df: Classification metrics dataframe
        seeds: Random seeds used
        model_dir: Directory to save model

    Returns:
        Tuple of (best_model, best_model_idx, best_seed)
    """
    logger.info("Selecting and saving best model...")

    # Find best model based on Cohen's Kappa
    best_model_idx = results_df["Cohen's Kappa"].idxmax()
    best_model = models[best_model_idx]
    best_seed = seeds[best_model_idx]
    best_kappa = results_df.iloc[best_model_idx]["Cohen's Kappa"]

    # Save model
    model_dir.mkdir(parents=True, exist_ok=True)
    model_path = model_dir / "xgboost_model.json"
    best_model.save_model(str(model_path))

    # Save metadata
    metadata = {
        "best_model_idx": int(best_model_idx),
        "best_seed": int(best_seed),
        "cohens_kappa": float(best_kappa),
    }
    pd.Series(metadata).to_csv(model_dir / "best_model_metadata.csv")

    logger.info(f"Best model (index {best_model_idx}, seed {best_seed}) saved")
    logger.info(f"Cohen's Kappa: {best_kappa:.4f}")

    return best_model, best_model_idx, best_seed


def classify_trajectories_from_files(trajectory_dir, model, selected_genes):
    """
    Apply classifier to trajectory files (memory-efficient version).

    Processes trajectories one file at a time to avoid loading all data into memory.

    Args:
        trajectory_dir: Directory containing trajectory CSV files
        model: Trained XGBoost model
        selected_genes: List of genes used for training

    Returns:
        DataFrame with predictions
    """
    csv_files = [f for f in trajectory_dir.glob("*.csv")]

    predictions_all = []

    pbar = tqdm(
        csv_files,
        desc="Classifying trajectories",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
    )

    for traj_id, csv_file in enumerate(pbar):
        # Load one trajectory at a time
        # Trajectory files have timepoints as rows and genes as columns
        df_traj = pd.read_csv(csv_file, index_col=0)

        # Check which selected genes are actually present in the trajectory data
        available_genes = [g for g in selected_genes if g in df_traj.columns]

        # On first trajectory, provide helpful diagnostics if no genes match
        if traj_id == 0:
            if len(available_genes) == 0:
                logger.error("=" * 80)
                logger.error("GENE MISMATCH ERROR")
                logger.error("=" * 80)
                logger.error(
                    "No common genes found between training and trajectory data!"
                )
                logger.error(f"Training data has {len(selected_genes)} genes")
                logger.error(f"Trajectory data has {len(df_traj.columns)} genes")
                logger.error(f"First 10 training genes: {selected_genes[:10]}")
                logger.error(f"First 10 trajectory genes: {list(df_traj.columns[:10])}")
                logger.error(f"Training gene types: {type(selected_genes[0])}")
                logger.error(f"Trajectory gene types: {type(df_traj.columns[0])}")
                logger.error("=" * 80)
                raise ValueError(
                    "No matching genes between trajectory and training data. Check gene name format!"
                )
            elif len(available_genes) < len(selected_genes):
                logger.warning(
                    f"Only {len(available_genes)}/{len(selected_genes)} genes match"
                )

        if len(available_genes) == 0:
            continue

        # Select genes in the correct order
        new_data = df_traj[available_genes]

        # Create DMatrix and predict
        dnew = xgb.DMatrix(new_data)
        predictions_i = model.predict(dnew)

        # Format predictions
        df_pred_i = pd.DataFrame(
            predictions_i, index=new_data.index, columns=["early", "late"]
        )
        df_pred_i = df_pred_i.melt(
            value_vars=["early", "late"], var_name="Stage", value_name="Probability"
        )
        df_pred_i.insert(0, "Traj_ID", traj_id)
        df_pred_i.insert(1, "Interpol_Index", list(range(n_timepoints)) * 2)
        predictions_all.append(df_pred_i)

    pbar.close()

    if len(predictions_all) == 0:
        raise ValueError(
            "No trajectories could be processed. Check gene name matching."
        )

    df_predictions = pd.concat(predictions_all, axis=0)
    logger.info(f"Generated predictions for {len(csv_files)} trajectories")

    return df_predictions


# Keep old functions for backward compatibility but mark as deprecated
def load_trajectory_data(trajectory_dir):
    """
    Load trajectory gene expression data from CSV files.

    DEPRECATED: This function loads all data into memory and can cause OOM errors.
    Use classify_trajectories_from_files() instead for memory efficiency.

    Args:
        trajectory_dir: Directory containing trajectory CSV files

    Returns:
        DataFrame with trajectory data
    """
    dfs = []
    t_id = 0

    csv_files = [f for f in trajectory_dir.glob("*.csv")]

    pbar = tqdm(
        csv_files,
        desc="Loading trajectories",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
    )

    for csv_file in pbar:
        df_gene_i = pd.read_csv(csv_file, index_col=0).T
        df_gene_i.insert(0, "ID", t_id)
        df_gene_i.insert(1, "Trajectory", csv_file.stem)
        dfs.append(df_gene_i)
        t_id += 1

    pbar.close()

    df_gene = pd.concat(dfs, axis=0, ignore_index=False)
    logger.info(f"Loaded {t_id} trajectories")

    return df_gene


def classify_trajectories(traj_gene, model, selected_genes):
    """
    Apply classifier to trajectory data.

    DEPRECATED: Use classify_trajectories_from_files() for better memory efficiency.

    Args:
        traj_gene: Trajectory gene expression data
        model: Trained XGBoost model
        selected_genes: List of genes used for training

    Returns:
        DataFrame with predictions
    """
    all_traj_ids = traj_gene["ID"].unique()
    predictions_all = []

    pbar = tqdm(
        all_traj_ids,
        desc="Classifying trajectories",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
    )

    for traj_id in pbar:
        df_traj = traj_gene[traj_gene["ID"] == traj_id]

        # Prepare data
        new_data = df_traj.drop(["ID", "Trajectory"], axis=1).T
        new_data = new_data[selected_genes]  # Ensure correct gene order

        # Create DMatrix and predict
        dnew = xgb.DMatrix(new_data)
        predictions_i = model.predict(dnew)

        # Format predictions
        df_pred_i = pd.DataFrame(
            predictions_i, index=new_data.index, columns=["early", "late"]
        )
        df_pred_i = df_pred_i.melt(
            value_vars=["early", "late"], var_name="Stage", value_name="Probability"
        )
        df_pred_i.insert(0, "Traj_ID", traj_id)
        df_pred_i.insert(1, "Interpol_Index", list(range(n_timepoints)) * 2)
        predictions_all.append(df_pred_i)

    pbar.close()

    df_predictions = pd.concat(predictions_all, axis=0)
    logger.info(f"Generated predictions for {len(all_traj_ids)} trajectories")

    return df_predictions


def plot_trajectory_classification(df_predictions, save_dir, traj_type):
    """
    Create visualization of trajectory classification results.

    Args:
        df_predictions: Trajectory predictions dataframe
        save_dir: Directory to save figure
        traj_type: Type of trajectories (train_to_train or test_to_test)
    """

    # Create boxplot
    fig = px.box(
        df_predictions,
        x="Interpol_Index",
        y="Probability",
        color="Stage",
        points=False,
        color_discrete_sequence=["rgba(55, 126, 184, 0.85)", "rgba(228, 26, 28, 0.85)"],
        template="ggplot2",
    )

    # Calculate and add medians
    medians = (
        df_predictions.groupby(["Interpol_Index", "Stage"])["Probability"]
        .median()
        .reset_index()
    )
    stage_colors = {"early": "rgba(55, 126, 184, 1)", "late": "rgba(228, 26, 28, 1)"}

    for stage in medians["Stage"].unique():
        stage_medians = medians[medians["Stage"] == stage]
        delta = -0.2 if stage == "early" else 0.2
        fig.add_trace(
            go.Scatter(
                x=stage_medians["Interpol_Index"] + delta,
                y=stage_medians["Probability"],
                mode="markers",
                marker=dict(
                    size=10,
                    color=stage_colors[stage],
                    symbol="diamond",
                    line=dict(color="black", width=1.5),
                ),
                name=f"{stage} (Median)",
                showlegend=False,
            )
        )

    # Layout customization
    fig.update_layout(
        width=1600,
        height=600,
        font=dict(family="Times New Roman", size=24),
        legend=dict(
            title="Stage", font=dict(size=22), bordercolor="black", borderwidth=1
        ),
        xaxis=dict(
            title=dict(text="Interpolation Index", font=dict(size=26)),
            tickfont=dict(size=20),
            showgrid=False,
            linecolor="black",
            mirror=True,
            tickvals=list(range(0, n_timepoints)),
            range=[-0.5, n_timepoints - 0.5],
        ),
        yaxis=dict(
            title=dict(text="Probability", font=dict(size=26)),
            tickfont=dict(size=20),
            showgrid=False,
            linecolor="black",
            mirror=True,
        ),
        margin=dict(l=80, r=40, t=20, b=80),
        title=None,
    )

    # Save figure
    save_dir.mkdir(parents=True, exist_ok=True)
    fig.write_html(save_dir / f"{traj_type}_classification.html")
    fig.write_image(save_dir / f"{traj_type}_classification.png")
    fig.write_image(save_dir / f"{traj_type}_classification.pdf", scale=2)
    fig.write_image(save_dir / f"{traj_type}_classification.svg")

    logger.info(f"Saved {traj_type} classification plots")


# ============================================================================
# MAIN PIPELINE
# ============================================================================

if __name__ == "__main__":
    logger.info("=" * 80)
    logger.info("KIRC CLASSIFICATION PIPELINE - STEP 5")
    logger.info("=" * 80)

    # Set random seed
    set_seed(2023)

    # ========================================================================
    # STEP 1: Load Data
    # ========================================================================
    logger.info("-" * 80)
    logger.info("STEP 1: Loading data")
    logger.info("-" * 80)

    preprocessed_dir = INTERIM_DATA_DIR / "preprocessed_KIRC_data"
    train_test_split_dir = INTERIM_DATA_DIR / "train_test_split" / cancer_type
    important_genes_path = EXTERNAL_DATA_DIR / "important_genes_shap.csv"

    X_data, y_data, train_patients, test_patients, selected_genes, description = (
        load_classification_data(
            preprocessed_dir,
            train_test_split_dir,
            important_genes_path,
            use_important_genes=USE_IMPORTANT_GENES,
            stage_col_name=stage_col_name
        )
    )

    # Setup output directories
    output_dir = (
        INTERIM_DATA_DIR / f"classification_{cancer_type.lower()}_{description}"
    )
    figures_dir = (
        REPORTS_DIR
        / "figures"
        / f"{today}_classification_{cancer_type.lower()}_{description}"
    )
    model_dir = (
        MODELS_DIR / f"classification_{cancer_type.lower()}_{description}"
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    # Save gene lists
    pd.Series(selected_genes).to_csv(output_dir / "used_genes.csv")

    # ========================================================================
    # STEP 2: Train Multiple Classifiers
    # ========================================================================
    models, results_df, seeds = train_multiple_classifiers(
        X_data, y_data, n_seeds, output_dir
    )

    # ========================================================================
    # STEP 3: Plot Metrics
    # ========================================================================
    logger.info("-" * 80)
    logger.info("STEP 3: Visualizing metrics")
    logger.info("-" * 80)

    plot_metrics(results_df, figures_dir)

    # ========================================================================
    # STEP 4: Save Best Model
    # ========================================================================
    logger.info("-" * 80)
    logger.info("STEP 4: Saving best model")
    logger.info("-" * 80)

    best_model, best_model_idx, best_seed = save_best_model(
        models, results_df, seeds, model_dir
    )

    # ========================================================================
    # STEP 5: Apply Model to Trajectories
    # ========================================================================
    logger.info("-" * 80)
    logger.info("STEP 5: Applying model to trajectories")
    logger.info("-" * 80)

    # NOTE: We use classify_trajectories_from_files() which processes one file
    # at a time to avoid loading all trajectory data into memory simultaneously.
    # This prevents OOM errors when processing thousands of trajectories.

    # Base trajectory directory
    synthetic_data_dir = (
        PROCESSED_DATA_DIR
        / f"synthetic_data"
        / cancer_type.lower()
        / "early_to_late"
    )

    # Process train-to-train trajectories
    train_traj_dir = synthetic_data_dir / "train_to_train"
    if train_traj_dir.exists():
        logger.info("Processing train-to-train trajectories")

        # Use memory-efficient version that processes one file at a time
        predictions_train = classify_trajectories_from_files(
            train_traj_dir, best_model, selected_genes
        )

        # Save predictions
        predictions_train.to_csv(output_dir / "predictions_train_to_train.csv")

        # Plot
        plot_trajectory_classification(
            predictions_train, figures_dir / "trajectories", "train_to_train"
        )
    else:
        logger.warning("Train trajectory directory not found")

    # Process test-to-test trajectories
    test_traj_dir = synthetic_data_dir / "test_to_test"
    # test_traj_dir = synthetic_data_dir / "test_to_test" / "early_to_late" # testing already existing trajectories
    if test_traj_dir.exists():
        logger.info("Processing test-to-test trajectories")

        # Use memory-efficient version that processes one file at a time
        predictions_test = classify_trajectories_from_files(
            test_traj_dir, best_model, selected_genes
        )

        # Save predictions
        predictions_test.to_csv(output_dir / "predictions_test_to_test.csv")

        # Plot
        plot_trajectory_classification(
            predictions_test, figures_dir / "trajectories", "test_to_test"
        )
    else:
        logger.warning(f"Test trajectory directory not found: {test_traj_dir}")

    # ========================================================================
    # Pipeline Complete
    # ========================================================================
    logger.info("=" * 80)
    logger.info("CLASSIFICATION PIPELINE COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Results saved to: {output_dir}")
    logger.info(f"Figures saved to: {figures_dir}")
    logger.info(f"Model saved to: {model_dir}")
    logger.info("=" * 80)
