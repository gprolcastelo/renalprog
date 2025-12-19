#!/usr/bin/env python3
"""
Quick start example for renalprog package.

This script demonstrates basic usage of the renalprog package for
preprocessing KIRC data and creating train/test splits.

Usage:
    python quickstart_example.py
"""

import logging
from pathlib import Path
import pandas as pd
import numpy as np

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def create_example_data():
    """Create synthetic example data for demonstration."""
    logger.info("Creating synthetic example data...")

    # Create synthetic gene expression data
    np.random.seed(2023)
    n_genes = 100
    n_samples = 20

    genes = [f"GENE_{i:04d}" for i in range(n_genes)]
    samples = [f"TCGA-TEST-{i:02d}" for i in range(n_samples)]

    # Simulate gene expression with some structure
    expression_data = pd.DataFrame(
        np.random.lognormal(mean=5, sigma=2, size=(n_genes, n_samples)),
        index=genes,
        columns=samples
    )

    # Create clinical data
    stages = (
        ["Stage I"] * 5 +
        ["Stage II"] * 5 +
        ["Stage III"] * 5 +
        ["Stage IV"] * 5
    )
    clinical_data = pd.DataFrame({
        "ajcc_pathologic_tumor_stage": stages
    }, index=samples)

    # Save to temporary location
    temp_dir = Path("data/raw/example")
    temp_dir.mkdir(parents=True, exist_ok=True)

    expression_path = temp_dir / "expression.csv"
    clinical_path = temp_dir / "clinical.csv"

    expression_data.to_csv(expression_path)
    clinical_data.to_csv(clinical_path)

    logger.info(f"Created example data:")
    logger.info(f"  Expression: {expression_data.shape}")
    logger.info(f"  Clinical: {clinical_data.shape}")

    return expression_path, clinical_path


def run_preprocessing_example():
    """Run preprocessing pipeline example."""
    logger.info("\n" + "="*60)
    logger.info("STEP 1: Data Preprocessing")
    logger.info("="*60)

    from renalprog import features
    from renalprog.config import get_dated_dir, INTERIM_DATA_DIR

    # Create example data
    expression_path, clinical_path = create_example_data()

    # Load data
    data = pd.read_csv(expression_path, index_col=0)
    logger.info(f"Loaded expression data: {data.shape}")

    # Preprocess
    logger.info("Running preprocessing...")
    processed_data, metadata = features.preprocess_rnaseq(
        data,
        filter_expression=True,
        detect_outliers=True,
        min_threshold=1.0,
        min_sample_fraction=0.1,
        alpha=0.10  # Less stringent for small example
    )

    logger.info(f"Preprocessing complete:")
    logger.info(f"  Original shape: {metadata['original_shape']}")
    logger.info(f"  Final shape: {metadata['final_shape']}")
    logger.info(f"  Steps applied: {', '.join(metadata['steps_applied'])}")

    # Save results
    output_dir = get_dated_dir(INTERIM_DATA_DIR, "example_preprocessed")
    features.save_preprocessing_results(processed_data, metadata, output_dir)
    logger.info(f"Saved results to: {output_dir}")

    return output_dir / "rnaseq_preprocessed.csv", clinical_path


def run_split_example(rnaseq_path, clinical_path):
    """Run train/test split example."""
    logger.info("\n" + "="*60)
    logger.info("STEP 2: Train/Test Split")
    logger.info("="*60)

    from renalprog.dataset import create_train_test_split
    from renalprog.config import get_dated_dir, INTERIM_DATA_DIR

    output_dir = get_dated_dir(INTERIM_DATA_DIR, "example_split")

    logger.info("Creating stratified train/test split...")
    X_train, X_test, y_train, y_test, full_data, full_clinical = create_train_test_split(
        rnaseq_path,
        clinical_path,
        test_size=0.3,  # 30% test for small dataset
        seed=2023,
        use_onehot=True,
        output_dir=output_dir
    )

    logger.info(f"Split complete:")
    logger.info(f"  Train samples: {X_train.shape[0]}")
    logger.info(f"  Test samples: {X_test.shape[0]}")
    logger.info(f"  Features: {X_train.shape[1]}")
    logger.info(f"  Classes: {y_train.shape[1]}")

    # Load and display statistics
    stats = pd.read_csv(output_dir / "split_statistics.csv", index_col=0)
    logger.info(f"\nStage distribution:")
    logger.info(f"\n{stats.to_string()}")

    return output_dir


def run_visualization_example(split_dir):
    """Run visualization example."""
    logger.info("\n" + "="*60)
    logger.info("STEP 3: Visualization")
    logger.info("="*60)

    from renalprog import plots
    from renalprog.config import FIGURES_DIR
    import matplotlib.pyplot as plt

    # Load split data
    X_train = pd.read_csv(split_dir / "X_train.csv", index_col=0)
    y_train = pd.read_csv(split_dir / "y_train.csv", index_col=0)

    logger.info("Creating visualizations...")

    # Create figures directory
    fig_dir = FIGURES_DIR / "example"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Plot heatmap of top genes
    logger.info("  - Gene expression heatmap")
    top_genes = X_train.var(axis=0).nlargest(20).index
    fig = plots.plot_gene_expression_heatmap(
        X_train[top_genes],
        cluster_rows=True,
        cluster_cols=True,
        save_path=fig_dir / "gene_heatmap.png"
    )
    plt.close()

    logger.info(f"Saved visualizations to: {fig_dir}")


def main():
    """Run the complete example pipeline."""
    logger.info("="*60)
    logger.info("RENALPROG QUICK START EXAMPLE")
    logger.info("="*60)
    logger.info("\nThis example demonstrates:")
    logger.info("  1. Data preprocessing (filtering + outlier detection)")
    logger.info("  2. Train/test split creation")
    logger.info("  3. Basic visualization")
    logger.info("\n" + "="*60 + "\n")

    try:
        # Step 1: Preprocessing
        rnaseq_path, clinical_path = run_preprocessing_example()

        # Step 2: Train/test split
        split_dir = run_split_example(rnaseq_path, clinical_path)

        # Step 3: Visualization
        run_visualization_example(split_dir)

        logger.info("\n" + "="*60)
        logger.info("EXAMPLE COMPLETE!")
        logger.info("="*60)
        logger.info("\nNext steps:")
        logger.info("  1. Replace example data with real TCGA-KIRC data")
        logger.info("  2. Train VAE model (when implemented)")
        logger.info("  3. Generate trajectories (when implemented)")
        logger.info("  4. Run enrichment analysis")

    except Exception as e:
        logger.error(f"Error running example: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()

