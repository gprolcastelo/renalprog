#!/usr/bin/env python
"""
Quick Start Example Script for renalprog

This script demonstrates the complete basic workflow:
1. Download preprocessed data from Zenodo
2. Create train/test splits
3. Train VAE with postprocessing network
4. Visualize training progress

Based on: docs/docs/tutorials/quickstart.md
"""

from renalprog import dataset
from renalprog.config import VAEConfig, INTERIM_DATA_DIR, MODELS_DIR
from renalprog.modeling.train import train_vae_with_postprocessing
from renalprog.plots import plot_training_history, plot_reconstruction_losses
from renalprog.utils import configure_logging
from pathlib import Path
import logging

# Configure logging
configure_logging(level=logging.INFO)

# ============================================================================
# STEP 1: Download Preprocessed Data from Zenodo
# ============================================================================
logging.info("=" * 80)
logging.info("STEP 1: Downloading Preprocessed Data from Zenodo")
logging.info("=" * 80)

# Download preprocessed KIRC data from Zenodo
rnaseq, clinical = dataset.download_preprocessed_from_zenodo(
    rnaseq_url='https://zenodo.org/records/17987300/files/Static_KIRC.csv?download=1',
    clinical_url='https://zenodo.org/records/17987300/files/nodes_metadata.csv?download=1',
    output_dir=Path('data/interim/preprocessed_KIRC_data')
)

logging.info(f"RNA-seq data shape: {rnaseq.shape}")
logging.info(f"Genes: {rnaseq.shape[0]:,}, Samples: {rnaseq.shape[1]:,}")
logging.info(f"Clinical data shape: {clinical.shape}")

# ============================================================================
# STEP 2: Create Train/Test Splits
# ============================================================================
logging.info("\n" + "=" * 80)
logging.info("STEP 2: Creating Train/Test Splits")
logging.info("=" * 80)

# Create train/test splits (80/20)
traintest_dir = INTERIM_DATA_DIR / "train_test_split_quickstart"
X_train, X_test, y_train, y_test, _, _ = dataset.create_train_test_split(
    rnaseq_path=Path('data/interim/preprocessed_KIRC_data/preprocessed_rnaseq.csv'),
    clinical_path=Path('data/interim/preprocessed_KIRC_data/clinical_data.csv'),
    test_size=0.2,
    seed=2023,
    output_dir=traintest_dir,
    stage_column="stage"
)

logging.info(f"Training samples: {X_train.shape[0]}")
logging.info(f"Test samples: {X_test.shape[0]}")
logging.info(f"Number of genes: {X_train.shape[1]}")

# ============================================================================
# STEP 3: Configure and Train VAE
# ============================================================================
logging.info("\n" + "=" * 80)
logging.info("STEP 3: Configuring and Training VAE")
logging.info("=" * 80)

# Configure VAE
vae_config = VAEConfig()
vae_config.INPUT_DIM = X_train.shape[1]  # Number of genes
vae_config.MID_DIM = 128
vae_config.LATENT_DIM = 16
vae_config.BETA_CYCLES = 1
vae_config.EPOCHS = 100 * vae_config.BETA_CYCLES  # 600 epochs total
vae_config.BETA_RATIO = 0.5
vae_config.BATCH_SIZE = 8

logging.info("VAE Configuration:")
logging.info(f"  Input dim: {vae_config.INPUT_DIM}")
logging.info(f"  Mid dim: {vae_config.MID_DIM}")
logging.info(f"  Latent dim: {vae_config.LATENT_DIM}")
logging.info(f"  Beta cycles: {vae_config.BETA_CYCLES}")
logging.info(f"  Total epochs: {vae_config.EPOCHS}")
logging.info(f"  Batch size: {vae_config.BATCH_SIZE}")

# Configure postprocessing Reconstruction Network
recnet_dims = [X_train.shape[1], 3_512, 824, 3_731, X_train.shape[1]]
logging.info(f"Reconstruction Network dims: {recnet_dims}")

# Create output directory
model_dir = MODELS_DIR / "quickstart_vae"
model_dir.mkdir(parents=True, exist_ok=True)
logging.info(f"Models will be saved to: {model_dir}")

# Train VAE with postprocessing network
logging.info("\nStarting training...")
vae_model, network, vae_history, reconstruction_history = train_vae_with_postprocessing(
    X_train=X_train,
    X_test=X_test,
    vae_config=vae_config,
    reconstruction_network_dims=recnet_dims,
    reconstruction_epochs=10,
    reconstruction_lr=0.1,
    batch_size_reconstruction=8,
    save_dir=model_dir,
    force_cpu=True,  # Set to False to use GPU if available
)

logging.info(f"\nTraining complete! Models saved to {model_dir}")

# ============================================================================
# STEP 4: Visualize Training Progress
# ============================================================================
logging.info("\n" + "=" * 80)
logging.info("STEP 4: Visualizing Training Progress")
logging.info("=" * 80)

# Plot VAE training history
vae_plot_path = model_dir / "vae_training_history.png"
plot_training_history(
    history=vae_history,
    save_path=vae_plot_path,
    title="VAE Training History"
)
logging.info(f"VAE training history saved to: {vae_plot_path}")

# Plot reconstruction network training history
recnet_plot_path = model_dir / "reconstruction_training_history.png"
plot_reconstruction_losses(
    loss_train=reconstruction_history["train_loss"],
    loss_test=reconstruction_history["test_loss"],
    save_path=recnet_plot_path,
    title="Reconstruction Network Training History"
)
logging.info(f"Reconstruction network training history saved to: {recnet_plot_path}")

# ============================================================================
# SUMMARY
# ============================================================================
logging.info("\n" + "=" * 80)
logging.info("QUICKSTART COMPLETE!")
logging.info("=" * 80)
logging.info("\nSummary:")
logging.info(f"  ✓ Downloaded and processed data: {rnaseq.shape[1]:,} samples, {rnaseq.shape[0]:,} genes")
logging.info(f"  ✓ Created train/test splits: {X_train.shape[0]} train, {X_test.shape[0]} test")
logging.info(f"  ✓ Trained VAE: {vae_config.EPOCHS} epochs")
logging.info("  ✓ Trained reconstruction network: 1000 epochs")
logging.info(f"  ✓ Saved models to: {model_dir}")
logging.info("  ✓ Saved training plots")
logging.info("\nNext steps:")
logging.info("  1. Check reconstruction quality: python scripts/pipeline_steps/3_check_reconstruction.py")
logging.info("  2. Generate trajectories: python scripts/pipeline_steps/4_generate_trajectories.py")
logging.info("  3. Run classification: python scripts/pipeline_steps/5_classification.py")
logging.info("  4. Perform enrichment: python scripts/pipeline_steps/6_enrichment.py")
logging.info("=" * 80)
