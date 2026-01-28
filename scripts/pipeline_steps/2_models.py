from renalprog import dataset
from renalprog.config import INTERIM_DATA_DIR, VAEConfig, get_dated_dir, MODELS_DIR
from renalprog.utils import configure_logging
from renalprog.plots import plot_training_history, plot_reconstruction_losses
from renalprog.modeling.train import train_vae_with_postprocessing
import logging

# Configure logging for the pipeline
configure_logging(level=logging.INFO)

cancer_type = "KIRC"

# Pre-processing data
# preprocessed_data_path = INTERIM_DATA_DIR / f"preprocessed_{cancer_type}_data"
# path_rnaseq = preprocessed_data_path / "preprocessed_rnaseq.csv"
# path_clinical = preprocessed_data_path / "stages.csv"
# Using these for testing:
path_rnaseq = "data/interim/preprocessed_KIRC/preprocessed_rnaseq.csv"
path_clinical = "data/interim/preprocessed_KIRC/clinical_data.csv"

# Create train/test splits
test_size = 0.2
logging.info(
    f"Splitting into train and test set with {int((1 - test_size) * 100)}/{int(test_size * 100)} ratio"
)
traintest_dir = get_dated_dir(INTERIM_DATA_DIR, "train_test_split")
X_train, X_test, y_train, y_test, _, _ = dataset.create_train_test_split(
    rnaseq_path=path_rnaseq,
    clinical_path=path_clinical,
    test_size=test_size,
    seed=2023,
    output_dir=traintest_dir,
)
logging.info(f"Samples in Train set: {X_train.shape[0]}")
logging.info(f"Samples in Test set: {X_test.shape[0]}")
logging.info(f"Train test split saved to {traintest_dir}")


# Configure VAE
vae_config = VAEConfig()
vae_config.INPUT_DIM = X_train.shape[1]  # Number of features (genes)
vae_config.MID_DIM = 512
vae_config.LATENT_DIM = 256
vae_config.BETA_CYCLES = 3
vae_config.EPOCHS = (
    200 * vae_config.BETA_CYCLES
)  # Total epochs = epochs per cycle * number of cycles
vae_config.BETA_RATIO = 0.5
vae_config.BATCH_SIZE = 8

# Configure postprocessing Reconstruction Network
recnet_dims = [X_train.shape[1], 3_512, 824, 3_731, X_train.shape[1]]
recnet_lr = 1e-4
recnet_epochs = 1_000
batch_recnet = 8
use_cpu = True


model_dir = get_dated_dir(MODELS_DIR, "models_KIRC")

# Train VAE followed by postprocessing Reconstruction Network
vae_model, network, vae_history, reconstruction_history = train_vae_with_postprocessing(
    X_train=X_train,
    X_test=X_test,
    vae_config=vae_config,
    reconstruction_network_dims=recnet_dims,
    reconstruction_epochs=recnet_epochs,
    reconstruction_lr=recnet_lr,
    batch_size_reconstruction=batch_recnet,
    save_dir=model_dir,
    force_cpu=use_cpu,
)

logging.info(f"Trained VAE and Reconstruction Network saved to {model_dir}")

# Plot training histories

logging.info("Plotting VAE training history...")
vae_plot = plot_training_history(
    history=vae_history,
    save_path=model_dir / "vae",
    title="VAE Training History",
    log_scale=False,
)
logging.info(
    f"VAE training history plot saved to {model_dir / 'vae' / 'vae_training_history'}"
)

logging.info("Plotting Reconstruction Network training history...")
recnet_plot = plot_reconstruction_losses(
    loss_train=reconstruction_history["train_loss"],
    loss_test=reconstruction_history["test_loss"],
    save_path=model_dir / "reconstruction_network_history",
    title="Reconstruction Network Training History",
)
logging.info(
    f"Reconstruction network history plot saved to {model_dir / 'reconstruction_network_history'}"
)

logging.info("Training and plotting complete!")
