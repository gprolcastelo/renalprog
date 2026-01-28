"""
Model training functionality for renalprog.

Contains VAE model architectures, loss functions, and training infrastructure.
"""

import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, Optional, Dict, Any, List
import logging
from datetime import datetime
from tqdm import tqdm
from sklearn.preprocessing import MinMaxScaler

from renalprog.config import VAEConfig
from renalprog.utils import set_seed, get_device
from renalprog.modeling.checkpointing import ModelCheckpointer, save_model_config

logger = logging.getLogger(__name__)

# Set default dtype
torch.set_default_dtype(torch.float32)


# ============================================================================
# MODEL ARCHITECTURES
# ============================================================================


class VAE(nn.Module):
    """Variational Autoencoder (VAE).

    Standard VAE implementation with encoder-decoder architecture and
    reparameterization trick for sampling from the latent space.

    Args:
        input_dim: Dimension of input data (number of genes)
        mid_dim: Dimension of hidden layer
        features: Dimension of latent space
        output_layer: Output activation function (default: nn.ReLU)
    """

    def __init__(
        self, input_dim: int, mid_dim: int, features: int, output_layer=nn.ReLU
    ):
        super().__init__()
        self.input_dim = input_dim
        self.mid_dim = mid_dim
        self.features = features
        self.output_layer = output_layer

        # Encoder: input -> mid_dim -> (mu, logvar)
        self.encoder = nn.Sequential(
            nn.Linear(in_features=input_dim, out_features=mid_dim),
            nn.ReLU(),
            nn.Linear(in_features=mid_dim, out_features=features * 2),
        )

        # Decoder: latent -> mid_dim -> reconstruction
        self.decoder = nn.Sequential(
            nn.Linear(in_features=features, out_features=mid_dim),
            nn.ReLU(),
            nn.Linear(in_features=mid_dim, out_features=input_dim),
            output_layer(),
        )

    def reparametrize(self, mu: torch.Tensor, log_var: torch.Tensor) -> torch.Tensor:
        """Reparameterization trick: sample from N(mu, var) using N(0,1).

        Args:
            mu: Mean of the latent distribution
            log_var: Log variance of the latent distribution

        Returns:
            Sampled latent vector
        """
        if self.training:
            std = torch.exp(0.5 * log_var)
            eps = torch.randn_like(std)
            return mu + eps * std
        else:
            # During evaluation, return mean directly
            return mu

    def forward(
        self, x: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass through VAE.

        Args:
            x: Input data (batch_size, input_dim)

        Returns:
            Tuple of (reconstruction, mu, log_var, z)
        """
        # Encode
        encoded = self.encoder(x)
        mu_logvar = encoded.view(-1, 2, self.features)
        mu = mu_logvar[:, 0, :]
        log_var = mu_logvar[:, 1, :]

        # Sample from latent space
        z = self.reparametrize(mu, log_var)

        # Decode
        reconstruction = self.decoder(z)

        return reconstruction, mu, log_var, z


class AE(nn.Module):
    """Standard Autoencoder (without variational inference).

    Similar architecture to VAE but without reparameterization trick.

    Args:
        input_dim: Dimension of input data
        mid_dim: Dimension of hidden layer
        features: Dimension of latent space
        output_layer: Output activation function
    """

    def __init__(
        self, input_dim: int, mid_dim: int, features: int, output_layer=nn.ReLU
    ):
        super().__init__()
        self.input_dim = input_dim
        self.mid_dim = mid_dim
        self.features = features
        self.output_layer = output_layer

        self.encoder = nn.Sequential(
            nn.Linear(in_features=input_dim, out_features=mid_dim),
            nn.ReLU(),
            nn.Linear(in_features=mid_dim, out_features=features),
        )

        self.decoder = nn.Sequential(
            nn.Linear(in_features=features, out_features=mid_dim),
            nn.ReLU(),
            nn.Linear(in_features=mid_dim, out_features=input_dim),
            output_layer(),
        )

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, None, None, torch.Tensor]:
        """Forward pass through AE.

        Args:
            x: Input data

        Returns:
            Tuple of (reconstruction, None, None, z)
            None values for mu and logvar to maintain consistency with VAE
        """
        z = self.encoder(x)
        reconstruction = self.decoder(z)
        return reconstruction, None, None, z


class CVAE(VAE):
    """Conditional Variational Autoencoder.

    VAE that conditions on additional information (e.g., clinical data).

    Args:
        input_dim: Dimension of input data
        mid_dim: Dimension of hidden layer
        features: Dimension of latent space
        num_classes: Number of condition classes
        output_layer: Output activation function
    """

    def __init__(
        self,
        input_dim: int,
        mid_dim: int,
        features: int,
        num_classes: int,
        output_layer=nn.ReLU,
    ):
        super().__init__(input_dim, mid_dim, features, output_layer)
        self.num_classes = num_classes

        # Modified encoder: accepts input + condition
        self.encoder = nn.Sequential(
            nn.Linear(in_features=input_dim + num_classes, out_features=mid_dim),
            nn.ReLU(),
            nn.Linear(in_features=mid_dim, out_features=features * 2),
        )

        # Modified decoder: accepts latent + condition
        self.decoder = nn.Sequential(
            nn.Linear(in_features=features + num_classes, out_features=mid_dim),
            nn.ReLU(),
            nn.Linear(in_features=mid_dim, out_features=input_dim),
        )

    def forward(
        self, x: torch.Tensor, condition: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass through CVAE.

        Args:
            x: Input data
            condition: Conditioning information (one-hot encoded)

        Returns:
            Tuple of (reconstruction, mu, log_var, z)
        """
        # Concatenate input with condition
        x_cond = torch.cat([x, condition], dim=1)

        # Encode
        encoded = self.encoder(x_cond)
        mu_logvar = encoded.view(-1, 2, self.features)
        mu = mu_logvar[:, 0, :]
        log_var = mu_logvar[:, 1, :]

        # Sample
        z = self.reparametrize(mu, log_var)

        # Concatenate latent with condition
        z_cond = torch.cat([z, condition], dim=1)

        # Decode
        reconstruction = self.decoder(z_cond)
        reconstruction = self.output_layer()(reconstruction)

        return reconstruction, mu, log_var, z


# ============================================================================
# LOSS FUNCTIONS
# ============================================================================


def vae_loss(
    reconstruction: torch.Tensor,
    x: torch.Tensor,
    mu: torch.Tensor,
    log_var: torch.Tensor,
    beta: float = 1.0,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Calculate VAE loss: reconstruction loss + KL divergence.

    Args:
        reconstruction: Reconstructed output
        x: Original input
        mu: Mean of latent distribution
        log_var: Log variance of latent distribution
        beta: Weight for KL divergence term (beta-VAE)

    Returns:
        Tuple of (total_loss, reconstruction_loss, kl_divergence)
    """
    # Reconstruction loss (MSE)
    recon_loss = nn.functional.mse_loss(reconstruction, x, reduction="sum")

    # KL divergence: -0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    kl_div = -0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())

    # Total loss
    total_loss = recon_loss + beta * kl_div

    return total_loss, recon_loss, kl_div


def reconstruction_loss(
    reconstruction: torch.Tensor, x: torch.Tensor, reduction: str = "sum"
) -> torch.Tensor:
    """Calculate reconstruction loss (MSE).

    Args:
        reconstruction: Reconstructed output
        x: Original input
        reduction: Reduction method ('sum' or 'mean')

    Returns:
        Reconstruction loss
    """
    return nn.functional.mse_loss(reconstruction, x, reduction=reduction)


def kl_divergence(mu: torch.Tensor, log_var: torch.Tensor) -> torch.Tensor:
    """Calculate KL divergence between approximate posterior and prior.

    Args:
        mu: Mean of approximate posterior
        log_var: Log variance of approximate posterior

    Returns:
        KL divergence
    """
    return -0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())


# ============================================================================
# BETA ANNEALING SCHEDULE
# ============================================================================


def frange_cycle_linear(
    start: float, stop: float, n_epoch: int, n_cycle: int = 4, ratio: float = 0.5
) -> np.ndarray:
    """
    Generate a linear cyclical schedule for beta hyperparameter.

    This creates a cyclical annealing schedule where beta increases linearly
    from start to stop over a portion of each cycle (controlled by ratio),
    then stays constant at stop for the remainder of the cycle.

    Args:
        start: Initial value of beta (typically 0.0)
        stop: Final/maximum value of beta (typically 1.0)
        n_epoch: Total number of epochs
        n_cycle: Number of cycles (default: 4)
        ratio: Ratio of cycle spent increasing beta (default: 0.5)
               - 0.5 means half cycle increasing, half constant
               - 1.0 means entire cycle increasing

    Returns:
        Array of beta values for each epoch

    Example:
        >>> # 3 cycles over 300 epochs, beta increases from 0 to 1 over first half of each cycle
        >>> beta_schedule = frange_cycle_linear(0.0, 1.0, 300, n_cycle=3, ratio=0.5)
        >>> # Epoch 0-50: beta increases 0.0 -> 1.0
        >>> # Epoch 50-100: beta stays at 1.0
        >>> # Epoch 100-150: beta increases 0.0 -> 1.0
        >>> # Epoch 150-200: beta stays at 1.0
        >>> # Epoch 200-250: beta increases 0.0 -> 1.0
        >>> # Epoch 250-300: beta stays at 1.0
    """
    L = np.ones(n_epoch) * stop  # Initialize all to stop value
    period = n_epoch / n_cycle
    step = (stop - start) / (period * ratio)  # Linear schedule

    for c in range(n_cycle):
        v, i = start, 0
        while v <= stop and (int(i + c * period) < n_epoch):
            L[int(i + c * period)] = v
            v += step
            i += 1

    return L


# ============================================================================
# DATA LOADING
# ============================================================================


def create_dataloader(
    X: np.ndarray,
    y: Optional[np.ndarray] = None,
    batch_size: int = 32,
    shuffle: bool = True,
) -> torch.utils.data.DataLoader:
    """Create DataLoader with MinMax normalization.

    Args:
        X: Input data (samples x features)
        y: Optional labels
        batch_size: Batch size
        shuffle: Whether to shuffle data

    Returns:
        DataLoader
    """
    # Normalize with MinMaxScaler
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X)

    # Convert to tensors
    X_tensor = torch.tensor(X_scaled, dtype=torch.float32)

    if y is not None:
        y_tensor = torch.tensor(y, dtype=torch.float32)
        dataset = torch.utils.data.TensorDataset(X_tensor, y_tensor)
    else:
        dataset = torch.utils.data.TensorDataset(X_tensor)

    dataloader = torch.utils.data.DataLoader(
        dataset, batch_size=batch_size, shuffle=shuffle
    )

    return dataloader


# ============================================================================
# TRAINING FUNCTIONS
# ============================================================================


def train_epoch(
    model: nn.Module,
    dataloader: torch.utils.data.DataLoader,
    optimizer: torch.optim.Optimizer,
    device: str,
    config: VAEConfig,
    beta: Optional[float] = None,
) -> Dict[str, float]:
    """Train model for one epoch.

    Args:
        model: VAE model
        dataloader: Training DataLoader
        optimizer: Optimizer
        device: Device to use
        config: Training configuration
        beta: Beta value for this epoch (if None, uses config.BETA)

    Returns:
        Dictionary with loss metrics
    """
    if beta is None:
        beta = config.BETA
    model.train()
    total_loss = 0.0
    total_recon = 0.0
    total_kl = 0.0

    # Add progress bar
    pbar = tqdm(dataloader, desc="Training", leave=False)
    for batch in pbar:
        if len(batch) == 2:
            data, _ = batch
        else:
            data = batch[0]

        data = data.to(device)

        # Forward pass
        optimizer.zero_grad()
        reconstruction, mu, log_var, z = model(data)

        # Calculate loss (use beta parameter instead of config.BETA)
        loss, recon, kl = vae_loss(reconstruction, data, mu, log_var, beta)

        # Backward pass
        loss.backward()
        optimizer.step()

        # Accumulate losses
        total_loss += loss.item()
        total_recon += recon.item()
        total_kl += kl.item()

        # Update progress bar
        pbar.set_postfix(
            {
                "loss": f"{loss.item() / len(data):.4f}",
                "recon": f"{recon.item() / len(data):.4f}",
                "kl": f"{kl.item() / len(data):.4f}",
            }
        )

    # Average losses
    n_samples = len(dataloader.dataset)
    metrics = {
        "loss": total_loss / n_samples,
        "recon_loss": total_recon / n_samples,
        "kl_loss": total_kl / n_samples,
    }

    return metrics


def evaluate_model(
    model: nn.Module,
    dataloader: torch.utils.data.DataLoader,
    device: str,
    config: VAEConfig,
    beta: Optional[float] = None,
) -> Dict[str, float]:
    """Evaluate model on validation/test set.

    Args:
        model: VAE model
        dataloader: Validation DataLoader
        device: Device to use
        config: Training configuration
        beta: Beta value for this epoch (if None, uses config.BETA)

    Returns:
        Dictionary with loss metrics
    """
    if beta is None:
        beta = config.BETA
    model.eval()
    total_loss = 0.0
    total_recon = 0.0
    total_kl = 0.0

    with torch.no_grad():
        # Add progress bar
        pbar = tqdm(dataloader, desc="Validation", leave=False)
        for batch in pbar:
            if len(batch) == 2:
                data, _ = batch
            else:
                data = batch[0]

            data = data.to(device)

            # Forward pass
            reconstruction, mu, log_var, z = model(data)

            # Calculate loss (use beta parameter instead of config.BETA)
            loss, recon, kl = vae_loss(reconstruction, data, mu, log_var, beta)

            # Accumulate losses
            total_loss += loss.item()
            total_recon += recon.item()
            total_kl += kl.item()

            # Update progress bar
            pbar.set_postfix(
                {
                    "loss": f"{loss.item() / len(data):.4f}",
                    "recon": f"{recon.item() / len(data):.4f}",
                    "kl": f"{kl.item() / len(data):.4f}",
                }
            )

    # Average losses
    n_samples = len(dataloader.dataset)
    metrics = {
        "loss": total_loss / n_samples,
        "recon_loss": total_recon / n_samples,
        "kl_loss": total_kl / n_samples,
    }

    return metrics


def train_vae(
    X_train: np.ndarray,
    X_test: np.ndarray,
    y_train: Optional[np.ndarray] = None,
    y_test: Optional[np.ndarray] = None,
    config: Optional[VAEConfig] = None,
    save_dir: Optional[Path] = None,
    resume_from: Optional[Path] = None,
    force_cpu: bool = False,
) -> Tuple[nn.Module, Dict[str, list]]:
    """Train a VAE model with full checkpointing support.

    Args:
        X_train: Training data (samples × features) - numpy array or pandas DataFrame
        X_test: Test data (samples × features) - numpy array or pandas DataFrame
        y_train: Optional training labels for CVAE
        y_test: Optional test labels for CVAE
        config: Training configuration
        save_dir: Directory to save checkpoints
        resume_from: Optional checkpoint path to resume training
        force_cpu: Force CPU usage even if CUDA is available (for compatibility)

    Returns:
        Tuple of (trained_model, training_history)
    """
    # Convert DataFrames to numpy arrays if needed
    if hasattr(X_train, "values"):  # Check if it's a DataFrame
        X_train = X_train.values
    if hasattr(X_test, "values"):  # Check if it's a DataFrame
        X_test = X_test.values
    if y_train is not None and hasattr(y_train, "values"):
        y_train = y_train.values
    if y_test is not None and hasattr(y_test, "values"):
        y_test = y_test.values

    if config is None:
        config = VAEConfig()
        config.INPUT_DIM = X_train.shape[1]

    set_seed(config.SEED)

    # Setup save directory
    if save_dir is None:
        timestamp = datetime.now().strftime("%Y%m%d")
        save_dir = Path(f"models/{timestamp}_VAE_KIRC")
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # Save config
    save_model_config(config, save_dir / "config.json")

    # Setup device
    device = get_device(force_cpu=force_cpu)
    logger.info(f"Using device: {device}")

    # Initialize model
    model = VAE(
        input_dim=config.INPUT_DIM,
        mid_dim=config.MID_DIM,
        features=config.LATENT_DIM,
    ).to(device)

    logger.info(
        f"Model: VAE(input_dim={config.INPUT_DIM}, mid_dim={config.MID_DIM}, latent_dim={config.LATENT_DIM})"
    )
    logger.info(f"Parameters: {sum(p.numel() for p in model.parameters()):,}")

    # Setup optimizer
    optimizer = optim.Adam(model.parameters(), lr=config.LEARNING_RATE)

    # Setup checkpointer
    checkpointer = ModelCheckpointer(
        save_dir=save_dir,
        monitor="val_loss",
        mode="min",
        save_freq=config.CHECKPOINT_FREQ,
        keep_last_n=3,
    )

    # Resume from checkpoint if provided
    start_epoch = 0
    if resume_from is not None:
        checkpoint_info = checkpointer.load_checkpoint(
            resume_from, model, optimizer, device=str(device)
        )
        start_epoch = checkpoint_info["epoch"] + 1
        logger.info(f"Resuming training from epoch {start_epoch}")

    # Create dataloaders
    train_loader = create_dataloader(X_train, y_train, config.BATCH_SIZE, shuffle=True)
    test_loader = create_dataloader(X_test, y_test, config.BATCH_SIZE, shuffle=False)

    # Training history
    history = {
        "train_loss": [],
        "val_loss": [],
        "train_recon_loss": [],
        "train_kl_loss": [],
        "val_recon_loss": [],
        "val_kl_loss": [],
        "beta_schedule": [],  # Track beta values
    }

    # Setup beta annealing schedule
    if config.USE_BETA_ANNEALING:
        beta_schedule = frange_cycle_linear(
            start=config.BETA_START,
            stop=config.BETA,
            n_epoch=config.EPOCHS,
            n_cycle=config.BETA_CYCLES,
            ratio=config.BETA_RATIO,
        )
        logger.info(
            f"Using cyclical beta annealing: "
            f"{config.BETA_START} -> {config.BETA} over {config.BETA_CYCLES} cycles"
        )
    else:
        # Constant beta
        beta_schedule = np.ones(config.EPOCHS) * config.BETA
        logger.info(f"Using constant beta: {config.BETA}")

    # Training loop
    logger.info(f"Starting training for {config.EPOCHS} epochs")

    # Add epoch progress bar
    epoch_pbar = tqdm(range(start_epoch, config.EPOCHS), desc="Epochs", position=0)
    for epoch in epoch_pbar:
        # Get beta for this epoch from schedule
        current_beta = beta_schedule[epoch]

        # Train
        train_metrics = train_epoch(
            model, train_loader, optimizer, device, config, beta=current_beta
        )

        # Validate
        val_metrics = evaluate_model(
            model, test_loader, device, config, beta=current_beta
        )

        # Update history
        history["train_loss"].append(train_metrics["loss"])
        history["val_loss"].append(val_metrics["loss"])
        history["train_recon_loss"].append(train_metrics["recon_loss"])
        history["train_kl_loss"].append(train_metrics["kl_loss"])
        history["val_recon_loss"].append(val_metrics["recon_loss"])
        history["val_kl_loss"].append(val_metrics["kl_loss"])
        history["beta_schedule"].append(float(current_beta))

        # Update epoch progress bar
        epoch_pbar.set_postfix(
            {
                "train_loss": f"{train_metrics['loss']:.4f}",
                "val_loss": f"{val_metrics['loss']:.4f}",
                "beta": f"{current_beta:.3f}",
            }
        )

        # Log progress
        if (epoch + 1) % 10 == 0 or epoch == 0:
            logger.info(
                f"Epoch {epoch + 1}/{config.EPOCHS} - "
                f"train_loss: {train_metrics['loss']:.4f}, "
                f"val_loss: {val_metrics['loss']:.4f}"
            )

        # Combine metrics for checkpointing
        current_metrics = {
            "train_loss": train_metrics["loss"],
            "val_loss": val_metrics["loss"],
            "train_recon": train_metrics["recon_loss"],
            "train_kl": train_metrics["kl_loss"],
            "val_recon": val_metrics["recon_loss"],
            "val_kl": val_metrics["kl_loss"],
        }

        # # Save periodic checkpoint
        # if checkpointer.should_save_checkpoint(epoch):
        #     checkpointer.save_checkpoint(
        #         epoch, model, optimizer, current_metrics, config
        #     )

    # Save final model
    checkpointer.save_checkpoint(
        config.EPOCHS - 1, model, optimizer, current_metrics, config, is_final=True
    )

    logger.info("Training complete!")

    return model, history


# ============================================================================
# POSTPROCESSING NETWORK (VAE_plus_bias alternative)
# ============================================================================


class NetworkReconstruction(nn.Module):
    """
    Deep neural network to adjust VAE reconstruction.

    This network is trained on top of VAE output to improve reconstruction
    quality by learning a mapping from VAE reconstruction to original data.

    Args:
        layer_dims: List of layer dimensions [input_dim, hidden1, hidden2, ..., output_dim]
    """

    def __init__(self, layer_dims: List[int]):
        super().__init__()
        layers = []
        for i in range(len(layer_dims) - 1):
            layers.append(nn.Linear(layer_dims[i], layer_dims[i + 1]))
            # Add ReLU after all layers including the final output layer
            # This ensures all outputs are non-negative (suitable for gene expression)
            layers.append(nn.ReLU())
        self.network = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Forward pass through network."""
        return self.network(x)


def train_reconstruction_network(
    network: nn.Module,
    vae_reconstructions: pd.DataFrame,
    original_data: pd.DataFrame,
    train_indices: List,
    test_indices: List,
    epochs: int = 200,
    lr: float = 1e-4,
    batch_size: int = 32,
    device: str = "cpu",
) -> Tuple[nn.Module, List[float], List[float]]:
    """
    Train reconstruction network to adjust VAE output.

    Args:
        network: NetworkReconstruction model
        vae_reconstructions: DataFrame with VAE reconstructions (samples x genes)
        original_data: DataFrame with original gene expression (samples x genes)
        train_indices: List of training sample indices
        test_indices: List of test sample indices
        epochs: Number of training epochs
        lr: Learning rate
        batch_size: Batch size
        device: Device to use

    Returns:
        Tuple of (trained_network, train_losses, test_losses)
    """
    network = network.to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(network.parameters(), lr=lr)

    # Create dataloaders
    train_dataset = torch.utils.data.TensorDataset(
        torch.tensor(
            vae_reconstructions.loc[train_indices].values, dtype=torch.float32
        ),
        torch.tensor(original_data.loc[train_indices].values, dtype=torch.float32),
    )
    test_dataset = torch.utils.data.TensorDataset(
        torch.tensor(vae_reconstructions.loc[test_indices].values, dtype=torch.float32),
        torch.tensor(original_data.loc[test_indices].values, dtype=torch.float32),
    )

    train_loader = torch.utils.data.DataLoader(
        train_dataset, batch_size=batch_size, shuffle=True
    )
    test_loader = torch.utils.data.DataLoader(
        test_dataset, batch_size=batch_size, shuffle=False
    )

    loss_train = []
    loss_test = []

    logger.info(f"Training reconstruction network for {epochs} epochs")

    # Add epoch progress bar
    epoch_pbar = tqdm(range(epochs), desc="Reconstruction Network Training", position=0)
    for epoch in epoch_pbar:
        # Training
        network.train()
        running_loss = 0.0

        # Add batch progress bar
        train_pbar = tqdm(train_loader, desc="Train", leave=False, position=1)
        for vae_recon, original in train_pbar:
            vae_recon = vae_recon.to(device)
            original = original.to(device)

            optimizer.zero_grad()
            output = network(vae_recon)
            loss = criterion(output, original)
            loss.backward()
            optimizer.step()

            running_loss += loss.item()

            # Update batch progress bar
            train_pbar.set_postfix({"batch_loss": f"{loss.item():.6f}"})

        train_loss = running_loss / len(train_loader)
        loss_train.append(train_loss)

        # Validation
        network.eval()
        running_loss = 0.0
        with torch.no_grad():
            # Add validation batch progress bar
            val_pbar = tqdm(test_loader, desc="Val", leave=False, position=1)
            for vae_recon, original in val_pbar:
                vae_recon = vae_recon.to(device)
                original = original.to(device)

                output = network(vae_recon)
                loss = criterion(output, original)
                running_loss += loss.item()

                # Update validation progress bar
                val_pbar.set_postfix({"batch_loss": f"{loss.item():.6f}"})

        test_loss = running_loss / len(test_loader)
        loss_test.append(test_loss)

        # Update epoch progress bar with current metrics
        epoch_pbar.set_postfix(
            {"train_loss": f"{train_loss:.6f}", "test_loss": f"{test_loss:.6f}"}
        )

        if (epoch + 1) % 20 == 0:
            logger.info(
                f"Epoch {epoch + 1}/{epochs} - Train Loss: {train_loss:.6f}, Test Loss: {test_loss:.6f}"
            )

    logger.info("Reconstruction network training complete")
    return network, loss_train, loss_test


def train_vae_with_postprocessing(
    X_train: np.ndarray,
    X_test: np.ndarray,
    vae_config: Optional[VAEConfig] = None,
    reconstruction_network_dims: Optional[List[int]] = None,
    reconstruction_epochs: int = 200,
    reconstruction_lr: float = 1e-4,
    batch_size_reconstruction: int = 8,
    save_dir: Optional[Path] = None,
    force_cpu: bool = False,
) -> Tuple[nn.Module, nn.Module, Dict[str, list], Dict[str, list]]:
    """
    Train VAE followed by postprocessing network (full pipeline).

    This implements the complete training pipeline as in train_vae.sh:
    1. Train VAE on gene expression data
    2. Get VAE reconstructions
    3. Train NetworkReconstruction to adjust VAE output

    Args:
        X_train: Training data (numpy array or pandas DataFrame)
        X_test: Test data (numpy array or pandas DataFrame)
        vae_config: VAE configuration
        reconstruction_network_dims: Architecture for reconstruction network
            If None, defaults to [input_dim, 4096, 1024, 4096, input_dim]
        reconstruction_epochs: Epochs for training reconstruction network
        reconstruction_lr: Learning rate for reconstruction network
        save_dir: Directory to save models
        force_cpu: Force CPU usage

    Returns:
        Tuple of (vae_model, reconstruction_network, vae_history, reconstruction_history)
    """
    logger.info("Starting full VAE + postprocessing pipeline")

    # Convert DataFrames to numpy arrays if needed
    if hasattr(X_train, "values"):  # Check if it's a DataFrame
        X_train = X_train.values
    if hasattr(X_test, "values"):  # Check if it's a DataFrame
        X_test = X_test.values

    # Setup
    if save_dir is None:
        timestamp = datetime.now().strftime("%Y%m%d")
        save_dir = Path(f"models/{timestamp}_VAE_with_reconstruction")
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Train VAE
    logger.info("Step 1: Training VAE")
    vae_model, vae_history = train_vae(
        X_train,
        X_test,
        config=vae_config,
        save_dir=save_dir / "vae",
        force_cpu=force_cpu,
    )

    # Step 2: Get VAE reconstructions
    logger.info("Step 2: Getting VAE reconstructions")
    device = get_device(force_cpu=force_cpu)
    vae_model.eval()

    # Normalize data before passing to VAE (same as during training)
    # The VAE was trained on normalized [0,1] data, so inference must use the same scale
    logger.info("Normalizing data for VAE inference (same as training)")
    scaler = MinMaxScaler()
    X_train_normalized = scaler.fit_transform(X_train)
    X_test_normalized = scaler.transform(X_test)
    logger.info(
        f"Data normalized: min={X_train_normalized.min():.4f}, max={X_train_normalized.max():.4f}"
    )

    with torch.no_grad():
        X_train_tensor = torch.tensor(X_train_normalized, dtype=torch.float32).to(device)
        X_test_tensor = torch.tensor(X_test_normalized, dtype=torch.float32).to(device)

        train_recon_normalized, _, _, _ = vae_model(X_train_tensor)
        test_recon_normalized, _, _, _ = vae_model(X_test_tensor)

        # Denormalize VAE output to match original data scale
        train_recon = scaler.inverse_transform(train_recon_normalized.cpu().numpy())
        test_recon = scaler.inverse_transform(test_recon_normalized.cpu().numpy())

    # Convert to DataFrames
    train_indices = [f"train_{i}" for i in range(len(X_train))]
    test_indices = [f"test_{i}" for i in range(len(X_test))]

    all_recon = np.vstack([train_recon, test_recon])
    all_original = np.vstack([X_train, X_test])
    all_indices = train_indices + test_indices

    df_reconstruction = pd.DataFrame(all_recon, index=all_indices)
    df_original = pd.DataFrame(all_original, index=all_indices)

    # Step 3: Train reconstruction network
    logger.info("Step 3: Training reconstruction network")
    input_dim = X_train.shape[1]

    if reconstruction_network_dims is None:
        reconstruction_network_dims = [input_dim, 4096, 1024, 4096, input_dim]

    network = NetworkReconstruction(reconstruction_network_dims)

    network, loss_train, loss_test = train_reconstruction_network(
        network=network,
        vae_reconstructions=df_reconstruction,
        original_data=df_original,
        train_indices=train_indices,
        test_indices=test_indices,
        epochs=reconstruction_epochs,
        lr=reconstruction_lr,
        batch_size=batch_size_reconstruction,
        device=str(device),
    )

    # Step 4: Save everything
    logger.info("Step 4: Saving models and results")

    # Save reconstruction network
    torch.save(network.state_dict(), save_dir / "reconstruction_network.pth")

    # Save network dimensions
    pd.DataFrame(
        [reconstruction_network_dims],
        columns=["in_dim", "layer1_dim", "layer2_dim", "layer3_dim", "out_dim"],
    ).to_csv(save_dir / "network_dims.csv", index=False)

    # Save losses
    pd.DataFrame({"train_loss": loss_train, "test_loss": loss_test}).to_csv(
        save_dir / "reconstruction_losses.csv", index=False
    )

    # Plot losses using Plotly
    from renalprog.plots import plot_reconstruction_losses

    plot_reconstruction_losses(
        loss_train, loss_test, save_path=save_dir / "reconstruction_losses"
    )

    reconstruction_history = {"train_loss": loss_train, "test_loss": loss_test}

    logger.info(f"Full pipeline complete! Models saved to {save_dir}")

    return vae_model, network, vae_history, reconstruction_history


# Placeholder for additional training functions (TODO)
def train_postprocessing_network(
    vae_model: nn.Module,
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    config: Optional[Dict] = None,
    save_dir: Optional[Path] = None,
) -> nn.Module:
    """
    Train postprocessing network to adjust VAE reconstruction.

    Args:
        vae_model: Trained VAE model
        X_train: Training data
        X_test: Test data
        config: Configuration dictionary
        save_dir: Directory to save model

    Returns:
        Trained postprocessing network model
    """
    logger.info("Training postprocessing network")

    # TODO: Migrate from src/adjust_reconstruction.py
    raise NotImplementedError(
        "train_postprocessing_network() needs implementation from src/adjust_reconstruction.py"
    )


# ============================================================================
# CLASSIFICATION FUNCTIONS
# ============================================================================


def classification_metrics(
    y_test_le: np.ndarray,
    y_pred: np.ndarray,
    num_classes: int,
    weights_test: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Calculate classification metrics for model evaluation.

    Computes precision, recall, F1-score, AUC-ROC, AUC-PR, balanced accuracy,
    and Cohen's Kappa for multi-class classification.

    Args:
        y_test_le: Actual class labels (label-encoded)
        y_pred: Predicted probabilities from XGBoost
        num_classes: Number of classes
        weights_test: Optional sample weights for test set

    Returns:
        DataFrame with classification metrics, one per column
    """
    from sklearn.metrics import (
        classification_report,
        roc_auc_score,
        average_precision_score,
        balanced_accuracy_score,
        cohen_kappa_score,
    )
    from sklearn.preprocessing import label_binarize

    # Initialize y_pred_prob for binary classification
    y_pred_prob = None

    # Binarize test labels for metrics that need it
    if num_classes > 2:
        y_test_bin = label_binarize(y_test_le, classes=np.arange(num_classes))
        y_pred_labels = np.argmax(y_pred, axis=1)
    elif num_classes == 2:
        y_test_bin = label_binarize(y_test_le, classes=[0, 1])
        y_pred_prob = y_pred[:, 1]  # Probability for positive class
        y_pred_labels = np.where(y_pred_prob > 0.5, 1, 0)
    else:
        raise ValueError("Number of classes must be at least 2")

    # Get classification report
    report = classification_report(
        y_true=y_test_le, y_pred=y_pred_labels, output_dict=True, zero_division=0
    )

    # Get macro metrics
    macro_precision = report["macro avg"]["precision"]
    macro_recall = report["macro avg"]["recall"]
    macro_f1_score = report["macro avg"]["f1-score"]

    # Calculate AUC-ROC and AUC-PR
    if num_classes == 2:
        auc_roc = roc_auc_score(y_test_bin, y_pred_prob, sample_weight=weights_test)
        auc_pr = average_precision_score(
            y_test_bin, y_pred_prob, sample_weight=weights_test
        )
    else:
        # For multi-class: calculate macro average
        auc_roc = roc_auc_score(
            y_test_bin,
            y_pred,
            multi_class="ovr",
            average="macro",
            sample_weight=weights_test,
        )
        auc_pr = average_precision_score(
            y_test_bin, y_pred, average="macro", sample_weight=weights_test
        )

    # Calculate balanced accuracy
    balanced_acc = balanced_accuracy_score(
        y_test_le, y_pred_labels, sample_weight=weights_test
    )

    # Calculate Cohen's Kappa
    cohens_kappa = cohen_kappa_score(
        y_test_le, y_pred_labels, sample_weight=weights_test
    )

    # Compile metrics into DataFrame
    # Note: We repeat metrics for each class to match original implementation
    # (this allows setting index to class names, even though metrics are aggregated)
    metrics_data = {
        "Precision": [macro_precision] * num_classes,
        "Recall": [macro_recall] * num_classes,
        "F1-score": [macro_f1_score] * num_classes,
        "AUC-ROC": [auc_roc] * num_classes,
        "AUC-PR": [auc_pr] * num_classes,
        "Balanced Accuracy": [balanced_acc] * num_classes,
        "Cohen's Kappa": [cohens_kappa] * num_classes,
    }

    metrics = pd.DataFrame(metrics_data)

    return metrics


def optimize_xgboost(
    trial,
    dtrain,
    dvalid,
    y_valid,
    num_classes: int,
    num_threads: int,
    progress_bar=None,
) -> float:
    """
    Optuna objective function for XGBoost hyperparameter optimization.

    Args:
        trial: Optuna trial object
        dtrain: Training DMatrix
        dvalid: Validation DMatrix
        y_valid: Validation labels
        num_classes: Number of classes
        num_threads: Number of threads for XGBoost
        progress_bar: Optional tqdm progress bar

    Returns:
        Log loss value (lower is better)
    """
    import xgboost as xgb
    import optuna
    from sklearn.metrics import log_loss

    # Define the parameter space
    param = {
        "verbosity": 0,
        "silent": 1,
        "booster": "gbtree",
        "lambda": trial.suggest_float("lambda", 1e-8, 1, log=True),
        "alpha": trial.suggest_float("alpha", 1e-8, 1, log=True),
    }

    if param["booster"] in ["gbtree", "dart"]:
        param["max_depth"] = trial.suggest_int("max_depth", 1, 9)
        param["eta"] = trial.suggest_float("eta", 1e-8, 1, log=True)
        param["gamma"] = trial.suggest_float("gamma", 1e-8, 1, log=True)
        param["grow_policy"] = trial.suggest_categorical(
            "grow_policy", ["depthwise", "lossguide"]
        )
        param["min_child_weight"] = trial.suggest_int("min_child_weight", 1, 10)

    my_params = {
        "objective": "multi:softprob",
        "tree_method": "exact",
        "num_class": num_classes,
        "nthread": num_threads,
    }

    param = {**param, **my_params}

    # Train the XGBoost model
    model = xgb.train(param, dtrain)

    # Predictions
    preds = model.predict(dvalid)

    # Check for NaN values
    if np.isnan(preds).any():
        logger.warning("Predictions contain NaN values, pruning trial")
        if progress_bar is not None:
            progress_bar.update(1)
        raise optuna.TrialPruned()

    # Calculate log loss
    logloss = log_loss(y_valid, preds)

    # Update progress bar
    if progress_bar is not None:
        progress_bar.update(1)

    return logloss


def classification_benchmark(
    X_data: pd.DataFrame,
    y_data: pd.Series,
    classification_type: str,
    num_classes: int,
    seed: int = 2023,
    test_size: float = 0.2,
    n_br: int = 100,
    num_threads: Optional[int] = None,
    n_trials: int = 100,
    optimization: bool = True,
    tree_method: str = "exact",
) -> Tuple:
    """
    Train and evaluate an XGBoost classifier with hyperparameter optimization.

    This function performs classification using XGBoost with optional Optuna-based
    hyperparameter optimization. It supports both balanced (weighted) and unbalanced
    classification scenarios.

    Args:
        X_data: Feature matrix (samples × genes)
        y_data: Target labels (sample stages)
        classification_type: 'weighted' for balanced classes or 'unbalanced'
        num_classes: Number of target classes
        seed: Random seed for reproducibility
        test_size: Fraction of data for testing (default: 0.2)
        n_br: Number of boosting rounds (default: 100)
        num_threads: Number of threads for XGBoost (default: cpu_count - 1)
        n_trials: Number of Optuna optimization trials (default: 100)
        optimization: Whether to perform hyperparameter optimization
        tree_method: XGBoost tree method ('exact', 'approx', 'hist')

    Returns:
        Tuple containing:
            - model: Trained XGBoost model
            - metrics: DataFrame with classification metrics
            - y_test_le: Label-encoded test labels
            - y_pred: Predicted probabilities
            - data_cv: Tuple of (X_train, X_test, y_train, y_test)
            - all_params: Dictionary of XGBoost parameters used
    """
    import xgboost as xgb
    import optuna
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import LabelEncoder
    from sklearn.utils.class_weight import compute_class_weight

    # Validate inputs
    assert classification_type in ["unbalanced", "weighted"], (
        "classification_type must be 'unbalanced' or 'weighted'"
    )
    assert tree_method in ["exact", "approx", "hist"], (
        "tree_method must be 'exact', 'approx', or 'hist'"
    )

    # Set default num_threads
    if num_threads is None:
        import os

        num_threads = max(1, os.cpu_count() - 1)

    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(
        X_data, y_data, stratify=y_data, test_size=test_size, random_state=seed
    )

    # Encode labels
    le_train = LabelEncoder()
    le_test = LabelEncoder()
    le_train.fit(y_train)
    le_test.fit(y_test)
    y_train_le = le_train.transform(y_train)
    y_test_le = le_test.transform(y_test)

    # Get class dictionary
    classes = np.unique(y_data)
    classes_dict = {key: value for value, key in enumerate(classes)}

    # Handle class weights for balanced classification
    if classification_type == "weighted":
        class_weights = compute_class_weight("balanced", classes=classes, y=y_data)

        # Create weight dictionary
        dict(zip(classes, class_weights))

        # Assign weights to samples
        weights_train = np.zeros(len(y_train_le))
        weights_test = np.zeros(len(y_test_le))

        for i in classes_dict.values():
            weights_train[y_train_le == i] = class_weights[i]
            weights_test[y_test_le == i] = class_weights[i]

        # Create DMatrix with weights
        dtrain_clf = xgb.DMatrix(X_train, label=y_train_le, weight=weights_train)
        dtest_clf = xgb.DMatrix(X_test, label=y_test_le, weight=weights_test)
    else:
        # Create DMatrix without weights
        dtrain_clf = xgb.DMatrix(X_train, y_train_le)
        dtest_clf = xgb.DMatrix(X_test, y_test_le)
        weights_test = None

    # Base parameters
    my_params = {
        "objective": "multi:softprob",
        "booster": "gbtree",
        "tree_method": tree_method,
        "num_class": num_classes,
        "nthread": num_threads,
    }

    # Hyperparameter optimization
    if optimization:
        # Run Optuna optimization silently
        optuna.logging.set_verbosity(optuna.logging.ERROR)
        study = optuna.create_study(direction="minimize")

        # Optimize without progress bar
        study.optimize(
            lambda trial: optimize_xgboost(
                trial,
                dtrain_clf,
                dtest_clf,
                y_test_le,
                num_classes,
                num_threads,
                progress_bar=None,
            ),
            n_trials=n_trials,
            show_progress_bar=False,
        )

        # Get best parameters
        params = study.best_params
        all_params = {**my_params, **params}
    else:
        all_params = my_params

    # Train the model (silently)
    model = xgb.train(all_params, dtrain_clf, num_boost_round=n_br, verbose_eval=False)

    # Make predictions
    y_pred = model.predict(dtest_clf)

    # Calculate metrics
    metrics = classification_metrics(y_test_le, y_pred, num_classes, weights_test)
    metrics.index = classes_dict.keys()

    return (
        model,
        metrics,
        y_test_le,
        y_pred,
        (X_train, X_test, y_train, y_test),
        all_params,
    )


def train_stage_classifier(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    X_test: pd.DataFrame,
    y_test: np.ndarray,
    model_type: str = "xgboost",
    config: Optional[Dict] = None,
    save_dir: Optional[Path] = None,
) -> Any:
    """
    Train a classifier for cancer stage prediction.

    Args:
        X_train: Training features
        y_train: Training labels
        X_test: Test features
        y_test: Test labels
        model_type: Type of classifier ("xgboost", "random_forest", etc.)
        config: Configuration dictionary
        save_dir: Directory to save model

    Returns:
        Trained classifier model
    """
    logger.info(f"Training {model_type} stage classifier")

    if model_type == "xgboost":
        # TODO: Migrate from src/xgb_classifier.py
        raise NotImplementedError(
            "XGBoost classifier training needs migration from src/xgb_classifier.py"
        )
    else:
        raise ValueError(f"Unknown model type: {model_type}")
