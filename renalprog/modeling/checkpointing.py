"""Model checkpointing utilities for saving and loading training state."""

import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional

import torch
import torch.nn as nn
from torch.optim import Optimizer

logger = logging.getLogger(__name__)


class ModelCheckpointer:
    """Handles saving and loading model checkpoints during training.

    Features:
    - Save best model based on validation metric
    - Save checkpoints every N epochs
    - Save final model after training
    - Save training history and configuration
    - Resume training from checkpoint

    Attributes:
        save_dir: Directory to save checkpoints
        monitor: Metric to monitor ('loss', 'val_loss', etc.)
        mode: 'min' for loss, 'max' for accuracy
        save_freq: Save checkpoint every N epochs (0 = only best)
        keep_last_n: Keep only last N checkpoints (0 = keep all)
    """

    def __init__(
        self,
        save_dir: Path,
        monitor: str = "val_loss",
        mode: str = "min",
        save_freq: int = 0,
        keep_last_n: int = 3,
    ):
        """Initialize checkpointer.

        Args:
            save_dir: Directory to save checkpoints
            monitor: Metric name to monitor
            mode: 'min' to minimize metric, 'max' to maximize
            save_freq: Save every N epochs (0 = only save best)
            keep_last_n: Keep only N most recent checkpoints (0 = all)
        """
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)

        self.monitor = monitor
        self.mode = mode
        self.save_freq = save_freq
        self.keep_last_n = keep_last_n

        # Track best metric
        self.best_metric = float("inf") if mode == "min" else float("-inf")
        self.best_epoch = 0

        # Track saved checkpoints for cleanup
        self.checkpoint_history = []

        logger.info(f"ModelCheckpointer initialized: {save_dir}")
        logger.info(f"Monitoring: {monitor} ({mode})")

    def save_checkpoint(
        self,
        epoch: int,
        model: nn.Module,
        optimizer: Optimizer,
        metrics: Dict[str, float],
        config: Any,
        is_best: bool = False,
        is_final: bool = False,
    ) -> None:
        """Save a training checkpoint.

        Args:
            epoch: Current epoch number
            model: PyTorch model to save
            optimizer: Optimizer state to save
            metrics: Dictionary of current metrics
            config: Training configuration object
            is_best: Whether this is the best model so far
            is_final: Whether this is the final model
        """
        # Create checkpoint dictionary
        checkpoint = {
            "epoch": epoch,
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "metrics": metrics,
            "config": self._config_to_dict(config),
            "best_metric": self.best_metric,
            "monitor": self.monitor,
        }

        # Determine checkpoint filename
        if is_final:
            filename = "final_model.pth"
        elif is_best:
            filename = "best_model.pth"
        else:
            filename = f"checkpoint_epoch_{epoch:04d}.pth"

        checkpoint_path = self.save_dir / filename

        # Save checkpoint
        torch.save(checkpoint, checkpoint_path)
        logger.info(f"Saved checkpoint: {checkpoint_path}")

        # Track for cleanup
        if not is_best and not is_final:
            self.checkpoint_history.append(checkpoint_path)
            self._cleanup_old_checkpoints()

        # Save training history as JSON
        if is_best or is_final:
            self._save_history(metrics, epoch)

    def load_checkpoint(
        self,
        checkpoint_path: Path,
        model: nn.Module,
        optimizer: Optional[Optimizer] = None,
        device: str = "cpu",
    ) -> Dict[str, Any]:
        """Load a checkpoint and restore model state.

        Args:
            checkpoint_path: Path to checkpoint file
            model: Model to load state into
            optimizer: Optional optimizer to restore state
            device: Device to map checkpoint to

        Returns:
            Dictionary with checkpoint information (epoch, metrics, config)
        """
        if not checkpoint_path.exists():
            raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")

        # Load checkpoint
        checkpoint = torch.load(checkpoint_path, map_location=device)

        # Restore model state
        model.load_state_dict(checkpoint["model_state_dict"])
        logger.info(f"Loaded model state from epoch {checkpoint['epoch']}")

        # Restore optimizer state if provided
        if optimizer is not None and "optimizer_state_dict" in checkpoint:
            optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
            logger.info("Loaded optimizer state")

        # Return checkpoint info
        return {
            "epoch": checkpoint["epoch"],
            "metrics": checkpoint.get("metrics", {}),
            "config": checkpoint.get("config", {}),
            "best_metric": checkpoint.get("best_metric", self.best_metric),
        }

    def update_best(self, epoch: int, metric_value: float) -> bool:
        """Check if current metric is the best and update if so.

        Args:
            epoch: Current epoch number
            metric_value: Current metric value

        Returns:
            True if this is a new best, False otherwise
        """
        is_better = (self.mode == "min" and metric_value < self.best_metric) or (
            self.mode == "max" and metric_value > self.best_metric
        )

        if is_better:
            self.best_metric = metric_value
            self.best_epoch = epoch
            logger.info(f"New best {self.monitor}: {metric_value:.6f} at epoch {epoch}")
            return True
        return False

    def should_save_checkpoint(self, epoch: int) -> bool:
        """Determine if checkpoint should be saved this epoch.

        Args:
            epoch: Current epoch number

        Returns:
            True if checkpoint should be saved
        """
        if self.save_freq == 0:
            return False
        return epoch % self.save_freq == 0

    def _config_to_dict(self, config: Any) -> Dict[str, Any]:
        """Convert config object to dictionary.

        Args:
            config: Configuration object (dataclass or class with __dict__)

        Returns:
            Dictionary representation of config
        """
        if hasattr(config, "__dict__"):
            return {k: v for k, v in config.__dict__.items() if not k.startswith("_")}
        return {}

    def _save_history(self, metrics: Dict[str, float], epoch: int) -> None:
        """Save training history to JSON file.

        Args:
            metrics: Current metrics dictionary
            epoch: Current epoch number
        """
        history_path = self.save_dir / "training_history.json"

        # Load existing history or create new
        if history_path.exists():
            with open(history_path, "r") as f:
                history = json.load(f)
        else:
            history = {"epochs": [], "metrics": {}}

        # Append current epoch
        history["epochs"].append(epoch)
        for key, value in metrics.items():
            if key not in history["metrics"]:
                history["metrics"][key] = []
            history["metrics"][key].append(float(value))

        # Save updated history
        with open(history_path, "w") as f:
            json.dump(history, f, indent=2)

    def _cleanup_old_checkpoints(self) -> None:
        """Remove old checkpoints keeping only last N."""
        if self.keep_last_n <= 0:
            return

        while len(self.checkpoint_history) > self.keep_last_n:
            old_checkpoint = self.checkpoint_history.pop(0)
            if old_checkpoint.exists():
                old_checkpoint.unlink()
                logger.debug(f"Removed old checkpoint: {old_checkpoint}")

    def get_best_checkpoint_path(self) -> Optional[Path]:
        """Get path to best model checkpoint.

        Returns:
            Path to best model, or None if not saved yet
        """
        best_path = self.save_dir / "best_model.pth"
        return best_path if best_path.exists() else None

    def get_final_checkpoint_path(self) -> Optional[Path]:
        """Get path to final model checkpoint.

        Returns:
            Path to final model, or None if not saved yet
        """
        final_path = self.save_dir / "final_model.pth"
        return final_path if final_path.exists() else None


def save_model_config(config: Any, save_path: Path) -> None:
    """Save model configuration to JSON file.

    Args:
        config: Configuration object
        save_path: Path to save JSON file
    """
    config_dict = {}
    if hasattr(config, "__dict__"):
        config_dict = {
            k: v
            for k, v in config.__dict__.items()
            if not k.startswith("_") and not callable(v)
        }

    save_path.parent.mkdir(parents=True, exist_ok=True)
    with open(save_path, "w") as f:
        json.dump(config_dict, f, indent=2)

    logger.info(f"Saved config: {save_path}")


def load_model_config(config_path: Path) -> Dict[str, Any]:
    """Load model configuration from JSON file.

    Args:
        config_path: Path to JSON config file

    Returns:
        Dictionary with configuration
    """
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")

    with open(config_path, "r") as f:
        config = json.load(f)

    logger.info(f"Loaded config: {config_path}")
    return config
