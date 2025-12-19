"""
Utility functions used across the renalprog package.
"""

import random
import numpy as np
import torch
import logging
import sys
from sklearn.preprocessing import MinMaxScaler

def set_seed(seed: int = 2023):
    """
    Set random seed for reproducibility across numpy, random, and torch.
    
    Args:
        seed: Random seed value (default: 2023)
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def configure_logging(
    level: int = logging.INFO,
    format_string: str = None,
    datefmt: str = "%Y-%m-%d %H:%M:%S",
    handlers: list = None
) -> None:
    """
    Configure logging for the renalprog package with scientific output formatting.

    This function sets up logging with appropriate formatting for both console output
    and optional file logging. It should be called at the start of scripts that use
    the renalprog package to ensure log messages are visible.

    Args:
        level: Logging level (e.g., logging.INFO, logging.DEBUG).
            Default: logging.INFO
        format_string: Custom format string for log messages.
            Default: "[%(levelname)s] %(message)s" for clean scientific output
        datefmt: Date format for timestamps if included in format_string.
            Default: "%Y-%m-%d %H:%M:%S"
        handlers: List of custom logging handlers. If None, configures stdout handler.
            Default: None (uses console output)

    Examples:
        >>> # Basic configuration (recommended for most scripts)
        >>> from renalprog.utils import configure_logging
        >>> configure_logging()

        >>> # Debug mode with timestamps
        >>> configure_logging(
        ...     level=logging.DEBUG,
        ...     format_string="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        ... )

        >>> # Custom handlers (e.g., file logging)
        >>> file_handler = logging.FileHandler("output.log")
        >>> configure_logging(handlers=[file_handler])

    Notes:
        - This function configures the root logger, which affects all package loggers
        - Call this ONCE at the beginning of your script, not in library code
        - For publication-quality output, use the default clean format
    """
    if format_string is None:
        # Clean format for scientific output (no timestamps cluttering the output)
        format_string = "[%(levelname)s] %(message)s"

    if handlers is None:
        # Configure console output to stdout
        handlers = [logging.StreamHandler(sys.stdout)]

    # Configure root logger
    logging.basicConfig(
        level=level,
        format=format_string,
        datefmt=datefmt,
        handlers=handlers,
        force=True  # Override any existing configuration
    )

    # Set level for all renalprog loggers
    for logger_name in ['renalprog', 'renalprog.dataset', 'renalprog.features',
                        'renalprog.modeling', 'renalprog.trajectories', 'renalprog.plots']:
        logger = logging.getLogger(logger_name)
        logger.setLevel(level)


def get_device(force_cpu: bool = False):
    """
    Get the appropriate device for computation (cuda if available, else cpu).
    
    Args:
        force_cpu: If True, force CPU usage even if CUDA is available

    Returns:
        torch.device: Device to use for tensor operations
    """
    if force_cpu:
        return torch.device("cpu")

    if torch.cuda.is_available():
        try:
            # Test CUDA by creating a small tensor
            test_tensor = torch.zeros(1).cuda()
            del test_tensor
            return torch.device("cuda")
        except Exception as e:
            print(f"Warning: CUDA is available but not functional: {e}")
            print("Falling back to CPU")
            return torch.device("cpu")
    else:
        return torch.device("cpu")


def count_parameters(model):
    """
    Count the number of trainable parameters in a PyTorch model.
    
    Args:
        model: PyTorch model
        
    Returns:
        int: Number of trainable parameters
    """
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def save_model_info(model, path, hyperparams=None):
    """
    Save model information including architecture and hyperparameters.
    
    Args:
        model: PyTorch model
        path: Path to save model info
        hyperparams: Dictionary of hyperparameters
    """
    import json
    
    info = {
        "num_parameters": count_parameters(model),
        "architecture": str(model),
    }
    
    if hyperparams:
        info["hyperparameters"] = hyperparams
    
    with open(path, 'w') as f:
        json.dump(info, f, indent=2)


def load_model_info(path):
    """
    Load model information from JSON file.
    
    Args:
        path: Path to model info file
        
    Returns:
        dict: Model information
    """
    import json
    
    with open(path, 'r') as f:
        return json.load(f)

def apply_VAE(data,model_here,y=None):
    """

    :param data:
    :param model_here:
    :param y:
    :return:
    """
    #############################################################################
    # NORMALIZATION:
    # MinMaxScaler
    scaler = MinMaxScaler()
    scaler.fit(data)
    data2 = scaler.transform(data)
    #############################################################################
#     print("torch no grad")
    with torch.no_grad():
        if y is None:
            data_latent, mu, logvar, z = model_here(torch.tensor(data2).float())
            # DE-NORMALIZE DATA:
            data_vae = scaler.inverse_transform(data_latent)
            #return data_vae,mu, logvar, z, None, scaler
        else:
            data_latent, mu, logvar, z = model_here(torch.tensor(data2).float(),torch.tensor(y).float())
            # DE-NORMALIZE DATA:
            data_vae = scaler.inverse_transform(data_latent)
            #return data_vae,mu, logvar, z, condition, scaler
    return data_vae, mu, logvar, z, scaler