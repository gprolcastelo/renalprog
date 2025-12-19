"""
renalprog: A Python package for simulating kidney cancer progression
with synthetic data generation and machine learning.

This package provides tools for:
- Preprocessing TCGA-KIRC gene expression data
- Training Variational Autoencoders (VAE) for data representation
- Generating synthetic patient trajectories
- Classifying cancer stages
- Performing enrichment analysis on progression pathways
"""

__version__ = "0.1.0"
__author__ = "Guillermo Prol Castelo"
__license__ = "Apache-2.0"

from renalprog import config
from renalprog import dataset
from renalprog import features
from renalprog import plots
from renalprog import enrichment
from renalprog.utils import configure_logging, set_seed, get_device

__all__ = [
    "config",
    "dataset",
    "features",
    "plots",
    "enrichment",
    "configure_logging",
    "set_seed",
    "get_device",
    "__version__",
]
