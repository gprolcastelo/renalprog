"""Modeling subpackage for training and prediction."""


# Use lazy imports to avoid pulling in heavy dependencies during testing
def __getattr__(name):
    if name == "train_vae":
        from renalprog.modeling.train import train_vae

        return train_vae
    elif name == "VAE":
        from renalprog.modeling.train import VAE

        return VAE
    elif name == "NetworkReconstruction":
        from renalprog.modeling.train import NetworkReconstruction

        return NetworkReconstruction
    elif name == "generate_trajectory_data":
        from renalprog.modeling.predict import generate_trajectory_data

        return generate_trajectory_data
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "train_vae",
    "VAE",
    "NetworkReconstruction",
    "generate_trajectory_data",
]
