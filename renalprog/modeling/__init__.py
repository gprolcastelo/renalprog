"""Modeling subpackage for training and prediction."""

# Use lazy imports to avoid pulling in heavy dependencies during testing
def __getattr__(name):
    if name == "train_vae":
        from renalprog.modeling.train import train_vae
        return train_vae
    elif name == "generate_trajectories":
        from renalprog.modeling.predict import generate_trajectories
        return generate_trajectories
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "train_vae",
    "generate_trajectories",
]
