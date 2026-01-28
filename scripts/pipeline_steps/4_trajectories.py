"""
Patient Trajectory Construction and Synthetic Data Generation Pipeline

EXACT COPY of link_patients.py workflow for KIRC
"""

import pandas as pd
import torch
import json
import logging
from pathlib import Path
from datetime import datetime

from renalprog.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR
from renalprog.utils import configure_logging, get_device, set_seed
from renalprog.modeling.train import VAE, NetworkReconstruction
from renalprog.modeling.predict import (
    calculate_all_possible_transitions,
    link_patients_random,
    build_trajectory_network,
    generate_trajectory_data,
)
from sklearn.preprocessing import MinMaxScaler

# Configure logging
configure_logging(level=logging.INFO)
logger = logging.getLogger(__name__)

# ============================================================================
# Configuration
# ============================================================================
cancer_type = "KIRC"
set_seed(2023)
early_late = True
start_with_first_stage = True
link_next = 5
distance_metric = "wasserstein"
n_timepoints = 50
interpolation_method = "linear"

# Get device
force_cpu = True
device = get_device(force_cpu=force_cpu)

# Setup output directory
today = datetime.now().strftime("%Y%m%d")
output_dir = PROCESSED_DATA_DIR / f"patient_trajectories_{cancer_type}"
output_dir.mkdir(parents=True, exist_ok=True)

logger.info("=" * 80)
logger.info("PATIENT TRAJECTORY CONSTRUCTION PIPELINE")
logger.info("=" * 80)
logger.info(f"Cancer type: {cancer_type}")
logger.info(f"Device: {device}")
logger.info(f"Output directory: {output_dir}")
logger.info("=" * 80)

# ============================================================================
# STEP 1: Load and Preprocess Data
# ============================================================================
logger.info("-" * 80)
logger.info("STEP 1: Loading and preprocessing data")
logger.info("-" * 80)

# Use EXACT same approach as original working implementation
data_path = str(INTERIM_DATA_DIR / "preprocessed_KIRC_data" / "preprocessed_rnaseq.csv")
metadata_path = str(
    INTERIM_DATA_DIR / "preprocessed_KIRC_data" / "clinical_data.csv"
)

# Load the data (EXACT copy from load_data function)
data = pd.read_csv(data_path, index_col=0)
metadata = pd.read_csv(metadata_path, index_col=0)

# Process stage information
if "stage" not in metadata.columns:
    metadata["stage"] = metadata["ajcc_pathologic_tumor_stage"].replace(
        {f"Stage {i}": i for i in ["I", "II", "III", "IV"]}
    )
    if early_late:
        metadata["stage"] = metadata["stage"].replace(
            {"I": "early", "II": "early", "III": "late", "IV": "late"}
        )
    metadata.drop(columns=["ajcc_pathologic_tumor_stage"], inplace=True)

# Process cancer-specific metadata for KIRC
metadata_selection = metadata[["histological_type", "race", "gender", "stage"]].copy()

logger.info(f"Loaded data: {data.shape[0]} genes × {data.shape[1]} patients")

logger.info(
    f"Metadata: {metadata_selection.shape[0]} patients × {metadata_selection.shape[1]} features"
)

# ============================================================================
# STEP 2: Generate Positive Trajectories
# ============================================================================
logger.info("-" * 80)
logger.info("STEP 2: Generating positive trajectories (disease progression)")
logger.info("-" * 80)

all_traj = calculate_all_possible_transitions(
    data=data,
    metadata_selection=metadata_selection,
    distance=distance_metric,
    early_late=early_late,
    negative_trajectory=False,
)

logger.info(f"Generated {len(all_traj)} possible patient transitions")

# Generate random connections (link_next > 1 means random sampling)
all_linked = link_patients_random(
    results_df=all_traj,
    start_with_first_stage=start_with_first_stage,
    link_next=link_next,
    transitions_possible=["early_to_late"],
)

logger.info(f"Linked {len(all_linked)} patient pairs")

# ============================================================================
# STEP 3: Build Complete Trajectories
# ============================================================================
logger.info("-" * 80)
logger.info("STEP 3: Building complete trajectories from patient links")
logger.info("-" * 80)

network, trajectories = build_trajectory_network(all_linked)

# Format trajectories DataFrame
if early_late:
    trajectories_df = pd.DataFrame(trajectories, columns=["early", "late"])
else:
    max_len = max(len(t) for t in trajectories) if trajectories else 2
    columns = ["I", "II", "III", "IV"][:max_len]
    trajectories_df = pd.DataFrame(trajectories, columns=columns)

trajectories_df["TrajectoryID"] = trajectories_df.index

logger.info(f"Built {len(trajectories)} complete trajectories")

# Initialize train/test trajectory variables
train_trajectories = []
test_trajectories = []

# ============================================================================
# STEP 4: Generate Train/Test Splits (Optional)
# ============================================================================
train_test_split_dir = INTERIM_DATA_DIR / "20260112_train_test_split"
X_train_path = train_test_split_dir / "X_train.csv"
X_test_path = train_test_split_dir / "X_test.csv"

# Initialize variables
use_train_test_split = False
train_patients = None
test_patients = None
random_connections_train = None
random_connections_test = None

if X_train_path.exists() and X_test_path.exists():
    logger.info("-" * 80)
    logger.info("STEP 4: Generating train/test splits")
    logger.info("-" * 80)

    use_train_test_split = True

    # Load X_train and X_test - we only need their indices (patient IDs)
    X_train = pd.read_csv(X_train_path, index_col=0)
    X_test = pd.read_csv(X_test_path, index_col=0)

    # Get patient IDs from the indices
    train_patients = X_train.index
    test_patients = X_test.index

    logger.info(f"Train set: {len(train_patients)} patients")
    logger.info(f"Test set: {len(test_patients)} patients")

    # Split trajectories into train and test
    all_traj_train = all_traj[
        (all_traj["source"].isin(train_patients))
        & (all_traj["target"].isin(train_patients))
    ].reset_index(drop=True)

    all_traj_test = all_traj[
        (all_traj["source"].isin(test_patients))
        & (all_traj["target"].isin(test_patients))
    ].reset_index(drop=True)

    logger.info(f"Train transitions: {len(all_traj_train)}")
    logger.info(f"Test transitions: {len(all_traj_test)}")

    # Generate random connections for early→late
    transitions_possible = ["early_to_late"]
    random_connections_train = link_patients_random(
        results_df=all_traj_train,
        start_with_first_stage=True,
        link_next=5,
        transitions_possible=transitions_possible,
    )
    logger.info(
        f"Created {random_connections_train.shape[0]} connections in train data"
    )

    random_connections_test = link_patients_random(
        results_df=all_traj_test,
        start_with_first_stage=True,
        link_next=5,
        transitions_possible=transitions_possible,
    )
    logger.info(f"Created {random_connections_test.shape[0]} connections in test data")

    # Build separate trajectory networks for train and test
    logger.info("Building separate trajectory networks for train and test sets")

    _, train_trajectories = build_trajectory_network(random_connections_train)
    logger.info(f"Train trajectories: {len(train_trajectories)}")

    _, test_trajectories = build_trajectory_network(random_connections_test)
    logger.info(f"Test trajectories: {len(test_trajectories)}")

else:
    logger.info("Train/test split files not found. Skipping split generation.")
    use_train_test_split = False

# ============================================================================
# STEP 5: Create Trajectory Mapping
# ============================================================================
logger.info("-" * 80)
logger.info("STEP 5: Creating trajectory mapping")
logger.info("-" * 80)

# Create mapping of patients to trajectory IDs
mapping_dict = {"patient": [], "TrajectoryID": []}
for traj_id, traj_i in enumerate(trajectories):
    for pat_i in traj_i:
        mapping_dict["patient"].append(pat_i)
        mapping_dict["TrajectoryID"].append(traj_id)

trajectory_mapping_df = pd.DataFrame(mapping_dict)
logger.info(f"Trajectory mapping: {trajectory_mapping_df.shape}")

# ============================================================================
# STEP 6: Save Main Results
# ============================================================================
logger.info("-" * 80)
logger.info("STEP 6: Saving main results")
logger.info("-" * 80)

# Save nodes metadata
metadata_selection.to_csv(output_dir / "nodes_metadata.csv")
logger.info("Saved: nodes_metadata.csv")

# Save all distances
all_traj.to_csv(output_dir / "all_distances.csv", index=False)
logger.info("Saved: all_distances.csv")

# Save selected links
all_linked_filename = f"random_connections_to_{link_next}_next.csv"
all_linked.to_csv(output_dir / all_linked_filename, index=False)
logger.info(f"Saved: {all_linked_filename}")

# Save trajectories
trajectories_df.to_csv(output_dir / "trajectories.csv", index=False)
logger.info("Saved: trajectories.csv")

# Save trajectory mapping
trajectory_mapping_df.to_csv(output_dir / "trajectory_mapping.csv", index=False)
logger.info("Saved: trajectory_mapping.csv")

# Save train/test splits if generated
if use_train_test_split:
    random_connections_train.to_csv(
        output_dir / f"random_connections_to_{link_next}_next_train.csv", index=False
    )
    random_connections_test.to_csv(
        output_dir / f"random_connections_to_{link_next}_next_test.csv", index=False
    )
    logger.info("Saved train/test trajectory splits")

logger.info(f"All main results saved to: {output_dir}")

# ============================================================================
# STEP 7: Generate Negative Control Trajectories
# ============================================================================
logger.info("-" * 80)
logger.info("STEP 7: Generating negative control trajectories")
logger.info("-" * 80)

# Generate all same-stage transitions for KIRC
all_traj_negatives = calculate_all_possible_transitions(
    data=data,
    metadata_selection=metadata_selection,
    distance=distance_metric,
    early_late=early_late,
    negative_trajectory=True,
)

logger.info(f"Generated {len(all_traj_negatives)} control transitions")

# Generate random connections for controls
random_connections_negatives = link_patients_random(
    results_df=all_traj_negatives,
    start_with_first_stage=True,
    link_next=5,
    transitions_possible=["early_to_early", "late_to_late"],
)

# Remove self-connections
random_connections_negatives = random_connections_negatives[
    random_connections_negatives["source"] != random_connections_negatives["target"]
].reset_index(drop=True)

logger.info(
    f"Created {len(random_connections_negatives)} control connections (self-connections removed)"
)

# Save negative trajectories
random_connections_negatives.to_csv(output_dir / "all_negatives.csv", index=False)
logger.info("Saved: all_negatives.csv")

# Split negatives into train/test if applicable
if use_train_test_split:
    all_traj_negatives_train = all_traj_negatives[
        (all_traj_negatives["source"].isin(train_patients))
        & (all_traj_negatives["target"].isin(train_patients))
    ].reset_index(drop=True)

    all_traj_negatives_test = all_traj_negatives[
        (all_traj_negatives["source"].isin(test_patients))
        & (all_traj_negatives["target"].isin(test_patients))
    ].reset_index(drop=True)

    # Generate random connections for controls
    controls_train = link_patients_random(
        results_df=all_traj_negatives_train,
        start_with_first_stage=True,
        link_next=5,
        transitions_possible=["early_to_early", "late_to_late"],
    )
    controls_train = controls_train[
        controls_train["source"] != controls_train["target"]
    ].reset_index(drop=True)

    controls_test = link_patients_random(
        results_df=all_traj_negatives_test,
        start_with_first_stage=True,
        link_next=5,
        transitions_possible=["early_to_early", "late_to_late"],
    )
    controls_test = controls_test[
        controls_test["source"] != controls_test["target"]
    ].reset_index(drop=True)

    # Save control splits
    controls_train.to_csv(output_dir / "controls_train.csv", index=False)
    controls_test.to_csv(output_dir / "controls_test.csv", index=False)
    logger.info(f"Saved: controls_train.csv ({len(controls_train)} connections)")
    logger.info(f"Saved: controls_test.csv ({len(controls_test)} connections)")

# ============================================================================
# STEP 8: Load Trained Models
# ============================================================================
logger.info("-" * 80)
logger.info("STEP 8: Loading trained VAE and Reconstruction Network models")
logger.info("-" * 80)

model_dir = Path("models/20260112_models_KIRC")

# Load VAE
vae_model_path = model_dir / "vae" / "final_model.pth"
vae_config_path = model_dir / "vae" / "config.json"

with open(vae_config_path, "r") as f:
    vae_config_dict = json.load(f)

vae_model = VAE(
    input_dim=vae_config_dict["INPUT_DIM"],
    mid_dim=vae_config_dict["MID_DIM"],
    features=vae_config_dict["LATENT_DIM"],
).to(device)

checkpoint = torch.load(vae_model_path, map_location=device, weights_only=False)
vae_model.load_state_dict(checkpoint["model_state_dict"])
vae_model.eval()
logger.info(f"VAE model loaded from: {vae_model_path}")

# Load Reconstruction Network
recnet_model_path = model_dir / "reconstruction_network.pth"
recnet_dims_path = model_dir / "network_dims.csv"

network_dims = pd.read_csv(recnet_dims_path).values.tolist()[0]
logger.info(f"Reconstruction network dimensions: {network_dims}")

recnet_model = NetworkReconstruction(layer_dims=network_dims).to(device)
recnet_model.load_state_dict(
    torch.load(recnet_model_path, map_location=device, weights_only=False)
)
recnet_model.eval()
logger.info(f"Reconstruction network loaded from: {recnet_model_path}")

# ============================================================================
# Create scaler fitted on TRAIN split only (avoid data leakage)
# ============================================================================
logger.info("-" * 80)
logger.info("Fitting scaler on TRAIN split only (same as VAE training)")
logger.info("-" * 80)

# Load train/test split from the same directory used for VAE training
train_test_split_dir = INTERIM_DATA_DIR / "20251211_train_test_split"
X_train_path = train_test_split_dir / "X_train.csv"

if not X_train_path.exists():
    logger.error(f"Train split not found at {X_train_path}")
    logger.error("Please run 2_models.py first to create train/test split")
    raise FileNotFoundError(f"Train split not found: {X_train_path}")

# Load ONLY the train split (genes × patients)
X_train = pd.read_csv(X_train_path, index_col=0)
logger.info(
    f"Loaded train split: {X_train.shape[0]} patients × {X_train.shape[1]} genes"
)

# Fit scaler ONLY on train data to avoid data leakage
# X_train is (patients × genes) which is what scaler expects
scaler = MinMaxScaler()
scaler.fit(X_train.values)
logger.info(
    f"Scaler fitted on TRAIN split: {X_train.shape[0]} patients × {X_train.shape[1]} genes"
)
logger.info(f"Scaler range: min={scaler.data_min_[:5]}, max={scaler.data_max_[:5]}")

# ============================================================================
# STEP 9: Generate Synthetic Trajectory Data
# ============================================================================
logger.info("-" * 80)
logger.info("STEP 9: Generating synthetic trajectory data")
logger.info("-" * 80)

if use_train_test_split:
    # Train and test trajectories were already built in STEP 4
    logger.info(f"Train trajectories: {len(train_trajectories)}")
    logger.info(f"Test trajectories: {len(test_trajectories)}")

    # Generate train trajectories
    logger.info("-" * 80)
    logger.info("Generating train set trajectories")
    logger.info("-" * 80)

    train_synthetic_dir = (
        PROCESSED_DATA_DIR
        / f"{today}_synthetic_data"
        / cancer_type.lower()
        / "recnet"
        / "early_to_late"
        / "train_to_train"
    )
    train_synthetic_dir.mkdir(parents=True, exist_ok=True)

    for traj_idx, trajectory in enumerate(train_trajectories):
        if (traj_idx + 1) % 50 == 0 or traj_idx == 0:
            logger.info(
                f"Train progress: {traj_idx + 1}/{len(train_trajectories)} trajectories generated"
            )

        trajectory_name = "_to_".join(trajectory)
        save_path = train_synthetic_dir / f"{trajectory_name}.csv"

        try:
            traj_data = generate_trajectory_data(
                vae_model=vae_model,
                recnet_model=recnet_model,
                trajectory=trajectory,
                gene_data=data,
                n_timepoints=n_timepoints,
                interpolation_method=interpolation_method,
                device=str(device),
                save_path=save_path,
                scaler=scaler,
            )
        except Exception as e:
            logger.error(f"Error generating train trajectory {trajectory_name}: {e}")
            continue

    logger.info(f"Train trajectories saved to: {train_synthetic_dir}")

    # Generate test trajectories
    logger.info("-" * 80)
    logger.info("Generating test set trajectories")
    logger.info("-" * 80)

    test_synthetic_dir = (
        PROCESSED_DATA_DIR
        / f"{today}_synthetic_data"
        / cancer_type.lower()
        / "recnet"
        / "early_to_late"
        / "test_to_test"
    )
    test_synthetic_dir.mkdir(parents=True, exist_ok=True)

    for traj_idx, trajectory in enumerate(test_trajectories):
        if (traj_idx + 1) % 50 == 0 or traj_idx == 0:
            logger.info(
                f"Test progress: {traj_idx + 1}/{len(test_trajectories)} trajectories generated"
            )

        trajectory_name = "_to_".join(trajectory)
        save_path = test_synthetic_dir / f"{trajectory_name}.csv"

        try:
            traj_data = generate_trajectory_data(
                vae_model=vae_model,
                recnet_model=recnet_model,
                trajectory=trajectory,
                gene_data=data,
                n_timepoints=n_timepoints,
                interpolation_method=interpolation_method,
                device=str(device),
                save_path=save_path,
                scaler=scaler,
            )
        except Exception as e:
            logger.error(f"Error generating test trajectory {trajectory_name}: {e}")
            continue

    logger.info(f"Test trajectories saved to: {test_synthetic_dir}")

else:
    # Generate all trajectories without train/test split
    logger.info("No train/test split - generating all trajectories together")

    synthetic_output_dir = (
        PROCESSED_DATA_DIR
        / f"{today}_synthetic_data"
        / cancer_type.lower()
        / "recnet"
        / "early_to_late"
    )
    synthetic_output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        f"Generating {len(trajectories)} trajectories with {n_timepoints} timepoints each"
    )

    for traj_idx, trajectory in enumerate(trajectories):
        if (traj_idx + 1) % 50 == 0 or traj_idx == 0:
            logger.info(
                f"Progress: {traj_idx + 1}/{len(trajectories)} trajectories generated"
            )

        trajectory_name = "_to_".join(trajectory)
        save_path = synthetic_output_dir / f"{trajectory_name}.csv"

        try:
            traj_data = generate_trajectory_data(
                vae_model=vae_model,
                recnet_model=recnet_model,
                trajectory=trajectory,
                gene_data=data,
                n_timepoints=n_timepoints,
                interpolation_method=interpolation_method,
                device=str(device),
                save_path=save_path,
                scaler=scaler,
            )
        except Exception as e:
            logger.error(f"Error generating trajectory {trajectory_name}: {e}")
            continue

    logger.info(f"Synthetic trajectory data saved to: {synthetic_output_dir}")

# ============================================================================
# TODO: generate control trajectories (noise in real space)
# ============================================================================

# ============================================================================
# Pipeline Complete
# ============================================================================
logger.info("=" * 80)
logger.info("TRAJECTORY PIPELINE COMPLETE")
logger.info("=" * 80)
logger.info(f"Patient connections and metadata saved to: {output_dir}")
logger.info("")
logger.info("Output files:")
logger.info(f"  - nodes_metadata.csv: {len(metadata_selection)} patients")
logger.info(f"  - all_distances.csv: {len(all_traj)} possible transitions")
logger.info(f"  - random_connections_to_{link_next}_next.csv: {len(all_linked)} links")
logger.info(f"  - trajectories.csv: {len(trajectories)} complete trajectories")
logger.info(
    f"  - trajectory_mapping.csv: {len(trajectory_mapping_df)} patient-trajectory mappings"
)
logger.info(
    f"  - all_negatives.csv: {len(random_connections_negatives)} control connections"
)
if use_train_test_split:
    logger.info("  - Train/test splits saved")
    logger.info(
        f"  - Synthetic data (train): {len(train_trajectories)} trajectory files in train_to_train/"
    )
    logger.info(
        f"  - Synthetic data (test): {len(test_trajectories)} trajectory files in test_to_test/"
    )
else:
    logger.info(f"  - Synthetic data: {len(trajectories)} trajectory files")
logger.info("=" * 80)
