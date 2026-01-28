from renalprog import dataset
from renalprog.config import (
    INTERIM_DATA_DIR,
    PROCESSED_DATA_DIR,
    FIGURES_DIR,
    MODELS_DIR,
    VAEConfig,
)
from renalprog.utils import configure_logging, set_seed, get_device, apply_VAE
from renalprog.modeling.train import VAE, NetworkReconstruction
from renalprog.modeling.predict import evaluate_reconstruction
from renalprog.dataset import load_rnaseq_data, load_clinical_data
from renalprog.plots import plot_umap_plotly as plot_umap
import logging
import os
import torch
import pandas as pd
import json
from datetime import datetime
import argparse
parser = argparse.ArgumentParser(description='Check reconstruction of models.')
parser.add_argument('--cancer_type',
                    type=str,
                    default='KIRC',
                    help='Cancer type to process (KIRC or BRCA). Default: KIRC. Note: for downloading models from Hugging Face, only KIRC and BRCA are available.')
parser.add_argument('--stage_column',
                    type=str,
                    default='ajcc_pathologic_tumor_stage',
                    help='Stages column name in clinical data. Default: ajcc_pathologic_tumor_stage. In this repo, it has sometimes been called simply "stage".')
parser.add_argument('--hf_models',
                    action='store_true',
                    help='Load pre-trained models from Hugging Face Hub.')
parser.add_argument('--sdmetrics',
                    action='store_true',
                    help='Use sdmetrics to obtain relevant metrics for the generated test set synthetic data. Takes very long!')
args = parser.parse_args()

# Configure logging for the pipeline
configure_logging(level=logging.INFO)

cancer_type = args.cancer_type
set_seed(2023)
# Get device
device = get_device(force_cpu=True)
logging.info(f"Using device: {device}")

# ============================================================================
# Load preprocessed data
# ============================================================================
logging.info("Loading preprocessed data...")

# Using these paths for consistency with previous pipeline steps
path_rnaseq = f"data/interim/preprocessed_{cancer_type}_data/preprocessed_rnaseq.csv"
path_clinical = f"data/interim/preprocessed_{cancer_type}_data/clinical_data.csv"

preprocessed = load_rnaseq_data(path_rnaseq)
subtypes = load_clinical_data(path_clinical,stage_column=args.stage_column)

# Ensure samples match between data and clinical info
if preprocessed.shape[0]!=subtypes.shape[0]:
    preprocessed=preprocessed.T
    if preprocessed.shape[0]!=subtypes.shape[0]:
        raise ValueError("Number of samples in preprocessed data and clinical data do not match.")

logging.info(f"Preprocessed data shape: {preprocessed.shape}")
logging.info(f"Clinical data shape: {subtypes.shape}")
logging.info(f"Stage distribution:\n{subtypes.value_counts().sort_index()}")

# ============================================================================
# Setup color mappings for UMAP visualization
# ============================================================================
possible_classes = subtypes.value_counts().sort_index().index.to_list()
dict_umap = {possible_classes[0]: "#b22222", possible_classes[1]: "#6495ed"}
logging.info(f"Color mapping for original classes: {dict_umap}")

# Create output directory for UMAP figures
date_umap = datetime.now().strftime("%Y%m%d")
umap_path = os.path.join(FIGURES_DIR, f"{date_umap}_{cancer_type}_umap_reconstruction")
os.makedirs(umap_path, exist_ok=True)
logging.info(f"UMAP figures will be saved to: {umap_path}")

# ============================================================================
# Visualize preprocessed data with UMAP
# ============================================================================
logging.info("Creating UMAP visualization of preprocessed data...")

plot_umap(
    data=preprocessed,
    clinical=subtypes,
    colors_dict=dict_umap,
    shapes_dict=None,
    n_components=2,
    save_fig=True,
    save_as=os.path.join(umap_path, "preprocessed"),
    seed=2023,
    title="Preprocessed data",
    show=False,
)
logging.info(
    f"Preprocessed data UMAP saved to: {os.path.join(umap_path, 'preprocessed')}"
)

# ============================================================================
# Load trained VAE model
# ============================================================================

if not args.hf_models:
    logging.info("Loading trained VAE model...")
    model_dir = MODELS_DIR / f"models_{cancer_type}"
    model_path = f"{model_dir}/vae/final_model.pth"
    vae_config_path = f"{model_dir}/vae/config.json"

    # Load VAE configuration from JSON
    with open(vae_config_path, "r") as f:
        vae_config_dict = json.load(f)
    logging.info(f"VAE configuration: {vae_config_dict}")

    # Initialize VAE model
    vae_config = VAEConfig()
    vae_config.INPUT_DIM = vae_config_dict["INPUT_DIM"]
    vae_config.MID_DIM = vae_config_dict["MID_DIM"]
    vae_config.LATENT_DIM = vae_config_dict["LATENT_DIM"]

    model_vae = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    ).to(device)

    # Load model weights
    checkpoint = torch.load(model_path, map_location=device, weights_only=False)
    model_vae.load_state_dict(checkpoint["model_state_dict"])
    logging.info(f"VAE model loaded from: {model_path}")

# ============================================================================
# Import VAE from Hugging Face
# ============================================================================
else:
    import huggingface_hub as hf
    logging.info("Loading VAE model from Hugging Face")
    model_dir = MODELS_DIR / "pretrained" / f"{cancer_type}"

    # Download config.json from Hugging Face
    logging.info(f"Downloading VAE config.json from Hugging Face for {cancer_type}...")
    vae_config_path = hf.hf_hub_download(
        repo_id="gprolcastelo/evenflow_models",
        filename=f"{cancer_type}/config.json",
        local_dir=model_dir / ".."
    )

    # Load VAE configuration from downloaded JSON
    with open(vae_config_path, "r") as f:
        vae_config_dict = json.load(f)
    logging.info(f"VAE configuration: {vae_config_dict}")

    # Initialize VAE model
    vae_config = VAEConfig()
    vae_config.INPUT_DIM = vae_config_dict["INPUT_DIM"]
    vae_config.MID_DIM = vae_config_dict["MID_DIM"]
    vae_config.LATENT_DIM = vae_config_dict["LATENT_DIM"]

    model_vae = VAE(
        input_dim=vae_config.INPUT_DIM,
        mid_dim=vae_config.MID_DIM,
        features=vae_config.LATENT_DIM,
    ).to(device)

    # Set filename for model based on cancer type
    if cancer_type=='KIRC':
        filename_model="KIRC/20250321_VAE_idim8516_md512_feat256mse_relu.pth"
    elif cancer_type=='BRCA':
        filename_model="BRCA/20251209_VAE_idim8954_md1024_feat512mse_relu.pth"
    else:
        raise ValueError("For Hugging Face model loading, only \'KIRC\' and \'BRCA\' cancer types are available.")

    # Download the model file
    logging.info(f"Downloading VAE model: {filename_model}...")
    model_path = hf.hf_hub_download(
        repo_id="gprolcastelo/evenflow_models",
        filename=filename_model,
        local_dir=model_dir / ".."
    )

    # Load the model
    checkpoint = torch.load(model_path, map_location=device, weights_only=False)
    model_vae.load_state_dict(checkpoint)
    logging.info("VAE model loaded from Hugging Face")

# ============================================================================
# Apply VAE to preprocessed data
# ============================================================================
logging.info("Applying VAE to preprocessed data...")
model_vae.eval()

reconstruction_x, _, _, z, scaler = apply_VAE(
    torch.tensor(preprocessed.values).to(torch.float32), model_vae, y=None
)

df_reconstruction_x = pd.DataFrame(
    reconstruction_x, index=preprocessed.index, columns=preprocessed.columns
)
df_z = pd.DataFrame(z, index=preprocessed.index)

logging.info(f"VAE reconstruction shape: {reconstruction_x.shape}")
logging.info(f"VAE latent representation shape: {z.shape}")

# Create metadata for reconstructed data
subtype_rec = pd.Series(
    [str(i) + "_rec" for i in subtypes.values],
    name="rec",
    index=[i + "_rec" for i in subtypes.index],
)
df_reconstruction_x.index = [i + "_rec" for i in df_reconstruction_x.index]

# Color mapping for reconstructed data
dict_umap_rec = {
    possible_classes[0] + "_rec": "#ffa07a",
    possible_classes[1] + "_rec": "#87ceeb",
}
logging.info(f"Color mapping for reconstructed classes: {dict_umap_rec}")

# ============================================================================
# Visualize VAE reconstruction with UMAP
# ============================================================================
logging.info("Creating UMAP visualization of VAE reconstruction...")

plot_umap(
    data=df_reconstruction_x,
    clinical=subtype_rec,
    colors_dict=dict_umap_rec,
    shapes_dict=None,
    n_components=2,
    save_fig=True,
    save_as=os.path.join(umap_path, "VAE_output"),
    seed=2023,
    title="VAE Reconstruction",
    show=False,
)
logging.info(f"VAE UMAP saved to: {os.path.join(umap_path, 'VAE_output')}")

# ============================================================================
# Load trained Reconstruction Network
# ============================================================================
if not args.hf_models:
    logging.info("Loading trained Reconstruction Network...")

    network_model_path = os.path.join(model_dir, "reconstruction_network.pth")
    networkdims_path_recnet = os.path.join(model_dir, "network_dims.csv")

    # Load network dimensions
    hyper = pd.read_csv(networkdims_path_recnet)
    hyper = hyper.values.tolist()[0]
    logging.info(f"Reconstruction Network hidden layer dimensions: {hyper}")

    # Initialize Reconstruction Network
    model_net = NetworkReconstruction(layer_dims=hyper).to(device)
    model_net.load_state_dict(
        torch.load(network_model_path, map_location=device, weights_only=False)
    )
    logging.info(f"Reconstruction Network loaded from: {network_model_path}")

# ============================================================================
# Import Reconstruction Network from Hugging Face
# ============================================================================
else:
    logging.info("Loading Reconstruction Network model from Hugging Face")

    # Download network_dims.csv from Hugging Face
    logging.info(f"Downloading network_dims.csv from Hugging Face for {cancer_type}...")
    networkdims_path_recnet = hf.hf_hub_download(
        repo_id="gprolcastelo/evenflow_models",
        filename=f"{cancer_type}/network_dims.csv",
        local_dir=model_dir / ".."
    )

    # Load network dimensions
    hyper = pd.read_csv(networkdims_path_recnet)
    hyper = hyper.values.tolist()[0]
    logging.info(f"Reconstruction Network hidden layer dimensions: {hyper}")

    # Initialize Reconstruction Network
    model_net = NetworkReconstruction(layer_dims=hyper).to(device)

    # Download the model file
    logging.info(f"Downloading Reconstruction Network model for {cancer_type}...")
    model_path_net = hf.hf_hub_download(
        repo_id="gprolcastelo/evenflow_models",
        filename=f"{cancer_type}/network_reconstruction.pth",
        local_dir=model_dir / ".."
    )

    # Load the model
    checkpoint_recnet = torch.load(model_path_net, map_location=device, weights_only=False)
    model_net.load_state_dict(checkpoint_recnet)
    logging.info("Reconstruction Network model loaded from Hugging Face")
# ============================================================================
# Apply Reconstruction Network to VAE output
# ============================================================================
logging.info("Applying Reconstruction Network to VAE output...")
model_net.eval()
rec_tensor = torch.tensor(reconstruction_x).to(torch.float32).to(device)
with torch.no_grad():
    net_output = model_net(rec_tensor)

df_net_output = pd.DataFrame(
    net_output.cpu().numpy(), index=preprocessed.index, columns=preprocessed.columns
)

logging.info(f"Reconstruction Network output shape: {net_output.shape}")

# Update index for concatenation
df_net_output.index = [i + "_rec" for i in df_net_output.index]

# ============================================================================
# Visualize Reconstruction Network output with UMAP
# ============================================================================
logging.info("Creating UMAP visualization of Reconstruction Network output...")

plot_umap(
    data=df_net_output,
    clinical=subtype_rec,
    colors_dict=dict_umap_rec,
    shapes_dict=None,
    n_components=2,
    save_fig=True,
    save_as=os.path.join(umap_path, "recnet_output"),
    seed=2023,
    title="Reconstruction Network Output",
    show=False,
)
logging.info(f"RecNet UMAP saved to: {os.path.join(umap_path, 'recnet_output')}")

# ============================================================================
# Compare original and VAE reconstruction
# ============================================================================
logging.info("Creating comparison UMAP of original vs VAE reconstruction...")

dict_umap_all = {**dict_umap, **dict_umap_rec}

data_together_vae = pd.concat([preprocessed, df_reconstruction_x], axis=0)
metadata_together = pd.concat([subtypes, subtype_rec], axis=0)

plot_umap(
    data=data_together_vae,
    clinical=metadata_together,
    colors_dict=dict_umap_all,
    shapes_dict=None,
    n_components=2,
    save_fig=True,
    save_as=os.path.join(umap_path, "preprocessed_and_vae"),
    seed=2023,
    title="Original vs VAE Reconstruction",
    show=False,
)
logging.info(
    f"Comparison UMAP (original vs VAE) saved to: {os.path.join(umap_path, 'preprocessed_and_vae')}"
)

# ============================================================================
# Compare original and Reconstruction Network output
# ============================================================================
logging.info("Creating comparison UMAP of original vs Reconstruction Network output...")

data_together_recnet = pd.concat([preprocessed, df_net_output], axis=0)

plot_umap(
    data=data_together_recnet,
    clinical=metadata_together,
    colors_dict=dict_umap_all,
    shapes_dict=None,
    n_components=2,
    save_fig=True,
    save_as=os.path.join(umap_path, "preprocessed_and_recnet"),
    seed=2023,
    title="Original vs Reconstruction Network Output",
    show=False,
)
logging.info(
    f"Comparison UMAP (original vs RecNet) saved to: {os.path.join(umap_path, 'preprocessed_and_recnet')}"
)

# ============================================================================
# Summary
# ============================================================================
logging.info("=" * 80)
logging.info("Reconstruction validation complete!")
logging.info(f"Original data shape: {preprocessed.shape}")
logging.info(f"VAE reconstruction shape: {df_reconstruction_x.shape}")
logging.info(f"RecNet reconstruction shape: {df_net_output.shape}")
logging.info(f"Latent representation shape: {df_z.shape}")
logging.info(f"All UMAP visualizations saved to: {umap_path}")
logging.info("=" * 80)

# # ============================================================================
# # Check reconstruction quality using evaluation metrics
# # ============================================================================

if not args.sdmetrics:
    logging.info("Skipping synthetic data evaluation with sdmetrics.")
    quit()

path_synth_data_evaluation = os.path.join(
    PROCESSED_DATA_DIR, f"{date_umap}_reconstruction_evaluation"
)
path_synth_plots_evaluation = os.path.join(
    FIGURES_DIR, f"{date_umap}_reconstruction_evaluation"
)
os.makedirs(path_synth_data_evaluation, exist_ok=True)
os.makedirs(path_synth_plots_evaluation, exist_ok=True)

# Load train/test split data
logging.info("Loading train/test split data for reconstruction evaluation...")
traintest_dir = INTERIM_DATA_DIR / "20251212_train_test_split"
X_train = dataset.load_rnaseq_data(os.path.join(traintest_dir, "X_train.csv"))
X_test = dataset.load_rnaseq_data(os.path.join(traintest_dir, "X_test.csv"))
y_train = subtypes[X_train.index]
y_test = subtypes[X_test.index]
logging.info(f"Train data shape: {X_train.shape}, Test data shape: {X_test.shape}")

# Run VAE on test data
logging.info("Evaluating reconstruction quality on test data...")
reconstruction_test, _, _, _, _ = apply_VAE(
    torch.tensor(X_test.values).to(torch.float32), model_vae, y=None
)
synthetic_vae = pd.DataFrame(
    reconstruction_test, index=X_test.index, columns=X_test.columns
)

# Run Reconstruction Network on VAE output
rec_tensor_test = torch.tensor(reconstruction_test).to(torch.float32).to(device)
with torch.no_grad():
    net_output_test = model_net(rec_tensor_test)
synthetic_postprocess = pd.DataFrame(
    net_output_test.cpu().numpy(), index=X_test.index, columns=X_test.columns
)

# Evaluate reconstruction quality

# Get metadata in SDMetrics format:
test_path = traintest_dir / "X_test.csv"

# On VAE output
sd_path_vae = os.path.join(path_synth_data_evaluation, "vae_reconstruction")
os.makedirs(sd_path_vae, exist_ok=True)

df_ba, df_ks = evaluate_reconstruction(
    real_data=X_test,
    synthetic_data=synthetic_vae,
    save_path_data=path_synth_data_evaluation,
    save_path_figures=sd_path_vae,
    metadata_path=test_path,
)

# On Reconstruction Network output
sd_path_recnet = os.path.join(path_synth_data_evaluation, "recnet_reconstruction")
os.makedirs(sd_path_recnet, exist_ok=True)
df_ba_rec, df_ks_rec = evaluate_reconstruction(
    real_data=X_test,
    synthetic_data=synthetic_postprocess,
    save_path_data=path_synth_data_evaluation,
    save_path_figures=sd_path_recnet,
    metadata_path=test_path,
)
