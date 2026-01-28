from renalprog import features, dataset
from renalprog.config import RAW_DATA_DIR, INTERIM_DATA_DIR, PROCESSED_DATA_DIR
from renalprog.utils import configure_logging
import logging
import os
import json
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description='Process TCGA data with option to download preprocessed data from Zenodo.')
parser.add_argument('--zenodo_preprocessed',
                    action='store_true',
                    help='Download preprocessed data from Zenodo instead of processing raw data.')
parser.add_argument('--cancer_type',
                    type=str,
                    default='KIRC',
                    choices=['KIRC', 'BRCA'],
                    help='Cancer type to process (KIRC or BRCA). Default: KIRC')
parser.add_argument('--download_raw',
                    action='store_true',
                    help='Download raw TCGA data (otherwise assumes data is already downloaded).')
args = parser.parse_args()


# Configure logging for the pipeline
configure_logging(level=logging.INFO)

# Zenodo URLs for preprocessed data
ZENODO_URLS = {
    'BRCA': {
        'repository': 'https://zenodo.org/records/17986123',
        'rnaseq': 'https://zenodo.org/records/17986123/files/Static_BRCA.csv?download=1',
        'clinical': 'https://zenodo.org/records/17986123/files/nodes_metadata.csv?download=1'
    },
    'KIRC': {
        'repository': 'https://zenodo.org/records/17987300',
        'rnaseq': 'https://zenodo.org/records/17987300/files/Static_KIRC.csv?download=1',
        'clinical': 'https://zenodo.org/records/17987300/files/nodes_metadata.csv?download=1'
    }
}
# NOTE: Controls data is included in the repository at data/processed/controls/{CANCER_TYPE}/
# - rnaseq_control.csv
# - clinical_control.csv

cancer_type = args.cancer_type
stage_col_name = "stage" if args.zenodo_preprocessed else "ajcc_pathologic_tumor_stage"
download_data = args.download_raw

# ============================================================================
# Option 1: Download preprocessed data from Zenodo
# ============================================================================
if args.zenodo_preprocessed:
    logging.info(f"Processing preprocessed {cancer_type} data from Zenodo...")
    logging.info(f"Repository: {ZENODO_URLS[cancer_type]['repository']}")

    # Create output directory for preprocessed data
    preprocessed_data_path = INTERIM_DATA_DIR / f"preprocessed_{cancer_type}_data"

    # Download preprocessed RNAseq and clinical data from Zenodo
    preprocessed_data, clinical_data = dataset.download_preprocessed_from_zenodo(
        rnaseq_url=ZENODO_URLS[cancer_type]['rnaseq'],
        clinical_url=ZENODO_URLS[cancer_type]['clinical'],
        output_dir=preprocessed_data_path,
        rnaseq_filename="preprocessed_rnaseq.csv",
        clinical_filename="clinical_data.csv"
    )

    # ========================================================================
    # Copy control data from repository: these will be needed for enrichment analysis
    # ========================================================================
    logging.info("\n[CONTROLS DATA]")

    repo_controls_dir = PROCESSED_DATA_DIR / "controls" / cancer_type
    preprocessed_controls_dir = preprocessed_data_path / "controls"
    os.makedirs(preprocessed_controls_dir, exist_ok=True)

    # Check if controls exist in repo
    repo_rnaseq_controls = repo_controls_dir / "rnaseq_control.csv"
    repo_clinical_controls = repo_controls_dir / "clinical_control.csv"

    if repo_rnaseq_controls.exists() and repo_clinical_controls.exists():
        logging.info(f"  ✓ Control data found in repository at {repo_controls_dir}")
        logging.info(f"    - {repo_rnaseq_controls.name}")
        logging.info(f"    - {repo_clinical_controls.name}")
    else:
        logging.warning(f"  ⚠ Control data not found in repository at {repo_controls_dir}")
        logging.warning("  Expected files:")
        logging.warning(f"    - {repo_rnaseq_controls}")
        logging.warning(f"    - {repo_clinical_controls}")
        logging.warning("  Controls will need to be generated separately if needed for enrichment analysis.")

    # Summary of results
    logging.info("\n" + "=" * 80)
    logging.info("PREPROCESSED DATA DOWNLOAD COMPLETED")
    logging.info("=" * 80)
    logging.info(f"RNAseq data shape:   {preprocessed_data.shape[0]:,} genes × {preprocessed_data.shape[1]:,} samples")
    logging.info(f"Clinical data shape: {clinical_data.shape[0]:,} samples × {clinical_data.shape[1]:,} features")
    if stage_col_name in clinical_data.columns:
        logging.info("\nStage distribution:")
        for stage, count in clinical_data[stage_col_name].value_counts().sort_index().items():
            logging.info(f"  {stage}: {count}")
    logging.info("\nOutput directory: " + str(preprocessed_data_path))
    logging.info("\nTo proceed to the next step, run: python scripts/pipeline_steps/2_models.py")
    logging.info("=" * 80)

    # Exit early since we don't need to process raw data
    exit(0)

# ============================================================================
# Option 2: Process raw TCGA data
# ============================================================================
logging.info(f"Processing raw {cancer_type} data...")

## Download TCGA bulk-RNASeq dataset to the raw data directory
if download_data:
    logging.info(f"Downloading data into {RAW_DATA_DIR}")
    raw_rnaseq_path, raw_clinical_path, raw_phenotype_path = dataset.download_data(
        destination=RAW_DATA_DIR,
        remove_gz=True
    )
else:
    # To avoid repeated downloads during testing, we assume the data is already downloaded:
    raw_rnaseq_path = "data/raw/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
    raw_clinical_path = "data/raw/Survival_SupplementalTable_S1_20171025_xena_sp"
    raw_phenotype_path = "data/raw/TCGA_phenotype_denseDataOnlyDownload.tsv"

# Process the downloaded data and save to interim data directory
procesed_data_path = INTERIM_DATA_DIR / f"{cancer_type}_data"
rnaseq, clinical, pheno = dataset.process_downloaded_data(
    rnaseq_path=raw_rnaseq_path,
    clinical_path=raw_clinical_path,
    phenotype_path=raw_phenotype_path,
    cancer_type=cancer_type,
    output_dir=procesed_data_path,
    extract_controls=True,
    early_late=True
)

stages = clinical[stage_col_name]

# Pre-processing data

transpose_for_mcd = True if rnaseq.shape[0] == clinical.shape[0] else False

preprocessed_data_path = INTERIM_DATA_DIR / f"preprocessed_{cancer_type}_data"
os.makedirs(preprocessed_data_path, exist_ok=True)
preprocessed_data, preprocessing_info = features.preprocess_rnaseq(
    data=rnaseq,
    filter_expression=True,
    detect_outliers=True,
    log_transform=False,
    alpha=0.05,
    mean_threshold = 0.5,
    var_threshold = 0.5,
    min_sample_fraction = 0.2,
    support_fraction=None,
    transpose=transpose_for_mcd,
    seed = 2023
)
preprocessed_data.to_csv(preprocessed_data_path / "preprocessed_rnaseq.csv")

# Save preprocessing info as a properly formatted DataFrame
# Convert the dictionary to a single-row or single-column DataFrame
info_df = pd.DataFrame.from_dict(
    {k: [v] if not isinstance(v, list) else [str(v)] for k, v in preprocessing_info.items()}
)
info_df.to_csv(preprocessed_data_path / "preprocessing_info.csv", index=False)

# Also save as JSON for easier reading
with open(preprocessed_data_path / "preprocessing_info.json", 'w') as f:
    json.dump(preprocessing_info, f, indent=2)

stages.to_csv(preprocessed_data_path / "stages.csv")
logging.info(f"Preprocessed data saved to {preprocessed_data_path}")

# Drop duplicate genes from preprocessed data
if preprocessed_data.index.duplicated().any():
    logging.info(f"Found {preprocessed_data.index.duplicated().sum()} duplicate genes in preprocessed data, removing duplicates...")
    preprocessed_data = preprocessed_data[~preprocessed_data.index.duplicated(keep='first')]
    preprocessed_data.to_csv(preprocessed_data_path / "preprocessed_rnaseq.csv")
# Select same genes for control data
logging.info("Filtering control RNASeq genes to match preprocessed genes...")
rnaseq_controls = pd.read_csv(
    procesed_data_path/"controls"/"KIRC_control_rnaseq.csv",
    index_col=0
)
# Drop duplicate genes from control data
if rnaseq_controls.index.duplicated().any():
    logging.info(f"Found {rnaseq_controls.index.duplicated().sum()} duplicate genes in controls, removing duplicates...")
    rnaseq_controls = rnaseq_controls[~rnaseq_controls.index.duplicated(keep='first')]
logging.info(f"RNASeq controls size before filtering: {rnaseq_controls.shape[0]} genes x {rnaseq_controls.shape[1]} samples")
# Filter to match preprocessed genes
rnaseq_controls_filtered = rnaseq_controls.loc[rnaseq_controls.index.intersection(preprocessed_data.index)]
logging.info(f"RNASeq controls size after filtering: {rnaseq_controls_filtered.shape[0]} genes x {rnaseq_controls_filtered.shape[1]} samples")
os.makedirs(preprocessed_data_path / "controls", exist_ok=True)
rnaseq_controls_filtered.to_csv(preprocessed_data_path / "controls" / "preprocessed_control_rnaseq.csv")
logging.info(f"Preprocessed control data saved to {preprocessed_data_path}")

logging.info("Data processing completed.")
logging.info("To proceed to the next step, run the script: 2_models.py")
