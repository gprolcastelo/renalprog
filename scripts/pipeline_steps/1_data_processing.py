from renalprog import features, dataset
from renalprog.config import RAW_DATA_DIR, INTERIM_DATA_DIR
from renalprog.utils import configure_logging
import logging
import os
import json
import pandas as pd

# Configure logging for the pipeline
configure_logging(level=logging.INFO)

cancer_type = "KIRC"

## Download TCGA bulk-RNASeq dataset to the raw data directory
# raw_rnaseq_path, raw_clinical_path, raw_phenotype_path = dataset.download_data(
#     destination=RAW_DATA_DIR,
#     remove_gz=True
# )
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
    early_late = True
)

stages = clinical['ajcc_pathologic_tumor_stage']

# Pre-processing data

transpose_for_mcd = True if rnaseq.shape[0] == clinical.shape[0] else False

preprocessed_data_path = INTERIM_DATA_DIR / f"preprocessed_{cancer_type}_data"
os.makedirs(preprocessed_data_path, exist_ok=True)
preprocessed_data, preprocessing_info = features.preprocess_rnaseq(
    data=rnaseq,
    filter_expression=True,
    detect_outliers=True,
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

