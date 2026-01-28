#!/bin/bash

########################
# USAGE:
# bash scripts/enrichment/pipeline.sh <cancer_type> [options]
#
# Required:
#   cancer_type: Cancer type identifier (e.g., kirc)
#
# Optional (can be set as environment variables):
#   CMD_DIR: Directory containing trajectory CSV files
#   ERR_OUT_DIR: Directory for error/output logs
#   ST_FILE: Source-target file path
#   DATA_DIR: RNASeq data file path
#   METADATA_DIR: Clinical metadata file path
#   CONTROL_DATA: Control RNASeq data file path
#   CONTROL_METADATA: Control clinical metadata file path
#   META_NODES: Nodes metadata file path
#   STAGE_TRANSITION: Stage transition identifier (e.g., early_to_late)
#   PATHWAYS_FILE: Pathways GMT file (default: data/external/ReactomePathways.gmt)
#   USE_REPO_CONTROLS: Set to "true" to use controls from repository root (controls/{cancer_type}/)
#
# Example with defaults (KIRC):
#   bash scripts/enrichment/pipeline.sh kirc
#
# Example with custom paths:
#   CMD_DIR=/path/to/trajectories bash scripts/enrichment/pipeline.sh kirc
#
# Example with repo controls:
#   USE_REPO_CONTROLS=true bash scripts/enrichment/pipeline.sh kirc
#
########################

cancer_type=$1

if [ -z "$cancer_type" ]; then
    echo "Error: Cancer type is required as first argument"
    echo "Usage: bash scripts/enrichment/enrichment.sh <cancer_type>"
    exit 1
fi

########################
# INPUTS:
# Number of interpolation samples expected for each trajectory:
export N_SAMPLES=50

# Define the directory containing the parallel trajectories scripts
export PARALLEL_SCRIPTS_DIR="$(pwd)"
cd ${PARALLEL_SCRIPTS_DIR}
# get date in YYYYMMDD format:
export TODAY=$(date +"%Y%m%d")

########################
# DO NOT MODIFY REMAINING CODE
#
#
export PYTHONPATH=${PARALLEL_SCRIPTS_DIR}
# Define a function to run when the script is terminated
cleanup() {
    echo "Job was cancelled. Cleaning up..."
    # Exit the script
    exit 1
}

# Set default directories if not provided as environment variables
if [ -z "$CMD_DIR" ] || [ -z "$DATA_DIR" ] || [ -z "$METADATA_DIR" ] || [ -z "$CONTROL_DATA" ] || [ -z "$CONTROL_METADATA" ]; then
    echo "Using default paths for cancer type: ${cancer_type}"

    if [ "$cancer_type" = "kirc" ]; then
        # Default KIRC paths (Test to test trajectories)
        : ${CMD_DIR:="${PARALLEL_SCRIPTS_DIR}/data/processed/synthetic_data/kirc/recnet/early_to_late/test_to_test/"}
        : ${ERR_OUT_DIR:="${PARALLEL_SCRIPTS_DIR}/data/processed/synthetic_data/kirc/recnet/early_to_late/err_out_kirc_test"}
        : ${ST_FILE:="${PARALLEL_SCRIPTS_DIR}/data/processed/patient_trajectories_KIRC/random_connections_to_5_next_test.csv"}
        : ${DATA_DIR:="${PARALLEL_SCRIPTS_DIR}/data/interim/preprocessed_KIRC_data/KIRC_rnaseq.csv"}
        : ${METADATA_DIR:="${PARALLEL_SCRIPTS_DIR}/data/interim/preprocessed_KIRC_data/KIRC_clinical.csv"}

        # Use repo controls if specified, otherwise use preprocessed controls
        if [ "${USE_REPO_CONTROLS}" = "true" ]; then
            : ${CONTROL_DATA:="${PARALLEL_SCRIPTS_DIR}/controls/KIRC/rnaseq_controls.csv"}
            : ${CONTROL_METADATA:="${PARALLEL_SCRIPTS_DIR}/controls/KIRC/clinical_controls.csv"}
            echo "Using control data from repository root"
        else
            : ${CONTROL_DATA:="${PARALLEL_SCRIPTS_DIR}/data/interim/preprocessed_KIRC_data/controls/KIRC_control_rnaseq.csv"}
            : ${CONTROL_METADATA:="${PARALLEL_SCRIPTS_DIR}/data/interim/preprocessed_KIRC_data/controls/KIRC_control_clinical.csv"}
        fi

        : ${META_NODES:="${PARALLEL_SCRIPTS_DIR}/data/processed/patient_trajectories_KIRC/nodes_metadata.csv"}
        : ${STAGE_TRANSITION:="early_to_late"}
    else
        # Default paths for other cancer types
        : ${CMD_DIR:="${PARALLEL_SCRIPTS_DIR}/data/interim/ynthetic_data/${cancer_type}"}
        : ${ERR_OUT_DIR:="${PARALLEL_SCRIPTS_DIR}/data/interim/synthetic_data/err_out_${cancer_type}"}
        : ${DATA_DIR:="${PARALLEL_SCRIPTS_DIR}/data/interim/preprocessed_${cancer_type}_data/${cancer_type}_rnaseq.csv"}
        : ${METADATA_DIR:="${PARALLEL_SCRIPTS_DIR}/data/interim/preprocessed_${cancer_type}_data/${cancer_type}_clinical.csv"}

        # Use repo controls if specified, otherwise use preprocessed controls
        if [ "${USE_REPO_CONTROLS}" = "true" ]; then
            cancer_type_upper=$(echo ${cancer_type} | tr '[:lower:]' '[:upper:]')
            : ${CONTROL_DATA:="${PARALLEL_SCRIPTS_DIR}/controls/${cancer_type_upper}/rnaseq_controls.csv"}
            : ${CONTROL_METADATA:="${PARALLEL_SCRIPTS_DIR}/controls/${cancer_type_upper}/clinical_controls.csv"}
            echo "Using control data from repository root"
        else
            : ${CONTROL_DATA:="${PARALLEL_SCRIPTS_DIR}/data/interim/preprocessed_${cancer_type}_data/controls/${cancer_type}_control_rnaseq.csv"}
            : ${CONTROL_METADATA:="${PARALLEL_SCRIPTS_DIR}/data/interim/preprocessed_${cancer_type}_data/controls/${cancer_type}_control_clinical.csv"}
        fi
    fi
else
    echo "Using custom paths provided via environment variables"
fi

# Export all directory variables
export CMD_DIR
export ERR_OUT_DIR
export ST_FILE
export DATA_DIR
export METADATA_DIR
export CONTROL_DATA
export CONTROL_METADATA
export META_NODES
export STAGE_TRANSITION

echo "Configuration:"
echo "  Cancer type: ${cancer_type}"
echo "  CMD_DIR: ${CMD_DIR}"
echo "  ERR_OUT_DIR: ${ERR_OUT_DIR}"
echo "  ST_FILE: ${ST_FILE}"
echo "  DATA_DIR: ${DATA_DIR}"
echo "  METADATA_DIR: ${METADATA_DIR}"
echo "  CONTROL_DATA: ${CONTROL_DATA}"
echo "  CONTROL_METADATA: ${CONTROL_METADATA}"
echo "  META_NODES: ${META_NODES}"
echo "  STAGE_TRANSITION: ${STAGE_TRANSITION}"
echo ""

mkdir -p ${CMD_DIR}
mkdir -p ${ERR_OUT_DIR}

########################
# DESeq
########################

echo "Generate GREASY file for DESeq job"
if [ -f src_deseq_and_gsea_NCSR/greasy_deseq_file_${cancer_type}.txt ]; then
    rm src_deseq_and_gsea_NCSR/greasy_deseq_file_${cancer_type}.txt
fi
find ${CMD_DIR} -type f -name "*.csv" | while read file
do
    echo "python scripts/enrichment/py_deseq.py \
        --cancer_type ${cancer_type} \
        --traj_dir ${CMD_DIR} \
        --source_target_file ${ST_FILE} \
        --patient_stage_file ${META_NODES} \
        --stage_transition ${STAGE_TRANSITION} \
        --dir ${PARALLEL_SCRIPTS_DIR} \
        --data_dir ${DATA_DIR} \
        --metadata_dir ${METADATA_DIR} \
        --control_data_dir ${CONTROL_DATA} \
        --control_metadata_dir ${CONTROL_METADATA} \
        --file ${file}" >> scripts/enrichment/greasy_deseq_file_${cancer_type}.txt
done

echo "Generating DESeq data and GSEA commands on ${cancer_type} cancer type"

# Use trap to call cleanup when the script receives a SIGTERM signal
trap 'cleanup' SIGTERM

echo ""
echo "============================================"
echo "IMPORTANT: Manual execution required!"
echo "============================================"
echo "The DESeq greasy file has been generated at:"
echo "  scripts/enrichment/greasy_deseq_file_${cancer_type}.txt"
echo ""
echo "You must execute this file manually using one of the following methods:"
echo "  1. Sequential execution: bash scripts/enrichment/greasy_deseq_file_${cancer_type}.txt"
echo "  2. GNU Parallel: parallel -j <N> < scripts/enrichment/greasy_deseq_file_${cancer_type}.txt"
echo "  3. SLURM or other job scheduler (see documentation)"
echo ""
echo "Please execute the greasy file before continuing to the next step."
echo "============================================"
echo ""

########################
# GSEA
########################


# Create greasy file for gsea
# Define the directory containing the *.cmd files
export GREASY_GSEA="${CMD_DIR}/greasy_${cancer_type}.sh"

# Check if the file exists and delete it if it does
if [ -f "${GREASY_GSEA}" ]; then
    rm "${GREASY_GSEA}"
fi

# Use find to locate all *.cmd files in the specified directory
CMD_FILES=$(find ${CMD_DIR} -type f -name "*.cmd")
#Check if any *.cmd files were found
if [ -z "${CMD_FILES}" ]; then
	echo "No *.cmd files found in directory ${CMD_DIR}"
  exit 1
fi
## Concatenate the contents of all *.cmd files into a single file
echo "Concatenate all .cmd files into a single file"
for CMD_FILE in $CMD_FILES; do echo "$(cat "$CMD_FILE")"; echo ""; done >> ${GREASY_GSEA}

### Grant execution permission
chmod 750 ${GREASY_GSEA}


## Run greasy job for gsea commands
echo "starting GSEA"

## Define the function gsea to be used by GNU Parallel:
gsea() {
    "${PARALLEL_SCRIPTS_DIR}/GSEA_4.3.2/gsea-cli.sh" "$@"
}

## Export the function
export -f gsea

## Use trap to call cleanup when the script receives a SIGTERM signal
trap 'cleanup' SIGTERM

echo ""
echo "============================================"
echo "IMPORTANT: Manual execution required!"
echo "============================================"
echo "The GSEA greasy file has been generated at:"
echo "  ${GREASY_GSEA}"
echo ""
echo "You must execute this file manually using one of the following methods:"
echo ""
echo "  1. Sequential execution:"
echo "     bash ${GREASY_GSEA}"
echo ""
echo "  2. GNU Parallel (recommended for faster execution):"
echo "     parallel -j <N> < ${GREASY_GSEA}"
echo "     Example with 8 parallel jobs: parallel -j 8 < ${GREASY_GSEA}"
echo ""
echo "  3. SLURM job scheduler:"
echo "     See documentation for SLURM submission examples"
echo ""
echo "Please execute the greasy file before continuing to the final step."
echo "============================================"
echo ""

#########################
## Final dataset
#########################
# Set default pathways file if not provided
: ${PATHWAYS_FILE:="data/external/ReactomePathways.gmt"}
export PATHWAYS_FILE

parent_dir=$(dirname "$CMD_DIR")

## Use trap to call cleanup when the script receives a SIGTERM signal
trap 'cleanup' SIGTERM

echo ""
echo "============================================"
echo "IMPORTANT: Manual execution required!"
echo "============================================"
echo "After completing the GSEA analysis, you need to run the trajectory formatting step."
echo ""
echo "Configuration:"
echo "  Parent directory: ${parent_dir}"
echo "  Pathways file: ${PATHWAYS_FILE}"
echo "  CMD_DIR: ${CMD_DIR}"
echo "  Cancer type: ${cancer_type}"
echo ""
echo "Execute the following command:"
echo ""
echo "  python scripts/enrichment/trajectory_formatting.py \\"
echo "    --path_synth ${CMD_DIR} \\"
echo "    --pathways_file ${PATHWAYS_FILE} \\"
echo "    --save_dir ${parent_dir} \\"
echo "    --cancer_type ${cancer_type}"
echo ""
echo "This will create the final trajectory dataset with enrichment results."
echo "Results will be saved to: ${parent_dir}"
echo "============================================"
echo ""
echo "Pipeline setup complete! Follow the instructions above to complete the enrichment analysis."

