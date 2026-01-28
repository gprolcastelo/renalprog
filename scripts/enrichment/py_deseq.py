'''
This script is used to process RNASeq data and obtain a DESeq foldchange dataset for each synthetic patient.
A bash file is created for each patient. The bash file will later be used to run the GSEA analysis on the cluster.
Author: Guillermo Prol-Castelo
Date: 2024-02-21
Last modified: 2026-01-26
License: Apache 2.0
'''

import os
# Set the NUMEXPR_MAX_THREADS environment variable to be able to use more than 64 threads
os.environ['NUMEXPR_MAX_THREADS'] = '112'  # or any number greater than 64
import time
import pandas as pd
from renalprog.enrichment import fun_synth_single_patient_and_gsea
import warnings
import argparse

parser = argparse.ArgumentParser(description='Process RNASeq data and obtain GSEA enrichment.')

parser.add_argument('--cancer_type',
                    type=str,
                    choices=['kirc','lobular','ductal'],
                    help='Cancer type to analyze: "kirc" for kidney, or for breast: "lobular", "ductal"')

parser.add_argument('--traj_dir', type=str,
                    help='Path to the directory containing the synthetic trajectories'
                    )

parser.add_argument('--source_target_file', type=str,
                    help='Path to the file containing the source and target patients for the trajectories.'
                    )

parser.add_argument('--patient_stage_file', type=str,
                    help='Path to the file containing the nodes (sources and targets) metadata used in the trajectories.'
                    )

parser.add_argument('--stage_transition', type=str,
                    choices=['I_to_II','II_to_III','III_to_IV','early_to_late','early_to_early','late_to_late'],
                    help='Type of stage transition. '
                    )

parser.add_argument('--dir', type=str,
                    default = '.',
                    help='Directory to execute script at. Default is current directory.')

parser.add_argument('--data_dir',type=str,
                    default='data/interim/20230907_hard_genes_with_Maha/RNASeq.csv',
                    help='Path to the transcriptomics data file to be analyzed.')

parser.add_argument('--metadata_dir',type=str,
                    default='data/interim/20230905_preprocessed_anew/CuratedClinicalData.csv',
                    help='Path to the metadata, containing the stage')

parser.add_argument('--control_data_dir',type=str,
                    default='data/processed/controls/rnaseq_control.csv',
                    help='Path to the transcriptomics data file to be analyzed.')

parser.add_argument('--control_metadata_dir',type=str,
                    default='data/processed/controls/clinical_control.csv',
                    help='Path to the metadata, containing the stage')

# GSEA parameters
parser.add_argument('--pathway_file', type=str,
                    default='data/external/ReactomePathways.gmt',
                    help='Path to the GMT file containing gene sets for GSEA')

parser.add_argument('--mode', type=str,
                    default='Max_probe',
                    choices=['Max_probe', 'Median_of_probes', 'Mean_of_probes', 'Abs_max_of_probes'],
                    help='GSEA collapse mode for handling multiple probes per gene')

parser.add_argument('--norm', type=str,
                    default='meandiv',
                    choices=['meandiv', 'none'],
                    help='GSEA normalization mode')

parser.add_argument('--nperm', type=int,
                    default=1000,
                    help='Number of permutations for GSEA significance testing')

parser.add_argument('--rnd_seed', type=str,
                    default='timestamp',
                    help='Random seed for GSEA (use "timestamp" for random or a number for reproducibility)')

parser.add_argument('--scoring_scheme', type=str,
                    default='weighted',
                    choices=['weighted', 'classic', 'weighted_p2'],
                    help='GSEA scoring scheme for enrichment calculation')

parser.add_argument('--set_max', type=int,
                    default=500,
                    help='Maximum size of gene sets to include in GSEA')

parser.add_argument('--set_min', type=int,
                    default=15,
                    help='Minimum size of gene sets to include in GSEA')

# Add a new argument for the file name
parser.add_argument('--file', type=str, help='CSV file to process')
# Unpack arguments:
args = parser.parse_args()
cancer_type = args.cancer_type
PARALLEL_SCRIPTS_DIR = args.dir
file = args.file
# Define subpath as parent path of file:
subpath = os.path.dirname(file)
# Rename file as basename of file:
file = os.path.basename(file)
# Path to the directory containing the synthetic trajectories,
# above the transition directory (I_II or II_III):
work_dir = args.traj_dir

###############################
start_time = time.time()
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", category=RuntimeWarning, message="^invalid value encountered in log")
###############################
# State directory:
###############################
save_path_high = work_dir
os.makedirs(save_path_high, exist_ok=True)

###############################
# Load preprocessed data:
###############################


rnaseq_prepro = pd.read_csv(
    os.path.join(PARALLEL_SCRIPTS_DIR,args.data_dir),
    index_col=0
)
clinical_prepro = pd.read_csv(
    os.path.join(PARALLEL_SCRIPTS_DIR,args.metadata_dir),
    index_col=0
)

clinical_prepro = pd.DataFrame(clinical_prepro['ajcc_pathologic_tumor_stage'])

# When analyzing KIRC, change the stage to early/late
if args.cancer_type=='kirc':
    clinical_prepro.replace({'Stage I':'early','Stage II':'early','Stage III':'late','Stage IV':'late'},inplace=True)


# Make sure rnaseq_prepro is in the correct shape:
if rnaseq_prepro.shape[1]!=clinical_prepro.shape[0]:
    rnaseq_prepro=rnaseq_prepro.T
    if rnaseq_prepro.shape[1]!=clinical_prepro.shape[0]:
        raise ValueError('RNASeq data and clinical data have different number of samples')

###############################
# Load control data:
###############################
rna_control = pd.read_csv(os.path.join(PARALLEL_SCRIPTS_DIR,args.control_data_dir),index_col=0)
clinical_control = pd.read_csv(os.path.join(PARALLEL_SCRIPTS_DIR,args.control_metadata_dir),index_col=0)

# Get list of genes used:
genes_list = rnaseq_prepro.index.values
# Make sure rna_control contains the same genes:
try:
    rna_control=rna_control[genes_list]
except (KeyError, IndexError):
    print('Missing genes in control file. Recheck that this file contains at least same genes from preprocessing.')

assert rnaseq_prepro.shape[0]==rna_control.shape[1], 'RNASeq data and control data have different number of genes'

# Load source-target and patient stage data:
if cancer_type in ('lobular', 'ductal'):
    source_target = pd.read_csv(
        f'data/interim/20231004_{cancer_type}_patients_network_v2/source_target_{cancer_type}.csv',
        index_col=0
    )

    # Stages data:
    patient_stage = pd.read_csv(
        f'data/interim/20231004_{cancer_type}_patients_network_v2/stages_{cancer_type}.csv',
        index_col=0
    )

else:
    source_target = pd.read_csv(args.source_target_file)[['source', 'target']]
    patient_stage = pd.read_csv(args.patient_stage_file,)['stage']

# Loop through files in subpath:
if file.endswith('.csv'):
    # Import synthetic trajectory data
    synthetic_traj_i = pd.read_csv(
        os.path.join(subpath, file),
        index_col=0
    ).T
    # Get name of directory: I_to_II or II_to_III
    dir = os.path.basename(os.path.dirname(os.path.join(subpath, file)))
    # print(f"Subpath: {subpath}\nDir: {dir}\nFile: {file}\n")

    # print('The shape of the synthetic trajectory is: ', synthetic_traj_i.shape)
    # Loop through interpolation indexes
    for interpol_index in range(synthetic_traj_i.shape[1]):
        # time at start of iteration
        start_time_i = time.time()
        print(f"Interpol index: {interpol_index}")
        bash_file = fun_synth_single_patient_and_gsea(
            patient_here=file.split('.')[0],
            stage_trans_i=args.stage_transition,
            rnaseq_data=pd.DataFrame(synthetic_traj_i.iloc[:, interpol_index]),
            path_above_in=save_path_high,
            clinical_controls=clinical_control,
            rnaseq_controls=rna_control,
            gene_list=genes_list,
            index_pat=interpol_index,
            # GSEA parameters from command-line arguments
            pathway_file=args.pathway_file,
            mode=args.mode,
            norm=args.norm,
            nperm=args.nperm,
            rnd_seed=args.rnd_seed,
            scoring_scheme=args.scoring_scheme,
            set_max=args.set_max,
            set_min=args.set_min
        )


        path_save_bash = os.path.join(
            save_path_high,
            f'bash_file_{cancer_type}_{file.replace(".csv", "")}_{interpol_index}.cmd'
        )
        with open(
            path_save_bash,
            'w',
            encoding='utf-8',
        ) as bsh:
            bsh.writelines(bash_file)
        bsh.close()

print()
print('*'*50)
print(f'Done processing trajectory {os.path.join(subpath,file)}')
print('Total time elapsed:')
print("--- %s seconds ---" % (time.time() - start_time))
print('--- %s minutes ---' % ((time.time() - start_time) / 60))
print('--- %s hours ---' % ((time.time() - start_time) / 3600))
