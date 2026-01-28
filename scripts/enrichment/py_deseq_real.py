'''
This script is used to process RNASeq data and obtain a DESeq foldchange dataset for each real patient.
A bash file is created for each patient. The bash file will later be used to run the GSEA analysis on the cluster.
Author: Guillermo Prol-Castelo
Date: 2024-02-21
Last modified: 2024-09-30
License: Apache 2.0
'''

import os
# Set the NUMEXPR_MAX_THREADS environment variable to be able to use more than 64 threads
os.environ['NUMEXPR_MAX_THREADS'] = '112'
import time
import pandas as pd
# from tqdm import tqdm
from renalprog.enrichment import fun_single_patient_and_gsea
import warnings
warnings.filterwarnings('ignore')
# import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="^invalid value encountered in log")



def main(args):
    ductal_or_lobular = args.DUCTAL_OR_LOBULAR
    PARALLEL_SCRIPTS_DIR = args.dir
    # file = args.file
    # Define subpath as parent path of file:
    # subpath = os.path.dirname(file)
    # Rename file as basename of file:
    # file = os.path.basename(file)
    # Path to the directory containing the synthetic trajectories,
    # above the transition directory (I_II or II_III):
    # work_dir = os.path.join(args.traj_dir, ductal_or_lobular)
    ###############################
    # print('Begin processing file: ', file)
    start_time = time.time()


    ###############################
    # State directory:
    ###############################

    os.makedirs(args.save_path, exist_ok=True)

    if ductal_or_lobular is None:
        save_path = args.save_path
    else:
        save_path = os.path.join(args.save_path,ductal_or_lobular)
        os.makedirs(save_path, exist_ok=True)

    ###############################
    # Load preprocessed data:
    ###############################
    print('Load preprocessed data')
    rnaseq_prepro = pd.read_csv(
        filepath_or_buffer=os.path.join(PARALLEL_SCRIPTS_DIR, args.rnaseq_path),
        index_col=0,
    )
    clinical_prepro = pd.read_csv(
        filepath_or_buffer=os.path.join(PARALLEL_SCRIPTS_DIR, args.clinical_path),
        index_col=0,
    )
    clinical_prepro = clinical_prepro['ajcc_pathologic_tumor_stage']
    # Make sure samples are rows:
    if rnaseq_prepro.shape[0] != clinical_prepro.shape[0]:
        rnaseq_prepro = rnaseq_prepro.T
        if rnaseq_prepro.shape[0] != clinical_prepro.shape[0]:
            raise ValueError('RNASeq preprocessed and clinical preprocessed have incompatible numbers of samples.')

    ###############################
    # Load control data:
    ###############################
    print('Load control data')
    rnaseq_controls = pd.read_csv(
        filepath_or_buffer=os.path.join(PARALLEL_SCRIPTS_DIR, args.rnaseq_control),
        index_col=0,
    )
    clinical_controls = pd.read_csv(
        filepath_or_buffer=os.path.join(PARALLEL_SCRIPTS_DIR, args.clinical_control),
        index_col=0,
    )
    # Make sure samples are rows:
    if rnaseq_controls.shape[0] != clinical_controls.shape[0]:
        rnaseq_controls = rnaseq_controls.T
        if rnaseq_controls.shape[0] != clinical_controls.shape[0]:
            raise ValueError('RNASeq controls and clinical controls have incompatible numbers of samples.')

    print('rnaseq_prepro shape:', rnaseq_prepro.shape)
    print('rnaseq_controls shape:', rnaseq_controls.shape)
    print('clinical_controls shape:', clinical_controls.shape)
    print('clinical_prepro shape:', clinical_prepro.shape)
    # Make sure columns are shared between preprocessed and control RNASeq data:
    shared_columns = rnaseq_prepro.columns.intersection(rnaseq_controls.columns)
    rnaseq_prepro = rnaseq_prepro[shared_columns]
    rnaseq_controls = rnaseq_controls[shared_columns]


    # Get list of genes used:
    gene_list = rnaseq_prepro.columns.values
    # Get list of patients
    list_of_patients = rnaseq_prepro.index.to_list()
    rnaseq_prepro=rnaseq_prepro.T

    print(gene_list.shape,
          rnaseq_controls.shape,
          rnaseq_prepro.shape,
          clinical_controls.shape
          )

    for patient_i in list_of_patients:
        start_time_i = time.time()
        stage_i = clinical_prepro.loc[patient_i]
        if stage_i.startswith('Stage'):
            stage_i = stage_i.split(' ')[1]
        # print(f'Patient: {patient_i}\nStage: {stage_i}')
        bash_file = fun_single_patient_and_gsea(
            patient_here=patient_i,
            patient_stage=stage_i,
            rnaseq_data=rnaseq_prepro.T,
            path_above_in=save_path,
            clinical_controls=clinical_controls,
            rnaseq_controls=rnaseq_controls,
            gene_list=gene_list,
            foldchange=True,
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
            save_path,
            f'bash_file_{ductal_or_lobular}_{patient_i.replace(".csv", "")}.cmd'
        )
        print(f'Saving bash file to {path_save_bash}')
        print()
        with open(
                path_save_bash,
                'w',
                encoding='utf-8',
        ) as bsh:
            bsh.writelines(bash_file)
        bsh.close()

        print(f'Time elapsed for patient {patient_i}:')
        print("--- %s seconds ---" % (time.time() - start_time_i))
        print('--- %s minutes ---' % ((time.time() - start_time_i) / 60))
        print('--- %s hours ---' % ((time.time() - start_time_i) / 3600))

    print()
    print('*'*50)
    # print(f'Done processing trajectory {os.path.join(subpath,file)}')
    print('Total time elapsed:')
    print("--- %s seconds ---" % (time.time() - start_time))
    print('--- %s minutes ---' % ((time.time() - start_time) / 60))
    print('--- %s hours ---' % ((time.time() - start_time) / 3600))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Process RNASeq data and obtain GSEA enrichment.')

    parser.add_argument('--rnaseq_path',
                        type=str,
                        help='Path to the preprocessed RNASeq data',
                        default='data/interim/20230907_hard_genes_with_Maha/RNASeq.csv')

    parser.add_argument('--clinical_path',
                        type=str,
                        help='Path to the preprocessed clinical data',
                        default='data/interim/20230905_preprocessed_anew/CuratedClinicalData.csv')

    parser.add_argument('--rnaseq_control',
                        type=str,
                        help='Path to the control RNASeq data',
                        default='data/processed/controls/rnaseq_control.csv')

    parser.add_argument('--clinical_control',
                        type=str,
                        help='Path to the control clinical data',
                        default='data/processed/controls/clinical_control.csv')

    parser.add_argument('--cancer_type',
                        type=str,
                        help='Cancer type to analyze',
                        default='BRCA',
                        choices=['BRCA', 'KIRC'])

    parser.add_argument(
        '--DUCTAL_OR_LOBULAR',
        type=str,
        help='If the cancer type is BRCA choose which ductal/lobular to analyze',
        default=None)

    parser.add_argument('--save_path', type=str,
                        help='Directory to save results',
                        )

    parser.add_argument('--today',
                        type=str,
                        help='Experiment\'s beginning date')

    parser.add_argument('--dir', type=str,
                        default='.',
                        help='Directory to execute script at. Default is current directory.')

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
    # parser.add_argument('--file', type=str, help='CSV file to process')
    # Unpack arguments:
    args = parser.parse_args()
    main(args)

# Example for KIRC:
# python src_deseq_and_gsea_NCSR/py_deseq_real.py --rnaseq_path data/interim/20240930_preprocessed_KIRC/rnaseq_maha.csv \
# --clinical_path data/interim/20240930_preprocessed_KIRC/CuratedClinicalData.csv \
# --rnaseq_control data/processed/controls/KIRC/rnaseq_control.csv \
# --clinical_control data/processed/controls/KIRC/clinical_control.csv \
# --cancer_type KIRC \
# --save_path data/interim/20251210_KIRC_real_enrichments