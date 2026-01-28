'''
This script processes all GSEA enrichment data sets and saves the results to a single csv file.
Author: Guillermo Prol-Castelo
Date: 2024-02-20
Last modified: 2025-02-25
License: Apache 2.0
'''

# Example run:
# python src_deseq_and_gsea_NCSR/trajectory_formatting.py --path_synth data/interim/20250219_synthetic_data/lobular \
# --pathways_file data/external/ReactomePathways.gmt \
# --save_dir data/interim/20250219_synthetic_data \
# --cancer_type lobular

# Import necessary libraries
import os
import glob
import re
import pandas as pd
import numpy as np
from tqdm import tqdm

def find_file(path_to_file:str)->pd.DataFrame:
    """
    Find a file in a given directory, and return a DataFrame with its contents.
    :param path_to_file:
    :return: df, pd.DataFrame of the file contents
    """
    # Check if the file was found
    if path_to_file:
        # Read the TSV file
        df = pd.read_csv(path_to_file[0], sep='\t')
    else:
        # Raise an error if the file was not found
        raise FileNotFoundError(f'File starting with "gsea_report_for_na_" not found in {path_to_file}')
    return df

def read_gsea_reports(path_here:str,patient_here:str,interpol_index:int, transition:str)->pd.DataFrame:
    """
    Read GSEA reports from a given path, and return a combined dataframe with the positive and negative reports.
    :param path_here: the path to the GSEA reports in tsv format
    :param patient_here: string with the name of the patient
    :param interpol_index: the index of the interpolation
    :param transition: transition of stage from source to target patient (e.g., '1_to_2')
    :return: combined_report, pd.DataFrame with the combined positive and negative reports and columns NAME, ES, NES, FDR q-val, Patient, Idx, Transition
    """

    # Get positive and negative reports
    ## Positive is a tsv that starts with 'gsea_report_for_na_pos
    # Find the file that starts with 'gsea_report_for_na_pos'
    file_pos = glob.glob(os.path.join(path_here, 'gsea_report_for_na_pos*.tsv'))
    file_neg = glob.glob(os.path.join(path_here, 'gsea_report_for_na_neg*.tsv'))

    positive_report=find_file(file_pos)
    ## Negative is a tsv that starts with 'gsea_report_for_na_neg'
    negative_report=find_file(file_neg)
    # Keep only relevant columns
    positive_report = positive_report[['NAME','ES','NES','FDR q-val']]
    negative_report = negative_report[['NAME','ES','NES','FDR q-val']]
    # Combine positive and negative reports
    combined_report=pd.concat([positive_report,negative_report],axis=0)
    # Add patient, interpolation index, and transition columns
    combined_report.insert(0,'Patient',patient_here)
    combined_report.insert(1,'Idx',interpol_index)
    combined_report.insert(2,'Transition',transition)

    return combined_report

def add_missing_pathways(patient:str,transition:str, idx:int, report:pd.DataFrame,pathways:list)->pd.DataFrame:
    """
    Add missing pathways to the report, with NaN values for the numeric columns.
    :param patient: str, name of the patient
    :param transition: str, transition of stage from source to target patient (e.g., '1_to_2')
    :param idx: int, index of the interpolation
    :param report: pd.DataFrame, the GSEA report imported from the tsv files with columns Patient, Idx, Transition, NAME, ES, NES, FDR q-val
    :param pathways: list, list of all considered pathways
    :return: pd.DataFrame of report including the missing pathways
    """
    missing_pathways = np.array(pathways)[~np.isin(pathways,report['NAME'])].tolist()
    if missing_pathways.__len__()==0:
        print('No missing pathways')
        return report
    else:
        missing_report = pd.DataFrame(
        data={
                'Patient':patient,
                'Idx':idx,
                'Transition':transition,
                'NAME':missing_pathways,
                'ES':np.nan,
                'NES':np.nan,
                'FDR q-val':np.nan
            }
        )

        return pd.concat([report,missing_report],axis=0)


def gsea_interpolation_concat(interpol_dir:str,transition:str,patient_here:str,pathways:list)->pd.DataFrame:
    """
    Concatenate the GSEA reports in a given directory.
    :param interpol_dir: str, path to the directory with the GSEA reports
    :param transition: str, the label of the transition (e.g., '1_to_2')
    :param patient_here: str, the name of the patient
    :param pathways: list, list of all considered pathways
    :return: pd.DataFrame with the concatenated GSEA reports, including missing pathways with NaN values for the numeric columns
    """
    interpol_i_df = pd.DataFrame()
    # Loop through child directories, which correspond to the interpolations
    for dir_i in os.listdir(interpol_dir):
        # Get index of interpolation from directory name
        try:
            idx_i = int(re.split(r'[_\.]', dir_i)[3])
        except IndexError:
            # If no indices found, this means we're analyzing real patients, not interpolations
            idx_i = 0
        # Read GSEA reports
        combined_report = read_gsea_reports(
            path_here=os.path.join(interpol_dir,dir_i),
            patient_here=patient_here,
            interpol_index=idx_i,
            transition=transition
        )
        # print('combined report from gsea:')
        # print(idx_i,combined_report.shape)
        # Concatenate the missing pathways into the report
        full_report = add_missing_pathways(
            patient=patient_here,
            transition=transition,
            idx=idx_i,
            report=combined_report,
            pathways=pathways)
        # Sort by pathway name and reset index
        full_report.sort_values(by='NAME',inplace=True,ignore_index=True)
        interpol_i_df=pd.concat([interpol_i_df,full_report],axis=0,ignore_index=True)

    return interpol_i_df


def gsea_patient_concat(path_transition:str,transition:str,pathways:list)->pd.DataFrame:
    """
    Loop through patients at a given transition, and concatenate the GSEA reports.
    :param path_transition: path to the transition
    :param transition: transition name
    :param pathways: list of all considered pathways
    :return: patient_df, pd.DataFrame with the combined GSEA reports
    """
    patient_df = pd.DataFrame()
    for patient_i in tqdm(os.listdir(path_transition),
                          desc=f'Loop through patients at transition {transition}'):
        path_i = os.path.join(path_transition, patient_i)
        if os.path.isdir(path_i):
            patient_df_i = gsea_interpolation_concat(
                interpol_dir=os.path.join(path_i, 'reports'),
                transition=transition,
                patient_here=patient_i,
                pathways=pathways,
            )
            patient_df = pd.concat([patient_df, patient_df_i], axis=0,ignore_index=True)

    return patient_df

def get_regulation(nes, q_val, thres=0.05, alpha=0.05):
    """
    Determine regulation based on NES and q-value.
    1 = upregulated
    -1 = downregulated
    0 = not regulated
    NaN = not significant
    :param nes:
    :param q_val:
    :param thres:
    :param alpha:
    :return:
    """
    if pd.isna(nes) or nes == '---':
        return np.nan
    elif float(nes) > thres and float(q_val) < alpha:
        return 1
    elif float(nes) < thres and float(q_val) < alpha:
        return -1
    elif float(q_val) >= alpha:
        return np.nan
    else:
        return 0

def get_significant_nes(nes, q_val, alpha=0.05):
    """
    Return NES if significant, else NaN.
    :param nes:
    :param q_val:
    :param alpha:
    :return: NES if significant, else NaN
    """
    if pd.isna(nes) or nes == '---':
        return np.nan
    elif float(q_val) < alpha:
        return nes
    elif float(q_val) >= alpha:
        return np.nan
    else:
        return 0




def main(path_synth:str,cancer_type:str,pathways_file:str,save_dir:str):
    with open(pathways_file, 'r') as file:
        # Convert pathway names to uppercase and store in 'pathways' list
        pathways = [line.split('\t')[0].upper() for line in file]
    # Loop through all available transitions
    transition_df = pd.DataFrame()
    for transition_i in os.listdir(path_synth):
        path_transition_i = os.path.join(path_synth, transition_i)
        if os.path.isdir(path_transition_i):
            patient_df = gsea_patient_concat(path_transition=path_transition_i, transition=transition_i,pathways=pathways)
            transition_df = pd.concat([transition_df, patient_df], axis=0, ignore_index=True)
    # Sort by transition, patient, index, and pathway name
    transition_df.sort_values(
        by=['Transition', 'Patient', 'Idx', 'NAME'],
        inplace=True,
        ignore_index=True
    )

    # Make sure ES, NES, and FDR q-val are numeric
    transition_df['ES'] = pd.to_numeric(transition_df['ES'], errors='coerce', downcast=None)
    transition_df['NES'] = pd.to_numeric(transition_df['NES'], errors='coerce', downcast=None)
    transition_df['FDR q-val'] = pd.to_numeric(transition_df['FDR q-val'], errors='coerce', downcast=None)

    # Rows = number of pathways * number of patients * number of interpolations = 2630 * 1225 * 50 = 161,087,500 in train set
    # Rows = number of pathways * number of patients * number of interpolations = 2630 * 269 * 50 = 35,373,500 in test set
    print("Shape of final dataset:",transition_df.shape)

    # Add a source and a target patient columns when analyzing synthetic trajectories:
    if "_to_" in transition_df['Patient'].iloc[0]:
        st_cols = pd.DataFrame([i.split("_to_") for i in transition_df['Patient']], columns=['Source', 'Target'])
        transition_df = pd.concat([st_cols, transition_df], axis=1)

    # Add a Regulation column based on NES and FDR q-val
    transition_df["Regulation"] = [get_regulation(i, j, thres=0.05, alpha=0.05) for i, j in
                                    zip(transition_df['NES'], transition_df['FDR q-val'])]
    # Add a significant NES column based on NES and FDR q-val
    transition_df["significant_NES"] = [get_significant_nes(i, j, alpha=0.05) for i, j in
                                    zip(transition_df['NES'], transition_df['FDR q-val'])]

    # Add TrajectoryID
    transition_df['TrajectoryID'] = pd.factorize(transition_df['Patient'])[0]

    # Save the concatenated GSEA reports to a CSV file
    save_as = os.path.join(save_dir, f'full_gsea_reports_{cancer_type}.csv')
    transition_df.to_csv(save_as, index=False)
    print(f'Report saved to {save_as}')

    ## --> Get heatmap for all pathways
    # Get unique pathways
    unique_pathways = transition_df["NAME"].unique()
    n_unique_pathways = unique_pathways.shape[0]
    print('Number of unique pathways:', n_unique_pathways)
    if "_to_" in transition_df['Patient'].iloc[0]:
        n_timepoints = transition_df["Idx"].max(skipna=True)+1
        print('Number of timepoints:', n_timepoints)
    else: # For real patients, get number of stages from unique transitions
        n_timepoints = transition_df["Transition"].nunique()
        print('Number of stages:', n_timepoints)
    heatmap_regulation = np.zeros((n_timepoints, n_unique_pathways)).T
    heatmap_nes = np.zeros((n_timepoints, n_unique_pathways)).T
    if "_to_" in transition_df['Patient'].iloc[0]:
        for traj_i in transition_df['TrajectoryID'].unique():
            heatmap_regulation = np.nansum(
                [heatmap_regulation,
                 transition_df[(transition_df['TrajectoryID'] == traj_i)]['Regulation'].values.reshape(n_timepoints, n_unique_pathways).T],
                axis=0)
            heatmap_nes = np.nansum(
                [heatmap_nes,
                 transition_df[(transition_df['TrajectoryID'] == traj_i)]['significant_NES'].values.reshape(n_timepoints, n_unique_pathways).T],
                axis=0)
    else:
        for transition_i, stage in enumerate(transition_df['Transition'].unique()):
            stage_data = transition_df[transition_df['Transition'] == stage]
            # Group by pathway and aggregate across patients at this stage
            stage_regulation = stage_data.groupby('NAME')['Regulation'].sum().values
            stage_nes = stage_data.groupby('NAME')['significant_NES'].sum().values

            # Add to the appropriate column in the heatmap
            heatmap_regulation[:, transition_i] = stage_regulation
            heatmap_nes[:, transition_i] = stage_nes
    print("heatmap_regulation.shape:",heatmap_regulation.shape)
    print("heatmap_nes.shape:",heatmap_nes.shape)
    # Create DataFrames for heatmaps
    heatmap_regulation_df = pd.DataFrame(
        data=heatmap_regulation,
        index=unique_pathways,
    )
    heatmap_nes_df = pd.DataFrame(
        data=heatmap_nes,
        index=unique_pathways,
    )
    # Save the heatmaps to csv files
    heatmap_regulation_df.to_csv(os.path.join(save_dir, f'heatmap_{cancer_type}_Regulation.csv'))
    heatmap_nes_df.to_csv(os.path.join(save_dir, f'heatmap_{cancer_type}_significantNES.csv'))
    print(f'Heatmaps saved to {save_dir}')

if __name__ == '__main__':
    import argparse

    # ArgumentParser object for command-line option and argument parsing.
    parser = argparse.ArgumentParser(description='Process all GSEA enrichment data sets.')

    # --path_synth argument: Specifies the path to the synthetic data, which contains the transitions, then the patients, and (within a reports folder) finally the interpolations.
    parser.add_argument('--path_synth', type=str,
                        help='Path to the synthetic data. No default.')

    # --pathways_file argument: Specifies the path to the file with the list of pathways.
    parser.add_argument('--pathways_file', type=str,
                        help='Path to the file with the list of pathways. No default.')
    # --save_dir argument: Specifies the working directory.
    parser.add_argument('--save_dir', type=str,
                        default='results',
                        help='Directory to save processed data. Default is ./results/')

    # --cancer_type argument: Specifies the type of cancer to process.
    parser.add_argument('--cancer_type', type=str,
                        choices=['lobular', 'ductal', 'negative', 'random', 'kirc'],
                        help='Cancer type to process. Choices: lobular or ductal for BRCA, kirc for KIRC,. No default.')


    args = parser.parse_args()

    # Call the main function
    main(
        path_synth=args.path_synth,
        cancer_type=args.cancer_type,
        pathways_file=args.pathways_file,
        save_dir=args.save_dir
    )