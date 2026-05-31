# src/chromatopy/utils/import_data.py
from concurrent.futures import ThreadPoolExecutor
import os
import re

import pandas as pd


def read_data_concurrently(folder_path, files, signal_columns, time_column="RT (min)"):
    """
    Reads and cleans data from multiple files concurrently.

    Parameters
    ----------
    folder_path : str
        The path to the folder containing the data files.
    files : list of str
        List of filenames to be read from the folder.
    signal_columns : list of str
        List of signal identifiers to filter the data by.

    Returns
    -------
    results : list of pandas.DataFrame
        List of dataframes containing cleaned data for each file.
    """
    def load_and_clean_data(file):
        full_path = os.path.join(folder_path, file)

        # Load the entire CSV (no usecols) to check for column name variations
        df = pd.read_csv(full_path)

        # Extracting sample name from filename and storing it in the DataFrame
        df["Sample Name"] = os.path.basename(file)[:-4]

        if time_column not in df.columns and "RT(minutes) - NOT USED BY IMPORT" in df.columns:
            df.rename(columns={"RT(minutes) - NOT USED BY IMPORT": time_column}, inplace=True)

        # Remove ".0" from column names
        cleaned_columns = [col[:-2] if col.endswith(".0") else col for col in df.columns]
        df.columns = cleaned_columns

        if time_column not in df.columns:
            raise ValueError(f"Expected time column '{time_column}' in {file}")

        for signal_column in signal_columns:
            if signal_column not in df.columns:
                trace_id_with_dot_zero = signal_column + ".0"
                if trace_id_with_dot_zero in df.columns:
                    df.rename(columns={trace_id_with_dot_zero: signal_column}, inplace=True)

        required_columns = ["Sample Name", time_column] + signal_columns
        df = df[[col for col in df.columns if col in required_columns]]
        return df

    # Execute concurrently using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(load_and_clean_data, files))

    return results


def numerical_sort_key(filename):
    """
    Extracts numbers from the filename and returns them as an integer for sorting purposes.

    Parameters
    ----------
    filename : str
        The filename from which to extract the numerical values for sorting.

    Returns
    -------
    int
        The first numerical value found in the filename as an integer. If no numbers are found, returns 0.

    Notes
    -----
    - This function is typically used to sort files in numerical order based on the number(s) in their names.
    - If multiple numbers are present in the filename, only the first one is considered.
    """
    numbers = re.findall(r"\d+", filename)
    return int(numbers[0]) if numbers else 0

def import_data(results_file_path, folder_path, csv_files, signal_columns, time_column="RT (min)"):
    """
    Imports and reads the chromatographic data and existing results for HPLC analysis.
    
    This function reads in previously processed results (if available) from the specified file path. It also reads and loads the chromatographic data from CSV files located in the provided folder path and processes the trace IDs.
    
    Parameters
    ----------
    results_file_path : str
        The file path to the CSV file containing previously saved results.
    folder_path : str
        The folder path where the input CSV files from openChrom are stored.
    csv_files : list
        List of CSV files to be processed, which contain chromatographic data.
    signal_columns : list
        List of signal identifiers used to extract relevant data from the CSV files.
    
    Returns
    -------
    dict
        A dictionary containing:
        - "data" (list): The processed chromatographic data for each sample.
        - "reference" (pandas.DataFrame): The first dataset, treated as the reference sample.
        - "results_df" (pandas.DataFrame): A dataframe containing previously processed results, if available, or an empty dataframe if none exist.
    
    Notes
    -----
    - The first dataset in the data list is treated as the reference sample and is stored in the "reference" key of the returned dictionary.
    - If the results file doesn't exist at the specified path, an empty dataframe is created and returned in the "results_df".
    - The data is read concurrently using the `read_data_concurrently` function for efficient data loading.
    """
    # get or read results path
    if os.path.exists(results_file_path):
        results_df = pd.read_csv(results_file_path)
        # results_rts_df = pd.read_csv(results_rts_path)
        # results_area_unc_df = pd.read_csv(results_area_unc_path)
    else:
        results_df = pd.DataFrame(columns=["Sample Name"])
        # results_rts_df = pd.DataFrame(columns=["Sample Name"])
        # results_area_unc_df = pd.DataFrame(columns=["Sample Name"])

    print("Reading data...")
    data = read_data_concurrently(folder_path, csv_files, signal_columns, time_column=time_column)
    reference = data[0]
        
    return {
        "data": data,
        "reference": reference,
        "results_df": results_df,
        "time_column": time_column,
        "signal_columns": signal_columns,
    }
        # "results_rts_df": results_rts_df,
        # "results_area_unc_df": results_area_unc_df}


