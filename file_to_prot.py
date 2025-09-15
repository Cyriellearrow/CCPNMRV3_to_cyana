#!/usr/bin/python3
import pandas as pd
import os
import sys
import monentlature


def read_prot_seq(filepath):
    """Read prot.seq file and create a residue mapping dictionary."""
    residue_map = {}
    try:
        with open(filepath, 'r') as file:
            content = file.read()
        
        for line in content.strip().split('\n'):
            if line.strip():
                parts = line.split()
                if len(parts) >= 2:
                    residue, seq_code = parts[0].upper(), parts[1]
                    residue_map[seq_code] = residue
        return residue_map
    
    except FileNotFoundError:
        print(f"Warning: {filepath} not found. Proceeding without residue mapping.")
        return {}
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return {}

def transform_sequence_code(code):
    """Transform SequenceCode by decrementing if it contains '-1'."""
    code_str = str(code)
    if '-1' in code_str:
        try:
            number = int(code_str.split('-')[0]) - 1
            return str(number)
        except ValueError:
            return code_str
    return code_str

def update_residue_type(row, residue_map):
    """Update ResidueType with uppercase three-letter code from residue_map if None."""
    if pd.isna(row['ResidueType']) or str(row['ResidueType']).lower() in ['none', 'nan']:
        seq_code = str(row['SequenceCode'])
        return residue_map.get(seq_code, row['ResidueType'])
    return str(row['ResidueType']).upper()

def process_protein_data(input_csv, prot_seq_filepath, output_csv):
    """Process protein NMR data with deduplication and residue mapping."""

    
    # Read the prot.seq file
    residue_map = read_prot_seq(prot_seq_filepath)
    
    # Read the CSV file
    df = pd.read_csv(input_csv)
    
    # Clean column names (remove newlines and extra spaces)
    df.columns = df.columns.str.replace('\n', ' ').str.strip()
    
    # Verify and filter rows with non-zero and non-NaN Total Peak Count
    if 'Total Peak Count' in df.columns:
        initial_rows = len(df)
        df = df[(df['Total Peak Count'] != 0) & (pd.notna(df['Total Peak Count']))]
        filtered_rows = len(df)
    
    # Transform SequenceCode
    df['SequenceCode'] = df['SequenceCode'].apply(transform_sequence_code)
    
    # Update ResidueType with uppercase three-letter codes
    df['ResidueType'] = df.apply(lambda row: update_residue_type(row, residue_map), axis=1)
    
    # Convert columns to numeric
    numeric_columns = ['Value (ppm)', 'Value Error (ppm)', 'Total Peak Count', 'uniqueId']
    
    for col in numeric_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            nan_count = df[col].isnull().sum()
    
    
    # Group by SequenceCode and AtomName
    df_grouped = df.groupby(['SequenceCode', 'AtomName']).agg({
        'Value (ppm)': 'mean',
        'Value Error (ppm)': 'mean',
        'ResidueType': 'first',
        'uniqueId': 'min',
        'Total Peak Count': 'sum'
    }).reset_index()
    
    
    # Ensure output column names match the desired structure
    column_order = ['uniqueId', 'Value (ppm)', 'Value Error (ppm)', 
                   'SequenceCode', 'ResidueType', 'AtomName', 'Total Peak Count']
    
    # Check if all required columns exist
    missing_cols = [col for col in column_order if col not in df_grouped.columns]
    if missing_cols:
        # Use only available columns
        column_order = [col for col in column_order if col in df_grouped.columns]
    
    df_final = df_grouped[column_order]
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Save the modified file
    df_final.to_csv(output_csv, index=False)
    
    # Summary statistics
    original_rows = len(df)
    final_rows = len(df_final)
    duplicates_removed = original_rows - final_rows
    

    
    return df_final

def additional_processing(df, output_file):
    """Additional processing step to prepare final output files."""
    
    # Create a copy for processing
    data_zero = df.copy()
    
    # Define columns to delete
    columns_to_delete = ['uniqueId', 'NmrAtom', 'ResidueType', 'Total Peak Count']
    
    # Drop specified columns if they exist
    data_columns = data_zero.drop(columns=[col for col in columns_to_delete if col in data_zero.columns])
    

    # Remove all attributions with -1 (duplicates)
    initial_rows = len(data_columns)
    data_columns = data_columns[~data_columns['SequenceCode'].astype(str).str.contains('-1')]
    filtered_rows = len(data_columns)
    
    # Create a new column with incremental numbers, beginning by 1
    data_columns.insert(0, 'Index', range(1, len(data_columns) + 1))
    
    # Reorder the columns - handle both cleaned and original column names
    possible_value_cols = ['Value (ppm)', 'Value\n(ppm)']
    possible_error_cols = ['Value Error (ppm)', 'Value Error\n(ppm)']
    
    value_col = next((col for col in possible_value_cols if col in data_columns.columns), None)
    error_col = next((col for col in possible_error_cols if col in data_columns.columns), None)
    
    # Use the format with newlines for consistency
    colonnes_ordre = ['Index', 'Value\n(ppm)', 'Value Error\n(ppm)', 'AtomName', 'SequenceCode']
    
    # Rename columns to match expected format if needed
    if value_col != 'Value\n(ppm)':
        data_columns = data_columns.rename(columns={value_col: 'Value\n(ppm)'})
    if error_col != 'Value Error\n(ppm)':
        data_columns = data_columns.rename(columns={error_col: 'Value Error\n(ppm)'})
    
    # Select and reorder columns
    available_columns = [col for col in colonnes_ordre if col in data_columns.columns]
    data_columns = data_columns[available_columns]
    
    # Replace all "None" values with 0
    data_columns = data_columns.replace("None", 0)
    data_columns = data_columns.replace("none", 0)
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_file):
        os.makedirs(output_file)
    
    # Save the file as .prot format with tab separation
    prot_file_path = os.path.join(output_file, 'attib_cyana.prot')
    data_columns.to_csv(prot_file_path, sep='\t', index=False)
    
    return data_columns


def lunch_all (input_csv, prot_seq_file,output_directory, output_csv='attib_cyana.prot'):
    # Process the data
    result = process_protein_data(input_csv, prot_seq_file, output_csv)

    if result is not None:

        # Additional processing
        final_result = additional_processing(result, output_directory)

        if final_result is not None:

            seq_file = prot_seq_file
            name_prot_file = os.path.join(output_directory, 'attib_cyana.prot')
            output_cya_file = output_csv

            trad = monentlature.nomenclature(seq=seq_file, prot=name_prot_file, output_file=output_directory)
            monentlature.file_transforme(trad, prot=name_prot_file, output_file=output_cya_file)


def main():
    cyana_attrib_input = sys.argv[1]
    seq_file = sys.argv[2]
    save_path = sys.argv[3]
    output_csv = sys.argv[4]
    
    lunch_all(cyana_attrib_input,seq_file,save_path,output_csv)

if __name__ == "__main__":
    main()
