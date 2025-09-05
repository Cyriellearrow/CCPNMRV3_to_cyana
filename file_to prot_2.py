#!/usr/bin/python3
import pandas as pd
import numpy as np
import os
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
        
        print(f"Loaded {len(residue_map)} residue mappings from {filepath}")
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
    
    # Check if input file exists
    if not os.path.exists(input_csv):
        print(f"Error: Input file {input_csv} not found.")
        return
    
    # Read the prot.seq file
    residue_map = read_prot_seq(prot_seq_filepath)
    
    # Read the CSV file
    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)
    print(f"Original data shape: {df.shape}")
    print("Original columns:", df.columns.tolist())
    
    # Clean column names (remove newlines and extra spaces)
    df.columns = df.columns.str.replace('\n', ' ').str.strip()
    print("Cleaned columns:", df.columns.tolist())
    
    # Verify and filter rows with non-zero and non-NaN Total Peak Count
    if 'Total Peak Count' in df.columns:
        initial_rows = len(df)
        df = df[(df['Total Peak Count'] != 0) & (pd.notna(df['Total Peak Count']))]
        filtered_rows = len(df)
        print(f"Filtered out {initial_rows - filtered_rows} rows with zero/NaN Total Peak Count")
    
    # Transform SequenceCode
    print("Transforming SequenceCode...")
    df['SequenceCode'] = df['SequenceCode'].apply(transform_sequence_code)
    
    # Update ResidueType with uppercase three-letter codes
    print("Updating ResidueType...")
    df['ResidueType'] = df.apply(lambda row: update_residue_type(row, residue_map), axis=1)
    
    # Convert columns to numeric
    print("Converting columns to numeric...")
    numeric_columns = ['Value (ppm)', 'Value Error (ppm)', 'Total Peak Count', 'uniqueId']
    
    for col in numeric_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            nan_count = df[col].isnull().sum()
            if nan_count > 0:
                print(f"Warning: {nan_count} NaN values in {col} after conversion")
    
    print(f"Data shape before grouping: {df.shape}")
    
    # Group by SequenceCode and AtomName
    print("Grouping and aggregating data...")
    df_grouped = df.groupby(['SequenceCode', 'AtomName']).agg({
        'Value (ppm)': 'mean',
        'Value Error (ppm)': 'mean',
        'ResidueType': 'first',
        'uniqueId': 'min',
        'Total Peak Count': 'sum'
    }).reset_index()
    
    print(f"Data shape after grouping: {df_grouped.shape}")
    
    # Ensure output column names match the desired structure
    column_order = ['uniqueId', 'Value (ppm)', 'Value Error (ppm)', 
                   'SequenceCode', 'ResidueType', 'AtomName', 'Total Peak Count']
    
    # Check if all required columns exist
    missing_cols = [col for col in column_order if col not in df_grouped.columns]
    if missing_cols:
        print(f"Warning: Missing columns in output: {missing_cols}")
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
    
    print(f"\n=== Processing Summary ===")
    print(f"Original rows: {original_rows}")
    print(f"Final rows: {final_rows}")
    print(f"Duplicates removed: {duplicates_removed}")
    print(f"Transformation completed. Modified file saved as '{output_csv}'")
    
    print("\n=== First 5 rows of processed data ===")
    print(df_final.head())
    
    return df_final

def additional_processing(df, proline_map, output_file):
    """Additional processing step to prepare final output files."""
    print("\n=== Starting additional processing ===")
    
    # Create a copy for processing
    data_zero = df.copy()
    
    # Define columns to delete
    columns_to_delete = ['uniqueId', 'NmrAtom', 'ResidueType', 'Total Peak Count']
    
    # Drop specified columns if they exist
    data_columns = data_zero.drop(columns=[col for col in columns_to_delete if col in data_zero.columns])
    print(f"Dropped columns: {[col for col in columns_to_delete if col in data_zero.columns]}")
    
    # Change the numeration of proline (assuming proline is a dictionary mapping)
    if proline_map:
        data_columns['SequenceCode'] = data_columns['SequenceCode'].replace(proline_map)
        print("Applied proline mapping")
    
    # Remove all attributions with -1 (duplicates)
    initial_rows = len(data_columns)
    data_columns = data_columns[~data_columns['SequenceCode'].astype(str).str.contains('-1')]
    filtered_rows = len(data_columns)
    print(f"Removed {initial_rows - filtered_rows} rows containing '-1' in SequenceCode")
    
    # Create a new column with incremental numbers, beginning by 1
    data_columns.insert(0, 'Index', range(1, len(data_columns) + 1))
    
    # Reorder the columns - handle both cleaned and original column names
    possible_value_cols = ['Value (ppm)', 'Value\n(ppm)']
    possible_error_cols = ['Value Error (ppm)', 'Value Error\n(ppm)']
    
    value_col = next((col for col in possible_value_cols if col in data_columns.columns), None)
    error_col = next((col for col in possible_error_cols if col in data_columns.columns), None)
    
    if value_col is None or error_col is None:
        print("Warning: Could not find Value or Value Error columns")
        print("Available columns:", data_columns.columns.tolist())
        return None
    
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
    prot_file_path = os.path.join(output_file, 'name.prot')
    data_columns.to_csv(prot_file_path, sep='\t', index=False)
    print(f"Saved processed data to: {prot_file_path}")
    
    print(f"Final data shape: {data_columns.shape}")
    print("Final data preview:")
    print(data_columns.head())
    
    return data_columns


# File paths
input_csv = 'attrib.csv'
prot_seq_file = 'prot.seq'
output_csv = 'attrib_wo_double.csv'
output_directory = 'final'

# Define proline mapping (you'll need to define this based on your needs)
# Example: proline = {'P1': 'PRO1', 'P2': 'PRO2'}  # Adjust as needed
proline_map = {}  # Define your proline mapping here

# Process the data
result = process_protein_data(input_csv, prot_seq_file, output_csv)

if result is not None:
    print("\nInitial processing completed successfully!")

    # Additional processing
    final_result = additional_processing(result, proline_map, output_directory)

    if final_result is not None:
        print("\nAdditional processing completed successfully!")

        try:
            seq_file = os.path.join(output_directory, 'prot.seq')
            name_prot_file = os.path.join(output_directory, 'name.prot')
            output_cya_file = os.path.join(output_directory, 'attrib_cya.prot')

            trad = monentlature.nomenclature(seq=seq_file, prot=name_prot_file, output_file=output_directory)
            monentlature.file_transforme(trad, prot=name_prot_file, output_file=output_cya_file)
            print(f"Nomenclature conversion completed: {output_cya_file}")
        except ImportError:
            print("Warning: nomenclature module not found. Skipping final conversion.")
        except Exception as e:
            print(f"Error in nomenclature conversion: {e}")
    else:
        print("Additional processing failed!")
else:
    print("Initial processing failed!")
