#!/usr/bin/python3
import os
import argparse
import pandas as pd

def create_count_matrix(directory):
    """
    Create a matrix from multiple headerless count files in a directory.
    
    Parameters:
    -----------
    directory : str
        Path to the directory containing count files
    
    Returns:
    --------
    pandas.DataFrame
        A matrix with unique IDs as index and file names as columns
    """
    # Get a list of all files in the directory
    files = [f for f in os.listdir(directory) if f.endswith('.txt')]
    
    # Check if any files found
    if not files:
        raise ValueError(f"No .txt files found in directory: {directory}")
    
    # Initialize an empty dictionary to store DataFrames
    dfs = {}
    
    # Read each file
    for file in files:
        # Remove .txt to use as column name
        column_name = os.path.splitext(file)[0].replace('_counts','')
        
        # Read the file with no header, specifying column names
        df = pd.read_csv(
            os.path.join(directory, file), 
            sep='\t',  # Tab-separated 
            header=None,  # No header in the file
            names=['GENEID', column_name]  # Assign column names
        )
        
        # Set GENEID as index
        df = df.set_index('GENEID')
        
        # Store the count column
        dfs[column_name] = df[column_name]
    # Sort column names alphanumerically
    sorted_columns = sorted(dfs.keys())
    
    # Combine all DataFrames, filling missing values with 0
    result_matrix = pd.concat([dfs[col] for col in sorted_columns], axis=1, keys=sorted_columns)
    result_matrix = result_matrix.fillna(0)
    
    return result_matrix

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='Create a count matrix from headerless count files in a directory.'
    )
    
    # Input directory argument
    parser.add_argument(
        '-i', '--input_dir', 
        type=str, 
        required=True, 
        help='Path to the directory containing input count files'
    )
    
    # Output file argument
    parser.add_argument(
        '-o', '--output', 
        type=str, 
        default='count_matrix.tsv', 
        help='Name of the output TSV file (default: count_matrix.tsv)'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Create the count matrix
    matrix = create_count_matrix(args.input_dir)
    
    # Save to TSV
    matrix.to_csv(args.output, sep="\t")
    
    print(f"Count matrix saved to {args.output}")
    print("\nMatrix Preview:")
    print(matrix.head())

if __name__ == '__main__':
    main()
