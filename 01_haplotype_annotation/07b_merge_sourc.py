#!/usr/bin/env python3

"""
Merge multiple files by common key in the first column.
Properly combines all the content for each key into a single line.
"""
import sys
import os

def read_file(filename):
    """Read a file and return a dictionary of key -> value mappings"""
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            parts = line.split('\t', 1)  # Split on first tab
            key = parts[0]
            value = parts[1] if len(parts) > 1 else ""
            
            # Store the full value (everything after the first tab)
            data[key] = value
    return data

def merge_files(files, output_file):
    """Merge multiple files by common key, combining data with semicolons"""
    # Dictionary to store all data
    all_data = {}
    
    # Process each file
    for file in files:
        file_data = read_file(file)
        
        # Update all_data with new info
        for key, value in file_data.items():
            if key not in all_data:
                all_data[key] = []
                
            # Add the value if not empty
            if value:
                # Split into individual attributes (assuming semicolon-separated)
                attrs = value.split(';')
                for attr in attrs:
                    if attr.strip() and attr not in all_data[key]:
                        all_data[key].append(attr.strip())
    
    # Write to output file
    with open(output_file, 'w') as out:
        for key in sorted(all_data.keys()):
            # Join all unique attributes with semicolons
            joined_value = ';'.join(all_data[key])
            out.write(f"{key}\t{joined_value}\n")
    
    return len(all_data)

def main():
    if len(sys.argv) < 4:
        print(f"Usage: {os.path.basename(sys.argv[0])} output_file file1 file2 [file3 ...]")
        print("Merges multiple files based on a common key in the first column.")
        print("Format of each file: key<tab>attributes")
        print("Attributes are assumed to be semicolon-separated.")
        return 1
    
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    
    try:
        # Verify input files exist
        for file in input_files:
            if not os.path.isfile(file):
                print(f"Error: File '{file}' not found", file=sys.stderr)
                return 1
        
        count = merge_files(input_files, output_file)
        print(f"Successfully merged {count} records to {output_file}")
        print(f"Input files: {', '.join(input_files)}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())