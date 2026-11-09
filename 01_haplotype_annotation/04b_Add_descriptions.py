#!/usr/bin/env python3

# This script takes the output of  "_OLEA9.dmd.out2parse.tsv" files and searches for their corrsponding description in the file
# Derived from the GFF3 of the OLEA9 annotation: "OLEA9_descriptions.tsv". If a description is found, it merges it with the 
# corresponding line in the  parse file.

import sys
import re

def extract_oela_id(line):
    """Extract OELA9Hom ID from a line"""
    match = re.search(r'(OELA9Hom=[^;]+)', line)
    return match.group(1) if match else None

def extract_description(line):
    """Extract Description from a line"""
    match = re.search(r'(Description=[^;]+)', line)
    return match.group(1) if match else None

def merge_files(file1_path, file2_path, output_path):
    # Create a dictionary mapping OELA9Hom IDs to descriptions
    id_to_description = {}
    
    # Read file2 and build the mapping
    with open(file2_path, 'r') as file2:
        for line in file2:
            line = line.strip()
            if not line:
                continue
                
            oela_id = extract_oela_id(line)
            description = extract_description(line)
            
            if oela_id and description:
                id_to_description[oela_id] = description
    
    # Process file1 and add descriptions
    with open(file1_path, 'r') as file1, open(output_path, 'w') as output:
        for line in file1:
            line = line.strip()
            if not line:
                continue
                
            oela_id = extract_oela_id(line)
            
            if oela_id and oela_id in id_to_description:
                # Add the description to the line
                output.write(f"{line};{id_to_description[oela_id]}\n")
            else:
                # No matching description found, just write the original line
                output.write(f"{line}\n")

    print(f"Merged files successfully. Output written to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} file1.txt file2.txt output.txt")
        sys.exit(1)
        
    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    output_path = sys.argv[3]
    
    merge_files(file1_path, file2_path, output_path)