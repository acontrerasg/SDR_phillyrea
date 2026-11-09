#!/usr/bin/env python3
"""
Merge a GFF annotation file with key-value data file.
Uses the ID field in the GFF attributes to match with the key in the data file.
"""
import sys
import re
import os

def parse_gff_attributes(attr_string):
    """Parse GFF attribute string into a dictionary"""
    attrs = {}
    # Split by semicolons
    for pair in attr_string.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            attrs[key] = value
    return attrs

def read_data_file(filename):
    """Read key-value data file and return a dictionary"""
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip()
                
                # If the key starts with "ID=", extract just the ID part
                if key.startswith("ID="):
                    key = key[3:]
                    
                data[key] = value
    return data

def merge_gff_with_data(gff_file, data_file, output_file):
    """Merge GFF file with data file using ID field, only for 'gene' features"""
    # Read data file into memory
    supplementary_data = read_data_file(data_file)
    
    # Process GFF file and write enriched output
    with open(gff_file, 'r') as gff, open(output_file, 'w') as out:
        for line in gff:
            line = line.strip()
            if not line or line.startswith('#'):
                # Write comment lines unchanged
                out.write(line + '\n')
                continue
                
            # Parse the GFF line
            fields = line.split('\t')
            if len(fields) < 9:
                # Not enough fields for a standard GFF entry, write as is
                out.write(line + '\n')
                continue
                
            # Check if column 3 (feature type) is "gene"
            if fields[2] == "gene":
                # The 9th field contains attributes
                attr_string = fields[8]
                attrs = parse_gff_attributes(attr_string)
                
                # Look for ID field to match with supplementary data
                if 'ID' in attrs:
                    id_value = attrs['ID']
                    
                    # Check if we have supplementary data for this ID
                    if id_value in supplementary_data:
                        # Append the supplementary data to the attributes
                        fields[8] = attr_string + ';' + supplementary_data[id_value]
            
            # Write the updated line
            out.write('\t'.join(fields) + '\n')
    
    return True

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {os.path.basename(sys.argv[0])} gff_file data_file output_file")
        print("Merges a GFF annotation file with a key-value data file using the ID field.")
        print("Merging is ONLY performed for features of type 'gene' (column 3)")
        print("  gff_file: Input GFF annotation file")
        print("  data_file: Key-value data file (format: 'ID=key data')")
        print("  output_file: Output GFF file with merged data")
        return 1
    
    gff_file = sys.argv[1]
    data_file = sys.argv[2]
    output_file = sys.argv[3]
    
    try:
        # Verify input files exist
        for file in [gff_file, data_file]:
            if not os.path.isfile(file):
                print(f"Error: File '{file}' not found", file=sys.stderr)
                return 1
        
        success = merge_gff_with_data(gff_file, data_file, output_file)
        if success:
            print(f"Successfully merged GFF file with supplementary data.")
            print(f"Output written to: {output_file}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
