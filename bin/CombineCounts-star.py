#!/usr/bin/env python3

# Author: Kristina Gagalova
# Description: combine counts results in one unique file

import argparse
import pandas as pd

class FileColumnSelectorStar:
    def __init__(self):
        self.data_frames = {}
    
    def combine_counts_star(self, files, strandedness):
        """
        Combine star ReadCounts output from multiple files
        """

        strand = self.get_column_from_strandedness_star(strandedness)
        try:
            # Read the first file to get names and detect rows to skip
            df_first = pd.read_csv(files[0], sep='\t')
            rows_to_skip = df_first[df_first.iloc[:, 0].astype(str).str.startswith('N_')].index
            
            # Process each file and extract the nth column, excluding rows to skip
            for file in sorted(files):
                df = pd.read_csv(file, sep='\t', comment='#', skiprows=rows_to_skip)  # Using tab as the delimiter
                column_name = f"Column_{strandedness}_from_{file}"
                self.data_frames[file] = df.iloc[:, strand]
            
            # Create a DataFrame using the first column from the first file as names
            names = pd.read_csv(files[0], sep='\t', skiprows=rows_to_skip).iloc[:, 0]
            final_df = pd.DataFrame({names.name: names})
            
            # Add columns from other files using file names as column names
            for file, column_data in self.data_frames.items():
                final_df[file.replace("_ReadsPerGene.out.tab", "")] = column_data
            
            # Output file
            final_df.columns.values[0] = "Genes"
            print(final_df.to_string(index=False))

        except Exception as e:
            print(f"Error processing files: {e}")
        
    def get_column_from_strandedness_star(self, strandedness):
        """
        Assign columns to parse based on the type of library
        """
        column_n = 0
        if strandedness == 'RF':
            column_n = 2
        elif strandedness == 'FR':
            column_n = 3
        elif strandedness == 'unstranded':
            column_n = 1
        return column_n

def main():

    parser = argparse.ArgumentParser(description="Select specified columns from input files")
    parser.add_argument("--input", nargs='+', help="List of input files", required=True)
    parser.add_argument("--strandedness", choices=['RF', 'FR', 'unstranded'], default="RF", help="Library preparation stradedness (default: RF)")
    args = parser.parse_args()

    selector = FileColumnSelectorStar()
    selector.combine_counts_star(args.input, args.strandedness)

if __name__ == "__main__":
    main()

