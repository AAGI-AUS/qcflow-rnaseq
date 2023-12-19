#!/usr/bin/env python3

# Author: Kristina Gagalova
# Description: combine contamination hits

import argparse
import os
import pandas as pd

class FileColumnSelectorContaminationRate:
    def __init__(self):
        self.data_frames = {}
    
    def combine_contamination_rate(self, files, out_file):
        """
        Combine BBT contamination rate from multiple files
        """

        try:
            # Read the first file to get names and detect rows to skip
            patterns = '^multiMatch|^noMatch'
            df_first = pd.read_csv(files[0], sep='\t')
            df_first = df_first[~df_first['filter_id'].str.contains(patterns)]
            # Process each file and extract the nth column, excluding rows to skip
            for file in sorted(files):
                df = pd.read_csv(file, sep='\t')  # Using tab as the delimiter
                df = df[~df['filter_id'].str.contains(patterns)]
                column_name = f"Column_counts_from_{file}"
                self.data_frames[file] = df.loc[:, "rate_hit"]
            # Create a DataFrame using the first column from the first file as names
            names = df_first.iloc[:, 0]

            final_df = pd.DataFrame({names.name: names})

            # Add columns from other files using file names as column names
            for file, column_data in self.data_frames.items():
                file_base = os.path.basename(file)
                final_df[file_base.replace("_summary.tsv", "")] = column_data

            final_df.columns.values[0] = "BBT filters"
            final_df.to_csv(out_file, index=False, sep='\t')

        except Exception as e:
            print(f"Error processing files: {e}")

def main():

    parser = argparse.ArgumentParser(description="Select rate hit from BBT input files.")
    parser.add_argument("--input", nargs='+', help="List of input files", required=True)
    parser.add_argument("--output", help="Name of the output file", required=True)
    args = parser.parse_args()

    selector = FileColumnSelectorContaminationRate()
    selector.combine_contamination_rate(args.input, args.output)

if __name__ == "__main__":
    main()
