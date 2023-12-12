#!/usr/bin/env python3

# Author: Kristina Gagalova
# Description: combine bam stat results in one unique file

import re
import argparse
import os

class FileColumnParserRe:
    def __init__(self):
        self.data_frames = {}
    
    def get_pattern_counts(self, files, out_file):
        """
        Combine star ReadCounts output from multiple files
        """
        patterns = [
            r"mapq >= mapq_cut \(unique\):\s+(\d+)",
            r"Reads mapped in proper pairs:\s+(\d+)",
            r"Splice reads:\s+(\d+)"
            ]

        try:
            
            # Dictionary to store counts for each file
            data = {}

            # Process each file
            for file in sorted(files):
                sample_name = file.split('.')[0]  # Extracting sample name from file name (adjust according to your file naming pattern)

                with open(os.path.join(file), 'r') as file:
                    content = file.read()
                    counts = []

                for pattern in patterns:
                    match = re.search(pattern, content)
                    if match:
                        count = match.group(1)
                        counts.append(count)

                data[sample_name] = counts

            # Write the counts to a final output file
            with open(out_file, 'w') as output:
            # Write header
                output.write("Sample Name\tMapq >= Mapq_cut(unique)\tReads mapped in proper pairs\tSplice reads\n")

                # Write data
                for sample_name, counts in data.items():
                    sample_clean = re.sub("_bam-stats", "", str(sample_name))
                    new_counts = '\t'.join(counts)
                    output.write(f"{sample_clean}\t{new_counts}\n")
        
        except Exception as e:
            print(f"Error processing files: {e}")

def main():

    parser = argparse.ArgumentParser(description="Select specified columns from input files")
    parser.add_argument("--input", nargs='+', help="List of input files", required=True)
    parser.add_argument("--output", help="Name of the output file", required=True)
    args = parser.parse_args()

    selector = FileColumnParserRe()
    selector.get_pattern_counts(args.input, args.output)

if __name__ == "__main__":
    main()
