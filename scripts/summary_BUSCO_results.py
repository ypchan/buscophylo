#!/usr/bin/env python3
'''
summary_BUSCO_results.py -- Summary multiple BUSCO result files into a table.

Date:
    2020-04-14
Bugs：
    Any bugs should be reported to yanpengch@qq.com

Usage:
    find . -name 'short_summary.specific*txt' | summary_BUSCO_results.py - > BUSCO_results_summary.txt
'''

import os
import sys
import re
import fileinput
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_busco_result2lst(busco_result_file: str) -> list:
    """
    Parse the BUSCO result file (short_summary.specific*).
    Extracts key information such as completeness, single-copy, and duplicated BUSCOs.
    Returns a list of values to represent the summary for that genome.
    """
    with open(busco_result_file, 'r') as buscofh:
        genome_label = ""
        per_C, per_S, per_D, per_F, per_M, num_Total = None, None, None, None, None, None
        num_C, num_S, num_D, num_F, num_M = None, None, None, None, None

        for line in buscofh:
            line = line.strip()
            
            if line.startswith('# Summarized benchmarking in BUSCO notation for file'):
                genome_label = os.path.basename(line.split(' ')[-1])

            # Extract completion percentages
            if line.startswith('C:'):
                # Extract relevant percentages and numbers using regular expressions
                index_lst = [1, 3, 5, 8, 10, 12]
                per_C, per_S, per_D, per_F, per_M, num_Total = [re.split(r'[:\[\],]', line)[pos] for pos in index_lst]

            # Extract the specific counts for each BUSCO category
            if 'Complete BUSCOs (C)' in line:
                num_C = line.split()[0]
            if 'Complete and single-copy BUSCOs (S)' in line:
                num_S = line.split()[0]
            if 'Complete and duplicated BUSCOs (D)' in line:
                num_D = line.split()[0]
            if 'Fragmented BUSCOs (F)' in line:
                num_F = line.split()[0]
            if 'Missing BUSCOs (M)' in line:
                num_M = line.split()[0]

    return [genome_label, per_C, per_S, per_D, per_F, per_M, num_C, num_S, num_D, num_F, num_M, num_Total]


def process_busco_files(files: list) -> list:
    """
    Process multiple BUSCO result files concurrently using ThreadPoolExecutor.
    Returns a list of lists where each inner list represents a parsed BUSCO result.
    """
    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(parse_busco_result2lst, file): file for file in files}
        results = []
        for future in as_completed(futures):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                sys.stderr.write(f"Error processing file {futures[future]}: {e}\n")
        return results


def write_summary_table(results: list, output_file: str):
    """
    Write the results into a tab-delimited summary table.
    """
    head_lst = [
        'genome_label', 'C', 'S', 'D', 'F', 'M', 
        'Complete BUSCOs (C)', 'Complete and single-copy BUSCOs (S)', 
        'Complete and duplicated BUSCOs (D)', 'Fragmented BUSCOs (F)', 
        'Missing BUSCOs (M)', 'Total BUSCO groups searched'
    ]
    
    with open(output_file, 'w') as ofh:
        # Write header
        ofh.write('\t'.join(head_lst) + '\n')
        
        # Write rows
        for row in results:
            ofh.write('\t'.join(row) + '\n')


def main():
    # Argument parsing
    parser = argparse.ArgumentParser(
        description="Summary of multiple BUSCO result files into a single table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        prog="summary_BUSCO_results.py"
    )
    
    parser.add_argument(
        'label_busco_full_table',
        metavar='<label_full_table_path.txt>',
        type=str,
        help='File containing labels and paths to BUSCO result files (e.g., short_summary.specific*txt).'
    )
    parser.add_argument(
        '-o', '--out_matrix',
        metavar='<out_matrix.tsv>',
        type=str,
        default='BUSCO_results_summary.txt',
        help='Output file for the summarized BUSCO results table.'
    )
    args = parser.parse_args()

    # Read BUSCO result files from the input list
    busco_result_files = [line.strip() for line in fileinput.input(files=args.label_busco_full_table)]
    
    if not busco_result_files:
        print("Error: No BUSCO result files found.", file=sys.stderr)
        sys.exit(1)

    # Process BUSCO files in parallel
    results = process_busco_files(busco_result_files)

    # Write the summary table to the output file
    write_summary_table(results, args.out_matrix)

    print("Summary table written to:", args.out_matrix)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Process interrupted by the user.", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)
