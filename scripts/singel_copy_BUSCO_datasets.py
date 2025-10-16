#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BUSCO_gene_matrix_for_phylogenomic_analysis.py — Select BUSCO genes based on their taxa coverage.

- Reads and filters BUSCO results.
- Outputs a filtered matrix of BUSCO genes and taxa coverage.
- Also generates a FASTA dataset of single-copy BUSCO genes.

Optimized with multi-threading for faster processing.
"""

import os
import sys
import argparse
import fileinput
import textwrap
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# ----------------------------- I/O helpers -----------------------------

def read_busco_results(args_busco_full_table: str) -> Dict[str, str]:
    """
    Read the list of BUSCO results (paths to full_table.tsv).
    Returns a dictionary with {label: full_table_path}.
    """
    label_busco_result_dict = {}
    for line in fileinput.input(args_busco_full_table):
        label, full_table_path = line.strip('\n').split()
        label_busco_result_dict[label] = full_table_path
    return label_busco_result_dict


def read_busco_description(args_busco_desc: str) -> Dict[str, List[str]]:
    """
    Read the BUSCO gene descriptions.
    Returns a dictionary {busco_id: [url, description]}.
    """
    busco_desc_dict = {}
    with open(args_busco_desc, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            busco_id, desc, url = line.split('\t')
            busco_desc_dict[busco_id] = [url, desc]
    return busco_desc_dict


def read_fasta(fasta_file: str) -> str:
    """
    Read a FASTA file and return the concatenated sequence as a string.
    """
    seq_lst = []
    with open(fasta_file, 'rt') as infh:
        for line in infh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                continue
            seq_lst.append(line)
    return ''.join(seq_lst)


# ----------------------------- Parallelized BUSCO Result Processing -----------------------------

def process_busco_file(label: str, full_table_path: str) -> Dict[str, str]:
    """
    Process each BUSCO full_table file to extract gene statuses (Complete/Partial/Missing).
    Returns a dictionary {busco_id: status}.
    """
    buscoid_status = {}
    with open(full_table_path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#'):
                continue
            busco_id, status = line.split('\t')[:2]
            buscoid_status[busco_id] = status
    return label, buscoid_status


def process_busco_results(label_busco_result_dict: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    """
    Process all BUSCO results in parallel using threading.
    Returns a dictionary {label: {busco_id: status}}.
    """
    busco_status_dict = {}
    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(process_busco_file, label, path): label for label, path in label_busco_result_dict.items()}
        for future in as_completed(futures):
            label, busco_status = future.result()
            busco_status_dict[label] = busco_status
    return busco_status_dict


# ----------------------------- Constructing the Matrix -----------------------------

def construct_matrix(busco_desc_dict: Dict[str, List[str]], label_busco_result_dict: Dict[str, str]) -> pd.DataFrame:
    """
    Construct the BUSCO gene matrix.
    The index represents BUSCO gene IDs, columns are taxa (labels), and we also add URL and description.
    """
    taxa_label_lst = list(label_busco_result_dict.keys())
    col_name_lst = ['OrthoDB_URL', 'Desc', 'No_taxa', 'Coverage%'] + taxa_label_lst
    row_name_lst = list(busco_desc_dict.keys())
    
    # Initialize DataFrame
    df = pd.DataFrame(columns=col_name_lst, index=row_name_lst)

    # Fill URL and Description columns
    df['OrthoDB_URL'] = [lst[0] for lst in busco_desc_dict.values()]
    df['Desc'] = [lst[1] for lst in busco_desc_dict.values()]

    return df


def add_busco_status_to_matrix(df: pd.DataFrame, busco_status_dict: Dict[str, Dict[str, str]], label_busco_result_dict: Dict[str, str]) -> pd.DataFrame:
    """
    Add the BUSCO gene status ('Complete', 'Partial', 'Missing') to the DataFrame.
    """
    for label, buscoid_status in busco_status_dict.items():
        df[label] = df.index.map(lambda busco_id: buscoid_status.get(busco_id, 'Missing'))
    return df


# ----------------------------- Filtering and Output -----------------------------

def filter_matrix(df: pd.DataFrame, taxa_coverage: int, out_matrix: str) -> pd.DataFrame:
    """
    Filter the matrix based on taxa coverage.
    Returns a filtered DataFrame and writes the full matrix to a file.
    """
    num_taxa = df.shape[1] - 4
    df['No_taxa'] = (df.iloc[:, 4:] == 'Complete').sum(axis=1)
    coverage_lst = df['No_taxa'] / num_taxa * 100
    df['Coverage%'] = [round(coverage, 2) for coverage in coverage_lst]

    # Filter based on taxa coverage threshold
    df_filtered = df[df['Coverage%'] >= taxa_coverage]
    num_busco_filtered = df.shape[0] - df_filtered.shape[0]
    print(f'[INFO] Number of BUSCO genes with low taxa coverage: {num_busco_filtered}\n', file=sys.stdout, flush=True)

    # Save the matrix to a file
    df.to_csv(out_matrix, sep='\t', index=True)
    
    return df_filtered


def construct_busco_gene_dataset(label_busco_result_dict: Dict[str, str], df_filtered: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    """
    Construct the dataset of single-copy BUSCO sequences.
    Returns a dictionary {busco_id: {taxon_id: sequence}}.
    """
    single_copy_busco_dict = {}

    busco_id_lst = df_filtered.index.tolist()
    pbar = tqdm.tqdm(label_busco_result_dict.keys(), desc="Processing BUSCO sequences")

    for label in pbar:
        full_table_path = label_busco_result_dict[label]
        single_copy_sequence_path = os.path.dirname(full_table_path) + '/busco_sequences/single_copy_busco_sequences'

        for busco_id in busco_id_lst:
            if df_filtered.loc[busco_id, label] == 'Complete':
                busco_sequence_path = f"{single_copy_sequence_path}/{busco_id}.faa"
                busco_sequence = read_fasta(busco_sequence_path)
                if busco_id not in single_copy_busco_dict:
                    single_copy_busco_dict[busco_id] = {}
                single_copy_busco_dict[busco_id][label] = busco_sequence.upper()

    return single_copy_busco_dict


def write_busco_single_copy_dataset(single_copy_busco_dict: Dict[str, Dict[str, str]], out_dir: str) -> None:
    """
    Write single-copy BUSCO sequences to FASTA files.
    """
    os.makedirs(out_dir, exist_ok=True)

    pbar = tqdm.tqdm(single_copy_busco_dict.keys(), desc="Writing single-copy BUSCO sequences")
    for busco_id in pbar:
        taxa_seq_dict = single_copy_busco_dict[busco_id]
        with open(f"{out_dir}/{busco_id}.faa", 'wt') as ofh:
            for taxa, seq in taxa_seq_dict.items():
                wrapped_seq = textwrap.fill(seq, width=80)
                ofh.write(f">{taxa}\n{wrapped_seq}\n")


# ----------------------------- CLI -----------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        prog="BUSCO_gene_matrix_for_phylogenomic_analysis.py",
    )
    p.add_argument(
        "label_busco_full_table",
        metavar="<label_full_table_path.txt>",
        type=str,
        help="A file with label and paths to BUSCO results (full_table.tsv)."
    )
    p.add_argument(
        "-B", "--busco_desc",
        metavar="<BUSCO_gene_description.txt>",
        type=str,
        required=True,
        help="BUSCO gene description file (ID, OrthoDB URL, description)."
    )
    p.add_argument(
        "-o", "--out_matrix",
        metavar="<out_matrix.tsv>",
        type=str,
        default="busco_full_matrix.tsv",
        help="Output matrix filename (default: busco_full_matrix.tsv)."
    )
    p.add_argument(
        "-t", "--threads",
        metavar="<num_threads>",
        type=int,
        default=4,
        help="Number of threads to use for parallel reading of alignments (default: 4)."
    )
    p.add_argument(
        "-O", "--out_dir",
        metavar="<out_directory>",
        type=str,
        default="single_copy_BUSCO_dataset",
        help="Output directory for single-copy BUSCO sequences."
    )
    p.add_argument(
        "-p", "--progress",
        action="store_true",
        help="Show a tqdm progress bar during processing."
    )
    p.add_argument(
        "-c", "--taxa_coverage",
        metavar="<int>",
        type=int,
        choices=range(1, 101),
        default=80,
        help="Minimum taxa coverage (default: 80%)."
    )
    return p


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    # Phase 1: Read all BUSCO results (with parallelism)
    label_busco_result_dict = read_busco_results(args.label_busco_full_table)
    busco_desc_dict = read_busco_description(args.busco_desc)

    # Phase 2: Construct the initial matrix
    df = construct_matrix(busco_desc_dict, label_busco_result_dict)

    # Phase 3: Process BUSCO statuses in parallel
    busco_status_dict = process_busco_results(label_busco_result_dict)
    df = add_busco_status_to_matrix(df, busco_status_dict, label_busco_result_dict)

    # Phase 4: Filter matrix based on taxa coverage
    df_filtered = filter_matrix(df, args.taxa_coverage, args.out_matrix)

    # Phase 5: Generate single-copy BUSCO sequences
    single_copy_busco_dict = construct_busco_gene_dataset(label_busco_result_dict, df_filtered)
    write_busco_single_copy_dataset(single_copy_busco_dict, args.out_dir)

    print("Done.", file=sys.stdout, flush=True)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # Allow piping to tools like `head` without stack traces
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(130)
