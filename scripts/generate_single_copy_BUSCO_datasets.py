#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BUSCO_gene_matrix_for_phylogenomic_analysis.py — Select BUSCO genes based on their taxa coverage.

- Reads and filters BUSCO results (full_table.tsv per label).
- Outputs a full matrix (URLs, Desc, counts, coverage, per-taxon status).
- Generates single-copy BUSCO FASTA datasets for genes passing coverage.

Optimized with multi-threading for faster I/O and processing.
"""

from __future__ import annotations

import os
import sys
import argparse
import fileinput
import textwrap
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Iterable

# Optional tqdm: if unavailable, silently no-op
try:
    import tqdm  # type: ignore

    def _tqdm(it, **kw):
        return tqdm.tqdm(it, **kw)
except Exception:
    def _tqdm(it, **kw):
        return it


# ----------------------------- I/O helpers -----------------------------

def read_busco_results(list_path: str) -> dict[str, str]:
    """
    Read label → full_table.tsv mapping from a file or stdin ('-').
    Each non-empty, non-comment line must contain two fields separated by whitespace:
        <label> <path/to/full_table.tsv>
    """
    mapping: dict[str, str] = {}
    for raw in fileinput.input(list_path, mode="r"):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            raise ValueError(f"Bad line (need 2 fields): {line}")
        label, full_table_path = parts[0], parts[1]
        if not os.path.exists(full_table_path):
            raise FileNotFoundError(f"full_table.tsv not found: {full_table_path} (label={label})")
        mapping[label] = full_table_path
    if not mapping:
        raise ValueError("No valid <label, full_table> pairs found.")
    return mapping


def read_busco_description(busco_desc_path: str) -> dict[str, list[str]]:
    """
    Read BUSCO gene descriptions.
    Expected TSV columns: BUSCO_ID \t Description \t OrthoDB_URL
    Returns: {busco_id: [url, desc]}
    """
    desc: dict[str, list[str]] = {}
    with open(busco_desc_path, "rt") as fh:
        for i, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(f"{busco_desc_path}:{i}: need 3 tab-separated fields (ID, Desc, URL)")
            busco_id, desc_txt, url = parts[0], parts[1], parts[2]
            desc[busco_id] = [url, desc_txt]
    if not desc:
        raise ValueError(f"No BUSCO descriptions parsed from {busco_desc_path}")
    return desc


def read_fasta(fasta_file: str) -> str:
    """
    Read a FASTA and return the concatenated sequence (headers ignored).
    """
    seq_chunks: list[str] = []
    with open(fasta_file, "rt") as infh:
        for raw in infh:
            line = raw.rstrip("\n")
            if not line or line.startswith(">"):
                continue
            seq_chunks.append(line)
    return "".join(seq_chunks)


# --------------------- Parallelized BUSCO result processing ---------------------

def _parse_full_table(full_table_path: str) -> dict[str, str]:
    """
    Parse one BUSCO full_table.tsv into {busco_id: status}.
    Columns vary with BUSCO versions; we only need first two fields:
      0: BUSCO id   1: status (Complete/Fragmented/Missing/etc.)
    Ignores comment lines starting with '#'.
    """
    res: dict[str, str] = {}
    with open(full_table_path, "rt") as fh:
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            bid, status = parts[0], parts[1]
            # First occurrence wins; later lines for same BUSCO id are typically redundant
            res.setdefault(bid, status)
    return res


def process_busco_results(label2table: dict[str, str], threads: int, show_progress: bool) -> dict[str, dict[str, str]]:
    """
    Read many full_table.tsv files in parallel.
    Returns: {label: {busco_id: status}}
    """
    out: dict[str, dict[str, str]] = {}
    items = list(label2table.items())
    iterator = _tqdm(items, desc="Reading full_table.tsv", disable=not show_progress)

    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        fut2label = {ex.submit(_parse_full_table, path): label for label, path in items}
        for fut in as_completed(fut2label):
            label = fut2label[fut]
            try:
                out[label] = fut.result()
            except Exception as e:
                # Keep going; mark label as empty so downstream can still proceed
                sys.stderr.write(f"[WARN] Failed to parse {label}: {e}\n")
                out[label] = {}
    return out


# ----------------------------- Matrix construction -----------------------------

def construct_matrix(busco_desc: dict[str, list[str]],
                     label2table: dict[str, str]) -> pd.DataFrame:
    """
    Build the scaffold DataFrame with rows=BUSCO IDs, columns=[meta + labels].
    Meta columns: OrthoDB_URL, Desc, No_taxa, Coverage%
    """
    labels = list(label2table.keys())
    colnames = ["OrthoDB_URL", "Desc", "No_taxa", "Coverage%"] + labels
    # Fix row order to all BUSCO ids in description file
    row_index = list(busco_desc.keys())

    df = pd.DataFrame(index=row_index, columns=colnames)
    # Fill URL/Desc using the same row order
    df["OrthoDB_URL"] = [busco_desc[bid][0] for bid in row_index]
    df["Desc"]        = [busco_desc[bid][1] for bid in row_index]
    return df


def add_busco_status_to_matrix(df: pd.DataFrame,
                               label2status: dict[str, dict[str, str]]) -> pd.DataFrame:
    """
    For each label (taxon), fill a column with BUSCO status per BUSCO id.
    Missing defaults to 'Missing'.
    """
    for label, status_map in label2status.items():
        # Fast vectorized map with fallback
        df[label] = df.index.to_series().map(lambda bid: status_map.get(bid, "Missing"))
    return df


# ----------------------------- Filtering & outputs -----------------------------

def filter_matrix(df: pd.DataFrame, min_coverage: int, out_path: str, show_progress: bool) -> pd.DataFrame:
    """
    Compute per-BUSCO counts / coverage, write the full matrix, and return filtered rows.
    Coverage is the % of taxa with status == 'Complete'.
    """
    # Number of taxa-columns (after 4 meta columns)
    n_taxa_cols = df.shape[1] - 4
    if n_taxa_cols <= 0:
        raise ValueError("No taxa columns found in the matrix (did you pass any labels?).")

    # Count 'Complete' across taxa columns
    df["No_taxa"] = (df.iloc[:, 4:] == "Complete").sum(axis=1)
    df["Coverage%"] = (df["No_taxa"] / n_taxa_cols * 100).round(2)

    df_filtered = df[df["Coverage%"] >= float(min_coverage)].copy()
    n_removed = df.shape[0] - df_filtered.shape[0]
    print(f"[INFO] BUSCO genes filtered out by coverage < {min_coverage}%: {n_removed}", file=sys.stdout, flush=True)

    # Write the full matrix (includes rows that were filtered below threshold)
    df.to_csv(out_path, sep="\t", index=True)
    return df_filtered


# --------------- Build single-copy BUSCO FASTA dataset (parallel) ---------------

def _load_single_copy_sequence(task: tuple[str, str, str]) -> tuple[str, str, str] | None:
    """
    Worker to load one single-copy FASTA (if exists).
    Task: (busco_id, label, fasta_path)
    Returns: (busco_id, label, sequence) or None if missing/failed.
    """
    busco_id, label, fasta_path = task
    try:
        if not os.path.exists(fasta_path):
            # It's normal if some taxa lack a given BUSCO; just skip
            return None
        seq = read_fasta(fasta_path).upper()
        if not seq:
            return None
        return (busco_id, label, seq)
    except Exception as e:
        sys.stderr.write(f"[WARN] Failed to read {fasta_path}: {e}\n")
        return None


def construct_busco_gene_dataset(label2table: dict[str, str],
                                 df_filtered: pd.DataFrame,
                                 threads: int,
                                 show_progress: bool) -> dict[str, dict[str, str]]:
    """
    Build a dict {busco_id: {label: sequence}} for BUSCO IDs passing coverage.
    Loads sequences in parallel across (busco_id, label) pairs with status == 'Complete'.
    """
    # Prepare tasks: only where status == 'Complete'
    tasks: list[tuple[str, str, str]] = []
    busco_ids = list(df_filtered.index)
    labels = [col for col in df_filtered.columns if col not in ("OrthoDB_URL", "Desc", "No_taxa", "Coverage%")]

    for label in labels:
        full_table_path = label2table[label]
        base = os.path.dirname(full_table_path)
        sc_path = os.path.join(base, "busco_sequences", "single_copy_busco_sequences")
        for busco_id in busco_ids:
            if df_filtered.at[busco_id, label] == "Complete":
                tasks.append((busco_id, label, os.path.join(sc_path, f"{busco_id}.faa")))

    dataset: dict[str, dict[str, str]] = {}
    iterator = _tqdm(tasks, desc="Loading single-copy sequences", disable=not show_progress)

    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        futs = [ex.submit(_load_single_copy_sequence, task) for task in iterator]
        for fut in as_completed(futs):
            res = fut.result()
            if res is None:
                continue
            busco_id, label, seq = res
            d = dataset.get(busco_id)
            if d is None:
                dataset[busco_id] = {label: seq}
            else:
                d[label] = seq
    return dataset


def write_busco_single_copy_dataset(single_copy: dict[str, dict[str, str]],
                                    out_dir: str,
                                    show_progress: bool) -> None:
    """
    Write FASTA files per BUSCO ID: <out_dir>/<BUSCO_ID>.faa
    Each file contains multiple sequences (one per label/taxon) wrapped at 80 columns.
    """
    os.makedirs(out_dir, exist_ok=True)
    ids = list(single_copy.keys())
    iterator = _tqdm(ids, desc="Writing single-copy FASTAs", disable=not show_progress)

    for busco_id in iterator:
        taxa2seq = single_copy[busco_id]
        out_path = os.path.join(out_dir, f"{busco_id}.faa")
        with open(out_path, "wt") as ofh:
            for taxa, seq in taxa2seq.items():
                ofh.write(f">{taxa}\n{textwrap.fill(seq, width=80)}\n")


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
        help="File or '-' with lines: <label> <path/to/run_*/full_table.tsv>"
    )
    p.add_argument(
        "-B", "--busco_desc",
        metavar="<BUSCO_gene_description.txt>",
        type=str,
        required=True,
        help="BUSCO description TSV with columns: BUSCO_ID, Desc, OrthoDB_URL"
    )
    p.add_argument(
        "-o", "--out_matrix",
        metavar="<out_matrix.tsv>",
        type=str,
        default="busco_full_matrix.tsv",
        help="Output full matrix (TSV). Default: busco_full_matrix.tsv"
    )
    p.add_argument(
        "-t", "--threads",
        metavar="<num_threads>",
        type=int,
        default=4,
        help="Threads for parallel I/O (full_table parsing and FASTA loading)."
    )
    p.add_argument(
        "-O", "--out_dir",
        metavar="<out_directory>",
        type=str,
        default="single_copy_BUSCO_dataset",
        help="Directory to write single-copy BUSCO FASTAs."
    )
    p.add_argument(
        "-p", "--progress",
        action="store_true",
        help="Show progress bars (requires tqdm)."
    )
    p.add_argument(
        "-c", "--taxa_coverage",
        metavar="<int>",
        type=int,
        choices=range(1, 101),
        default=80,
        help="Minimum taxa coverage percentage to keep a BUSCO gene. Default: 80"
    )
    return p


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    # 1) Read inputs
    label2table = read_busco_results(args.label_busco_full_table)
    busco_desc = read_busco_description(args.busco_desc)

    # 2) Build matrix scaffold
    df = construct_matrix(busco_desc, label2table)

    # 3) Parse full_table.tsv files in parallel and fill statuses
    label2status = process_busco_results(label2table, threads=args.threads, show_progress=args.progress)
    df = add_busco_status_to_matrix(df, label2status)

    # 4) Compute coverage, write full matrix, filter by threshold
    df_filtered = filter_matrix(df, args.taxa_coverage, args.out_matrix, show_progress=args.progress)

    # 5) Build single-copy dataset (only retained BUSCO ids) with parallel FASTA loading
    single_copy = construct_busco_gene_dataset(label2table, df_filtered, threads=args.threads, show_progress=args.progress)

    # 6) Write per-BUSCO FASTAs
    write_busco_single_copy_dataset(single_copy, args.out_dir, show_progress=args.progress)

    print("Done.", file=sys.stdout, flush=True)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # Allow piping like `| head` without a traceback
        try:
            sys.stdout.close()
        except Exception:
            pass
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(130)
