#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
concatenate_BUSCO_gene_msa.py — Concatenate trimmed BUSCO gene alignments into a supermatrix
and output a 0/1 matrix indicating the presence/absence of genes per taxon.

- Reads alignment list (MSA paths) from stdin or file.
- Outputs a FASTA supermatrix and a 0/1 presence-absence matrix (in CSV format).
- Optimized with multi-threading to read large MSA files concurrently.
"""

import os
import sys
import argparse
import fileinput
import gzip
from typing import Dict, List, Tuple, Set, Optional, Generator
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue

# Optional tqdm: silently disable if not available
try:
    import tqdm  # type: ignore
    def _tqdm(iterable, **kwargs):
        return tqdm.tqdm(iterable, **kwargs)
except Exception:  # pragma: no cover
    def _tqdm(iterable, **kwargs):
        return iterable


# ----------------------------- I/O helpers -----------------------------

def open_maybe_gzip(path: str):
    """Open a text file, transparently handling gzip files by extension (.gz)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def fasta_records(path: str) -> Generator[Tuple[str, str], None, None]:
    """
    Parse a FASTA file (possibly gzipped) and yield (id, seq) tuples.
    Sequences are concatenated across multiple lines; whitespace in sequence lines is removed.
    """
    with open_maybe_gzip(path) as fh:
        cur_id: Optional[str] = None
        chunks: List[str] = []
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Flush previous record
                if cur_id is not None:
                    yield cur_id, "".join(chunks)
                cur_id = line[1:].strip()
                chunks.clear()
            else:
                # Remove spaces/tabs just in case
                if " " in line or "\t" in line:
                    line = "".join(ch for ch in line if ch not in (" ", "\t"))
                chunks.append(line)
        # Flush last
        if cur_id is not None:
            yield cur_id, "".join(chunks)


def iter_alignment_paths(list_source: str) -> Iterable[str]:
    """
    Iterate alignment file paths from a list file (or '-' for stdin).
    Skips blank lines and lines starting with '#'.
    """
    for raw in fileinput.input(files=list_source, mode="r"):
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        yield s


# ----------------------------- Core logic -----------------------------

def read_alignment(path: str) -> Dict[str, str]:
    """
    Read a single alignment file and return a dictionary {taxon_id: sequence}.
    """
    taxa_to_seq: Dict[str, str] = {}
    for tid, seq in fasta_records(path):
        if tid in taxa_to_seq:
            raise ValueError(f"Duplicate taxon ID in alignment '{path}': {tid}")
        taxa_to_seq[tid] = seq
    return taxa_to_seq


def read_all_alignments(alignment_list_path: str, show_progress: bool = False, num_threads: int = 4) -> Dict[str, Dict[str, str]]:
    """
    Read all alignments into a nested dict:
        all_alignment_dict[alignment_basename] = { taxon_id: aligned_sequence, ... }
    """
    all_alignment_dict: Dict[str, Dict[str, str]] = {}

    paths = list(iter_alignment_paths(alignment_list_path))
    iterator = _tqdm(paths, desc="Reading alignments", disable=not show_progress)

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = {executor.submit(read_alignment, path): path for path in paths}
        for future in as_completed(futures):
            path = futures[future]
            try:
                alignment_data = future.result()
                aln_name = os.path.basename(path).split('.')[0]
                all_alignment_dict[aln_name] = alignment_data
            except Exception as e:
                sys.stderr.write(f"Error reading {path}: {e}\n")

    return all_alignment_dict


def alignment_length_and_taxa_list(all_alignment_dict: Dict[str, Dict[str, str]]) -> Tuple[Dict[str, int], List[str]]:
    """
    Compute alignment length per locus and global taxa list (unique).
    Returns:
        align_len: { alignment_basename: alignment_length }
        taxa_list: list of all unique taxon IDs (unsorted)
    """
    align_len: Dict[str, int] = {}
    taxa_set: Set[str] = set()

    for aln_name, fa_dict in all_alignment_dict.items():
        # Take length from any sequence (validated equal)
        first_len = next((len(seq) for seq in fa_dict.values()), 0)
        align_len[aln_name] = first_len
        taxa_set.update(fa_dict.keys())

    return align_len, list(taxa_set)


def placeholder_for_missing_taxa(align_len: Dict[str, int], taxa_list: List[str], all_alignment_dict: Dict[str, Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    """
    Ensure every alignment dictionary has entries for all taxa.
    Missing taxa get a gap-only placeholder ('-') of the alignment length.
    """
    # Pre-build gap strings per alignment to avoid repeated allocations
    gap_cache: Dict[str, str] = {aln: "-" * L for aln, L in align_len.items()}

    for aln_name, fa_dict in all_alignment_dict.items():
        needed = set(taxa_list) - fa_dict.keys()
        if not needed:
            continue
        gap_seq = gap_cache[aln_name]
        # Fill in missing taxa
        for taxon in needed:
            fa_dict[taxon] = gap_seq

    return all_alignment_dict


def concatenate_alignment(all_alignment_dict: Dict[str, Dict[str, str]], taxa_list: List[str], sort_taxa: bool = False) -> Dict[str, str]:
    """
    Concatenate all alignments per taxon ID.
    Returns:
        concatenated_dict: { taxon_id: concatenated_sequence }
    """
    if sort_taxa:
        taxa_iter = sorted(taxa_list)
    else:
        taxa_iter = taxa_list  # preserve discovery order

    # Order of loci in concatenation follows insertion order in dict (Py3.7+ preserves)
    loci_iter = list(all_alignment_dict.keys())

    concatenated: Dict[str, List[str]] = {t: [] for t in taxa_iter}

    # Append locus blocks per taxon
    for locus in loci_iter:
        fa_dict = all_alignment_dict[locus]
        for t in taxa_iter:
            concatenated[t].append(fa_dict[t])

    # Join to final strings (single pass; avoids repeated string concatenation)
    return {t: "".join(chunks) for t, chunks in concatenated.items()}


def write_supermatrix(out_path: str, concatenated: Dict[str, str], sort_taxa: bool = False) -> None:
    """
    Write supermatrix to FASTA.
    """
    taxa = sorted(concatenated.keys()) if sort_taxa else concatenated.keys()
    with open(out_path, "wt") as ofh:
        for t in taxa:
            ofh.write(f">{t}\n{concatenated[t]}\n")


def write_presence_absence_matrix(out_path: str, all_alignment_dict: Dict[str, Dict[str, str]], taxa_list: List[str]) -> None:
    """
    Write a 0/1 matrix showing presence (1) or absence (0) of each gene (alignment) for each taxon.
    """
    with open(out_path, "wt") as ofh:
        # Header: genes (loci)
        ofh.write("taxon_id," + ",".join(all_alignment_dict.keys()) + "\n")
        for taxon in taxa_list:
            row = [taxon]
            for aln_name in all_alignment_dict.keys():
                # Mark '1' if sequence is non-gap (else '0')
                sequence = all_alignment_dict[aln_name].get(taxon, '-' * len(next(iter(all_alignment_dict[aln_name].values()))))
                row.append("1" if sequence != "-" * len(sequence) else "0")
            ofh.write(",".join(row) + "\n")


# ----------------------------- CLI -----------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        prog="concatenate_BUSCO_gene_msa.py",
    )
    p.add_argument(
        "alignment_list",
        metavar="<trimmed_alignment.list>",
        type=str,
        help="A file with one alignment path per line (use '-' to read from stdin)."
    )
    p.add_argument(
        "-o", "--out",
        metavar="<super_BUSCO_matrix.faa>",
        type=str,
        default="super_BUSCO_matrix.faa",
        help="Output FASTA (supermatrix)."
    )
    p.add_argument(
        "-t", "--threads",
        metavar="<num_threads>",
        type=int,
        default=4,
        help="Number of threads to use for parallel reading of alignments."
    )
    p.add_argument(
        "--progress",
        action="store_true",
        help="Show a tqdm progress bar while reading alignments."
    )
    p.add_argument(
        "--sort-taxa",
        action="store_true",
        help="Sort taxa IDs alphabetically in output (deterministic order)."
    )
    p.add_argument(
        "--presence-absence",
        type=str,
        help="Write a 0/1 matrix indicating the presence/absence of each gene per taxon."
    )
    return p


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    # Phase 1: read all alignments (with validation)
    all_alignment_dict = read_all_alignments(args.alignment_list, show_progress=args.progress, num_threads=args.threads)

    # Phase 2: gather alignment lengths and global taxa
    align_len, taxa_list = alignment_length_and_taxa_list(all_alignment_dict)

    # Phase 3: fill gaps for missing taxa (per alignment)
    all_alignment_dict = placeholder_for_missing_taxa(align_len, taxa_list, all_alignment_dict)

    # Phase 4: concatenate per taxon
    concatenated = concatenate_alignment(all_alignment_dict, taxa_list, sort_taxa=args.sort_taxa)

    # Phase 5: write supermatrix to FASTA
    write_supermatrix(args.out, concatenated, sort_taxa=args.sort_taxa)

    # Phase 6: write 0/1 presence-absence matrix if requested
    if args.presence_absence:
        write_presence_absence_matrix(args.presence_absence, all_alignment_dict, taxa_list)

    print("Done", file=sys.stdout, flush=True)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # Allow piping to tools like `head` without stack traces
        try:
            sys.stdout.close()
        except Exception:
            pass
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(130)
