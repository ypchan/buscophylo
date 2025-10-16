#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_assembly_pieces.py — Merge multiple FASTA pieces per assembly into a single FASTA.

- Input: 2-column list (assembly_id, piece_path), read from file or stdin ("-").
- Output: one merged FASTA per assembly in the specified output directory.
- Delimiter: whitespace by default (robust for typical paths). Use --delimiter '\t' for tab-separated lists.

Example:
    cat file_list.txt | python merge_assembly_pieces.py - -o outdir

Why this rewrite?
- No shell 'cat' (avoids shell injection and command-length limits).
- Streams bytes with minimal memory usage.
- Ensures a newline between pieces to prevent header gluing.
- Clear error reporting per assembly; continues processing others.
"""

import sys
import os
import argparse
import fileinput
from collections import defaultdict
from typing import Dict, List, Optional

# Optional tqdm progress bar: fall back to a no-op if not installed
try:
    from tqdm import tqdm  # type: ignore
except Exception:  # pragma: no cover
    def tqdm(x, **kwargs):  # type: ignore
        return x


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        prog="merge_assembly_pieces.py",
    )
    p.add_argument(
        "list",
        type=str,
        help="Two-column list file (assembly_id, piece_path). Use '-' to read from stdin.",
    )
    p.add_argument(
        "-o", "--outdir",
        type=str,
        default=".",
        help="Output directory. Each assembly becomes <outdir>/<assembly_id>.fna",
    )
    p.add_argument(
        "--delimiter",
        type=str,
        default=None,
        help=r"Field delimiter. Default: any whitespace. Use '\t' for tab-separated.",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files (default: skip if exists).",
    )
    p.add_argument(
        "--no-join-newline",
        action="store_true",
        help="Do NOT insert a newline between pieces (default: insert one).",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be merged without writing files.",
    )
    return p.parse_args()


def read_mapping(list_source: str, delimiter: Optional[str]) -> Dict[str, List[str]]:
    """
    Read the two-column mapping (assembly_id -> list of piece paths) while preserving order.
    Skips blank lines and lines starting with '#'.
    """
    mapping: Dict[str, List[str]] = defaultdict(list)
    with fileinput.input(files=list_source, mode="r") as infh:
        for i, line in enumerate(infh, start=1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue

            if delimiter:
                parts = s.split(delimiter, 1)
                if len(parts) != 2:
                    raise ValueError(f"Line {i}: expected 2 columns separated by {repr(delimiter)}; got: {s}")
                asm, piece = parts[0].strip(), parts[1].strip()
            else:
                # Split on any whitespace. If the path itself contains spaces, prefer using --delimiter '\t'.
                parts = s.split()
                if len(parts) < 2:
                    raise ValueError(f"Line {i}: expected at least 2 fields; got: {s}")
                asm, piece = parts[0], " ".join(parts[1:])  # join the rest as path (best-effort)

            if not asm or not piece:
                raise ValueError(f"Line {i}: empty assembly ID or piece path.")

            mapping[asm].append(piece)
    return mapping


def ensure_outdir(path: str) -> None:
    """Create output directory if needed."""
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as e:
        raise RuntimeError(f"Failed to create output directory '{path}': {e}")


def merge_one(assembly: str, pieces: List[str], outdir: str,
             join_newline: bool, overwrite: bool, dry_run: bool) -> Optional[str]:
    """
    Merge all pieces for a single assembly into <outdir>/<assembly>.fna.

    Returns:
        None on success; error message string on failure.
    """
    out_path = os.path.join(outdir, f"{assembly}.fna")

    if dry_run:
        sys.stdout.write(f"[DRY-RUN] {assembly} -> {out_path}\n")
        for p in pieces:
            sys.stdout.write(f"  - {p}\n")
        return None

    if os.path.exists(out_path) and not overwrite:
        # Skip to avoid unintended overwrite
        return None

    # Stream copy each piece into the output. Insert a newline between pieces if requested.
    try:
        with open(out_path, "wb") as outf:
            first = True
            for piece in pieces:
                if not os.path.exists(piece):
                    return f"Piece not found: {piece}"
                if not os.path.isfile(piece):
                    return f"Not a regular file: {piece}"

                # Insert a newline separator between pieces (safer for FASTA concatenation)
                if not first and join_newline:
                    outf.write(b"\n")
                first = False

                # Copy bytes in chunks; no shelling out
                with open(piece, "rb") as inf:
                    # Use a moderate buffer size; Python's default is fine, but explicit loop is clear
                    while True:
                        chunk = inf.read(1024 * 1024)  # 1 MiB
                        if not chunk:
                            break
                        outf.write(chunk)

    except OSError as e:
        return f"I/O error while writing '{out_path}': {e}"
    except Exception as e:
        return f"Unexpected error for '{assembly}': {e}"

    return None


def main() -> None:
    args = parse_args()

    # Prepare output directory
    ensure_outdir(args.outdir)

    # Build mapping assembly -> piece list
    try:
        mapping = read_mapping(args.list, args.delimiter)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read mapping: {e}")

    if not mapping:
        sys.exit("[ERROR] No assemblies found in input list.")

    # Iterate with progress bar
    errors = 0
    iterator = tqdm(mapping.items(), desc="Merging assemblies")  # type: ignore
    for assembly, pieces in iterator:
        err = merge_one(
            assembly=assembly,
            pieces=pieces,
            outdir=args.outdir,
            join_newline=not args.no_join_newline,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
        )
        if err:
            errors += 1
            sys.stderr.write(f"[WARN] {assembly}: {err}\n")

    # Exit code indicates if any assembly failed
    sys.exit(1 if (errors > 0) else 0)


if __name__ == "__main__":
    main()
