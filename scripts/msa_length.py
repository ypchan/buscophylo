#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
msa_length.py — Fast threaded alignment-length scanner for FASTA or gzipped FASTA files.

Features:
- Reads only the first sequence in each MSA (until the next '>' header).
- Constant memory footprint (streaming line-by-line).
- Supports .gz files automatically.
- Optional: ignore gap characters ('-') with -u / --ungap.
- Uses native Python threads (-t / --threads).
- Robust: catches file errors individually, continues on others.
"""

import sys
import os
import argparse
import gzip
import threading
from queue import Queue
from typing import Optional, Tuple, Dict, List


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description="Compute MSA alignment length from the first sequence (supports .gz).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog="msa_length.py",
    )
    p.add_argument("msa", nargs="+", help="One or more MSA files (FASTA or gzipped).")
    p.add_argument(
        "-u", "--ungap", action="store_true",
        help="Ignore gaps ('-') and count only non-gap residues."
    )
    p.add_argument(
        "-t", "--threads", type=int, default=min(8, os.cpu_count() or 4),
        help="Number of threads for concurrent file processing."
    )
    p.add_argument(
        "-s", "--sep", type=str, default="\t",
        help="Field separator for output table."
    )
    p.add_argument(
        "-S", "--strict-first", action="store_true",
        help="Stop after the first sequence line (if each record is single-line)."
    )
    return p.parse_args()


def open_maybe_gzip(path: str):
    """Open a text file, automatically handling gzip if needed."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def scan_alignment_length(path: str, ungap: bool, strict_first: bool) -> Tuple[Optional[int], Optional[str]]:
    """
    Scan a FASTA file and compute the length of its first sequence.
    Stops reading at the next '>' header.
    Returns: (length, error_message)
    """
    length = 0
    in_first_seq = False

    try:
        with open_maybe_gzip(path) as fh:
            for raw in fh:
                line = raw.rstrip("\n\r")

                # Header line
                if line.startswith(">"):
                    if in_first_seq:
                        # Encountered the second header → stop
                        break
                    in_first_seq = True
                    continue

                # Sequence lines
                if in_first_seq:
                    if not line:
                        if strict_first:
                            # Stop immediately if single-line FASTA mode
                            break
                        continue

                    # Remove spaces or tabs (rare in FASTA)
                    if " " in line or "\t" in line:
                        line = "".join(ch for ch in line if ch not in (" ", "\t"))

                    # Count characters
                    if ungap:
                        length += sum(1 for ch in line if ch != "-")
                    else:
                        length += len(line)

                    if strict_first:
                        break

        if not in_first_seq:
            return None, "No FASTA header ('>') found."
        return length, None

    # Error handling
    except FileNotFoundError:
        return None, "File not found."
    except PermissionError:
        return None, "Permission denied."
    except UnicodeDecodeError:
        return None, "Decode error (not a text FASTA?)."
    except OSError as e:
        return None, f"I/O error: {e}"
    except Exception as e:
        return None, f"Unexpected error: {e}"


def worker(q: Queue, results: Dict[str, Tuple[Optional[int], Optional[str]]],
           lock: threading.Lock, ungap: bool, strict_first: bool):
    """Worker thread: consume file paths from queue, compute alignment length, save results."""
    while True:
        path = q.get()
        if path is None:
            q.task_done()
            break
        length, err = scan_alignment_length(path, ungap, strict_first)
        with lock:
            results[path] = (length, err)
        q.task_done()


def main():
    args = parse_args()

    # Deduplicate file list while preserving order
    seen = set()
    files: List[str] = []
    for f in args.msa:
        if f not in seen:
            seen.add(f)
            files.append(f)

    # Thread-safe containers
    results: Dict[str, Tuple[Optional[int], Optional[str]]] = {}
    lock = threading.Lock()
    q: Queue = Queue()

    # Start threads
    n_threads = max(1, args.threads)
    threads: List[threading.Thread] = []
    for _ in range(n_threads):
        th = threading.Thread(
            target=worker,
            args=(q, results, lock, args.ungap, args.strict_first),
            daemon=True,
        )
        th.start()
        threads.append(th)

    # Enqueue tasks
    for f in files:
        q.put(f)

    # Send termination signals
    for _ in range(n_threads):
        q.put(None)

    # Wait until all tasks complete
    q.join()

    # Output results
    sep = args.sep
    out = sys.stdout
    err = sys.stderr

    print(sep.join(["file", "length"]), file=out)
    had_error = False
    for f in files:
        length, e = results.get(f, (None, "No result (threading error)."))
        if e is None and length is not None:
            print(f"{f}{sep}{length}", file=out)
        else:
            had_error = True
            print(f"[WARN] {f}: {e}", file=err)

    sys.exit(1 if had_error else 0)


if __name__ == "__main__":
    main()
