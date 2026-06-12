"""Standalone: Kyte-Doolittle hydrophobicity profile plot for a single sequence.

Usage:
    python scripts/hydrophobicity_analysis.py <fasta_file> [--record_id ID] [--window N]
    python scripts/hydrophobicity_analysis.py --seq "MKT..." [--window N]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.report import find_record
from aaseq.visualization import plot_hydrophobicity


def main():
    parser = argparse.ArgumentParser(description="Kyte-Doolittle hydrophobicity profile plot.")
    parser.add_argument("fasta_file", nargs="?", help="Path to the FASTA file (.fasta). Not required when using --seq.")
    parser.add_argument("--seq", type=str, help="Directly analyze a single amino acid sequence string (no FASTA file needed)")
    parser.add_argument("--record_id", type=str, default=None, help="Record ID to plot (default: first record)")
    parser.add_argument("--window", type=int, default=9, help="Sliding window size (default: 9)")
    args = parser.parse_args()

    if args.seq:
        sequence = args.seq.replace("\n", "").strip().upper()
        label = "Input Sequence"
    elif args.fasta_file:
        record = find_record(args.fasta_file, args.record_id)
        if record is None:
            parser.error(f"Record '{args.record_id}' not found in '{args.fasta_file}'.")
        sequence = str(record.seq)
        label = record.id
    else:
        parser.error("Either fasta_file or --seq must be provided.")

    plot_hydrophobicity(sequence, label, window=args.window)


if __name__ == "__main__":
    main()
