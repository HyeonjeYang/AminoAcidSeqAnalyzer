"""Standalone: intrinsic-disorder and low-complexity analysis.

Usage:
    python scripts/idr_analysis.py <fasta_file> [--record_id ID] [--save]
    python scripts/idr_analysis.py --seq "MKT..." --save
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.disorder import analyze_idr
from aaseq.report import find_record
from aaseq.visualization import plot_idr_profile


def main():
    parser = argparse.ArgumentParser(description="IDR and low-complexity analysis.")
    parser.add_argument("fasta_file", nargs="?", help="Path to FASTA file. Not required with --seq.")
    parser.add_argument("--seq", type=str, help="Directly analyze one amino acid sequence")
    parser.add_argument("--record_id", type=str, default=None)
    parser.add_argument("--disorder_window", type=int, default=15)
    parser.add_argument("--complexity_window", type=int, default=12)
    parser.add_argument("--disorder_threshold", type=float, default=0.2)
    parser.add_argument("--low_complexity_threshold", type=float, default=0.55)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()

    if args.seq:
        sequence = args.seq.replace("\n", "").strip().upper()
        label = "Input_Sequence"
    elif args.fasta_file:
        record = find_record(args.fasta_file, args.record_id)
        if record is None:
            parser.error(f"Record '{args.record_id}' not found in '{args.fasta_file}'.")
        sequence = str(record.seq)
        label = record.id
    else:
        parser.error("Either fasta_file or --seq must be provided.")

    scores, regions, summary = analyze_idr(
        sequence,
        args.disorder_window,
        args.complexity_window,
        args.disorder_threshold,
        args.low_complexity_threshold,
    )
    print("\nIDR / low-complexity summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    print("\nRegions:")
    print(regions if not regions.empty else "  None")
    if args.save:
        scores.to_csv(f"{label}_idr_scores.csv", index=False)
        regions.to_csv(f"{label}_idr_regions.csv", index=False)
    plot_idr_profile(scores, label, regions)


if __name__ == "__main__":
    main()
