"""Standalone: IDR/LLPS-related sequence feature analysis.

This is a heuristic candidate-prioritization tool, not a trained LLPS predictor.
It uses sequence features commonly discussed for IDR-mediated condensates:
low-complexity/prion-like composition, aromatic/cationic stickers, RGG motifs,
and charge-pattern descriptors.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.llps import analyze_llps
from aaseq.report import find_record
from aaseq.visualization import plot_llps_profile


def main():
    parser = argparse.ArgumentParser(description="IDR/LLPS-related feature analysis.")
    parser.add_argument("fasta_file", nargs="?", help="Path to FASTA file. Not required with --seq.")
    parser.add_argument("--seq", type=str, help="Directly analyze one amino acid sequence")
    parser.add_argument("--record_id", default=None)
    parser.add_argument("--window", type=int, default=31)
    parser.add_argument("--threshold", type=float, default=0.6)
    parser.add_argument("--min_region", type=int, default=20)
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

    summary, profile, regions = analyze_llps(
        sequence,
        window=args.window,
        threshold=args.threshold,
        min_region=args.min_region,
    )
    print("\nLLPS / condensate feature summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    print("\nCandidate regions:")
    print(regions if not regions.empty else "  None")

    if args.save:
        profile.to_csv(f"{label}_llps_profile.csv", index=False)
        regions.to_csv(f"{label}_llps_regions.csv", index=False)
    plot_llps_profile(profile, label, regions)


if __name__ == "__main__":
    main()
