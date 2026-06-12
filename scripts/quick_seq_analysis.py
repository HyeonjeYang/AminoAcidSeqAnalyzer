"""Standalone: quick repeat-motif and physicochemical property analysis for a single sequence.

Usage:
    python scripts/quick_seq_analysis.py "MKTAYIAKQRQ..." [--min_len N] [--max_len N] [--min_reps N]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.motifs import getrep
from aaseq.report import print_quick_analysis


def main():
    parser = argparse.ArgumentParser(description="Quick repeat-motif and physicochemical property analysis for a single sequence.")
    parser.add_argument("seq", help="Amino acid sequence string")
    parser.add_argument("--min_len", type=int, default=3, help="Minimum motif length (default: 3)")
    parser.add_argument("--max_len", type=int, default=10, help="Maximum motif length (default: 10)")
    parser.add_argument("--min_reps", type=int, default=2, help="Minimum repetition count (default: 2)")
    args = parser.parse_args()

    sequence = args.seq.replace("\n", "").strip().upper()
    print_quick_analysis(sequence)

    repeat_counts = getrep(sequence, args.min_len, args.max_len, args.min_reps)
    print("\nRepeated motifs (min_reps reached):")
    if repeat_counts:
        for motif, count in sorted(repeat_counts.items(), key=lambda x: -x[1]):
            print(f"  {motif}: {count}")
    else:
        print("  None found.")


if __name__ == "__main__":
    main()
