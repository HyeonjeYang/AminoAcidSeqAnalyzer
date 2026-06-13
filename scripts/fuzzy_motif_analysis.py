"""Standalone: exact/fuzzy repeated motif discovery and shuffle enrichment.

Usage:
    python scripts/fuzzy_motif_analysis.py <fasta_file> [--record_id ID] [--fuzzy] [--stats] [--save]
    python scripts/fuzzy_motif_analysis.py --seq "MKT..." --fuzzy --stats
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.motifs import get_fuzzy_repeats, getrep
from aaseq.report import find_record
from aaseq.statistics import motif_enrichment


def main():
    parser = argparse.ArgumentParser(description="Exact/fuzzy repeated motif analysis.")
    parser.add_argument("fasta_file", nargs="?", help="Path to FASTA file. Not required with --seq.")
    parser.add_argument("--seq", type=str, help="Directly analyze one amino acid sequence")
    parser.add_argument("--record_id", type=str, default=None, help="Record ID to analyze (default: first record)")
    parser.add_argument("--min_len", type=int, default=3)
    parser.add_argument("--max_len", type=int, default=10)
    parser.add_argument("--min_reps", type=int, default=2)
    parser.add_argument("--fuzzy", action="store_true", help="Run fuzzy motif family detection")
    parser.add_argument("--fuzzy_mismatches", type=int, default=1)
    parser.add_argument("--stats", action="store_true", help="Run shuffle-based motif enrichment")
    parser.add_argument("--num_shuffles", type=int, default=200)
    parser.add_argument("--seed", type=int, default=None)
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

    exact = getrep(sequence, args.min_len, args.max_len, args.min_reps)
    print("\nExact repeated motifs:")
    print(exact if exact else "  None")

    if args.fuzzy:
        fuzzy = get_fuzzy_repeats(sequence, args.min_len, args.max_len, args.min_reps, args.fuzzy_mismatches)
        print("\nFuzzy repeated motif families:")
        print(fuzzy if fuzzy else "  None")

    if args.stats:
        stats = motif_enrichment(
            sequence,
            args.min_len,
            args.max_len,
            args.min_reps,
            args.num_shuffles,
            seed=args.seed,
            observed_counts=exact,
        )
        print("\nShuffle-based motif enrichment:")
        print(stats.head(20))
        if args.save:
            stats.to_csv(f"{label}_motif_enrichment.csv", index=False)


if __name__ == "__main__":
    main()
