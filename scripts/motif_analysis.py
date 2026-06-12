"""Standalone: build the repeat-motif/property matrix for a FASTA file and plot a heatmap.

Usage:
    python scripts/motif_analysis.py <fasta_file> [--min_len N] [--max_len N] [--min_reps N] [--save]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.processing import process_fasta_to_matrix
from aaseq.visualization import visualize_matrix


def main():
    parser = argparse.ArgumentParser(description="Build a repeat-motif/property matrix and plot a heatmap.")
    parser.add_argument("fasta_file", help="Path to the FASTA file (.fasta)")
    parser.add_argument("--min_len", type=int, default=3, help="Minimum motif length (default: 3)")
    parser.add_argument("--max_len", type=int, default=10, help="Maximum motif length (default: 10)")
    parser.add_argument("--min_reps", type=int, default=2, help="Minimum repetition count (default: 2)")
    parser.add_argument("--workers", type=int, default=4, help="Number of CPU workers (default: 4)")
    parser.add_argument("--chunk_size", type=int, default=5000, help="Chunk size for parsing (default: 5000)")
    parser.add_argument("--save", action="store_true", help="Save the matrix to 'repeat_matrix.csv'")
    args = parser.parse_args()

    matrix = process_fasta_to_matrix(
        args.fasta_file, args.min_len, args.max_len, args.min_reps, args.workers, args.chunk_size
    )
    if matrix is None:
        return

    print("\nSequence matrix summary (Motifs & Properties):")
    print(matrix.head())

    motif_cols = [c for c in matrix.columns if not c.startswith("Prop_")]
    if motif_cols:
        visualize_matrix(matrix[motif_cols])

    if args.save:
        matrix.to_csv("repeat_matrix.csv")
        print("Saved successfully as 'repeat_matrix.csv'!")


if __name__ == "__main__":
    main()
