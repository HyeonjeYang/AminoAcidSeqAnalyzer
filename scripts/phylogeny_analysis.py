"""Standalone: phylogenetic analysis of sequences in a FASTA file.

Two independent modes (pick one or both):
  --feature  : feature-based dendrogram (clusters by motif/property similarity)
  --sequence : alignment-based dendrogram (real BLOSUM62 pairwise alignment distances)
If neither flag is given, --sequence is used by default.

Usage:
    python scripts/phylogeny_analysis.py <fasta_file> [--feature] [--sequence] [--save]
"""
import argparse
import os
import sys

from Bio import SeqIO

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.processing import process_fasta_to_matrix
from aaseq.visualization import perform_phylogenetic_analysis, perform_sequence_phylogeny


def main():
    parser = argparse.ArgumentParser(description="Phylogenetic analysis of sequences in a FASTA file.")
    parser.add_argument("fasta_file", help="Path to the FASTA file (.fasta)")
    parser.add_argument("--feature", action="store_true", help="Feature-based dendrogram (motif/property similarity)")
    parser.add_argument("--sequence", action="store_true", help="Alignment-based dendrogram (BLOSUM62 pairwise alignment)")
    parser.add_argument("--min_len", type=int, default=3, help="Minimum motif length (default: 3)")
    parser.add_argument("--max_len", type=int, default=10, help="Maximum motif length (default: 10)")
    parser.add_argument("--min_reps", type=int, default=2, help="Minimum repetition count (default: 2)")
    parser.add_argument("--workers", type=int, default=4, help="Number of CPU workers (default: 4)")
    parser.add_argument("--chunk_size", type=int, default=5000, help="Chunk size for parsing (default: 5000)")
    parser.add_argument("--save", action="store_true", help="Save the feature matrix to 'repeat_matrix.csv' (only with --feature)")
    args = parser.parse_args()

    if not args.feature and not args.sequence:
        args.sequence = True

    if args.feature:
        matrix = process_fasta_to_matrix(
            args.fasta_file, args.min_len, args.max_len, args.min_reps, args.workers, args.chunk_size
        )
        if matrix is not None:
            perform_phylogenetic_analysis(matrix)
            if args.save:
                matrix.to_csv("repeat_matrix.csv")
                print("Saved successfully as 'repeat_matrix.csv'!")

    if args.sequence:
        records = list(SeqIO.parse(args.fasta_file, "fasta"))
        ids = [r.id for r in records]
        sequences = [str(r.seq) for r in records]
        perform_sequence_phylogeny(ids, sequences)


if __name__ == "__main__":
    main()
