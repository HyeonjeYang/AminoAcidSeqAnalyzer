"""Standalone: rule-based PTM and functional motif annotation."""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.annotations import annotate_functional_motifs
from aaseq.report import find_record


def main():
    parser = argparse.ArgumentParser(description="Rule-based protein motif annotation.")
    parser.add_argument("fasta_file", nargs="?", help="Path to FASTA file. Not required with --seq.")
    parser.add_argument("--seq", type=str, help="Directly analyze one amino acid sequence")
    parser.add_argument("--record_id", default=None)
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

    annotations = annotate_functional_motifs(sequence)
    print(annotations if not annotations.empty else "No rule-based annotations found.")
    if args.save:
        annotations.to_csv(f"{label}_annotations.csv", index=False)


if __name__ == "__main__":
    main()
