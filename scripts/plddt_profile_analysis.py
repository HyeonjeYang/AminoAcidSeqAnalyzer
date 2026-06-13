"""Standalone: ESMFold per-residue pLDDT profile.

Long sequences are split into overlapping chunks so a full-length pLDDT profile
can be estimated even when the public ESMFold endpoint rejects >400aa inputs.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.report import find_record
from aaseq.visualization import plot_plddt_profile


def main():
    parser = argparse.ArgumentParser(description="ESMFold pLDDT profile analysis.")
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

    profile = plot_plddt_profile(sequence, label)
    if profile is not None and args.save:
        profile.to_csv(f"{label}_plddt_profile.csv", index=False)


if __name__ == "__main__":
    main()
