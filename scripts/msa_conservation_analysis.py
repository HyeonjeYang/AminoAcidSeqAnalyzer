"""Standalone: reference-guided MSA and conservation profile."""
import argparse
import os
import sys

from Bio import SeqIO

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.alignment import conservation_dataframe, consensus_sequence, reference_msa, write_alignment_fasta
from aaseq.visualization import plot_conservation


def main():
    parser = argparse.ArgumentParser(description="Reference-guided MSA conservation analysis.")
    parser.add_argument("fasta_file", help="Path to FASTA file")
    parser.add_argument("--reference_id", default=None, help="Reference record ID (default: longest sequence)")
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.fasta_file, "fasta"))
    ids = [record.id for record in records]
    sequences = [str(record.seq) for record in records]
    alignment = reference_msa(ids, sequences, reference_id=args.reference_id)
    conservation = conservation_dataframe(alignment)
    consensus = consensus_sequence(alignment)
    print(f"MSA records: {len(alignment)}")
    print(f"Alignment length: {len(alignment[0][1]) if alignment else 0}")
    print(f"Consensus preview: {consensus[:120]}{'...' if len(consensus) > 120 else ''}")
    print(conservation.head(20))
    if args.save:
        write_alignment_fasta(alignment, "msa_alignment.fasta")
        conservation.to_csv("msa_conservation.csv", index=False)
    plot_conservation(conservation)


if __name__ == "__main__":
    main()
