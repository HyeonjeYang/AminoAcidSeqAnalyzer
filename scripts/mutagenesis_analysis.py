"""Standalone: in silico random point-mutagenesis analysis.

Pick one or more analysis types:
  --heatmap   : motif/property heatmap of wild type vs mutants
  --phylogeny : alignment-based phylogenetic tree of wild type vs mutants
  --plddt     : ESMFold mean pLDDT comparison of wild type vs mutants (requires internet, <=400aa)

Usage:
    python scripts/mutagenesis_analysis.py <fasta_file> --record_id ID --heatmap --phylogeny --plddt
    python scripts/mutagenesis_analysis.py --seq "MKT..." --phylogeny
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.report import find_record
from aaseq.visualization import in_silico_mutagenesis, mutagenesis_phylogeny, mutagenesis_plddt


def main():
    parser = argparse.ArgumentParser(description="In silico random point-mutagenesis analysis.")
    parser.add_argument("fasta_file", nargs="?", help="Path to the FASTA file (.fasta). Not required when using --seq.")
    parser.add_argument("--seq", type=str, help="Directly analyze a single amino acid sequence string (no FASTA file needed)")
    parser.add_argument("--record_id", type=str, default=None, help="Record ID to mutate (default: first record)")
    parser.add_argument("--num_mutations", type=int, default=20, help="Number of random point mutations to generate (default: 20)")
    parser.add_argument("--min_len", type=int, default=3, help="Minimum motif length for --heatmap (default: 3)")
    parser.add_argument("--max_len", type=int, default=10, help="Maximum motif length for --heatmap (default: 10)")
    parser.add_argument("--min_reps", type=int, default=2, help="Minimum repetition count for --heatmap (default: 2)")
    parser.add_argument("--heatmap", action="store_true", help="Motif/property heatmap of wild type vs mutants")
    parser.add_argument("--phylogeny", action="store_true", help="Alignment-based phylogenetic tree of wild type vs mutants")
    parser.add_argument("--plddt", action="store_true", help="ESMFold mean pLDDT comparison of wild type vs mutants")
    parser.add_argument("--save", action="store_true", help="Save results to 'mutation_matrix.csv' / 'mutation_plddt.csv'")
    args = parser.parse_args()

    if args.seq:
        sequence = args.seq.replace("\n", "").strip().upper()
    elif args.fasta_file:
        record = find_record(args.fasta_file, args.record_id)
        if record is None:
            parser.error(f"Record '{args.record_id}' not found in '{args.fasta_file}'.")
        sequence = str(record.seq)
    else:
        parser.error("Either fasta_file or --seq must be provided.")

    if not (args.heatmap or args.phylogeny or args.plddt):
        parser.error("Specify at least one of --heatmap, --phylogeny, --plddt.")

    if args.heatmap:
        mut_matrix = in_silico_mutagenesis(sequence, args.min_len, args.max_len, args.min_reps, args.num_mutations)
        if args.save:
            mut_matrix.to_csv("mutation_matrix.csv")
            print("Saved successfully as 'mutation_matrix.csv'!")

    if args.phylogeny:
        mutagenesis_phylogeny(sequence, args.num_mutations)

    if args.plddt:
        plddt_df = mutagenesis_plddt(sequence, args.num_mutations)
        if plddt_df is not None and args.save:
            plddt_df.to_csv("mutation_plddt.csv")
            print("Saved successfully as 'mutation_plddt.csv'!")


if __name__ == "__main__":
    main()
