import argparse
import os

from Bio import SeqIO

from aaseq.motifs import getrep
from aaseq.processing import process_fasta_to_matrix
from aaseq.report import find_record, generate_report, print_quick_analysis
from aaseq.visualization import (
    in_silico_mutagenesis,
    mutagenesis_phylogeny,
    mutagenesis_plddt,
    perform_pca_clustering,
    perform_phylogenetic_analysis,
    perform_sequence_phylogeny,
    plot_hydrophobicity,
    visualize_matrix,
)


def main():
    parser = argparse.ArgumentParser(description="Analyze Amino Acid Sequences for Exact Repeats and Properties.")
    parser.add_argument("fasta_file", nargs="?", help="Path to the FASTA file (.fasta). Not required when using --seq.")
    parser.add_argument("--seq", type=str, help="Directly analyze a single amino acid sequence string (no FASTA file needed)")
    parser.add_argument("--min_len", type=int, default=3, help="Minimum motif length (default: 3)")
    parser.add_argument("--max_len", type=int, default=10, help="Maximum motif length (default: 10)")
    parser.add_argument("--min_reps", type=int, default=2, help="Minimum repetition count (default: 2)")
    parser.add_argument("--workers", type=int, default=4, help="Number of CPU workers for parallel processing (default: 4)")
    parser.add_argument("--chunk_size", type=int, default=5000, help="Chunk size for parsing large FASTA files (default: 5000)")
    parser.add_argument("--pca", action="store_true", help="Perform PCA and Clustering on the analyzed data")
    parser.add_argument("--clusters", type=int, default=3, help="Number of clusters for KMeans (default: 3)")
    parser.add_argument("--phylogeny", action="store_true", help="Perform feature-based phylogenetic analysis (Dendrogram)")
    parser.add_argument("--seq_phylogeny", action="store_true", help="Perform alignment-based phylogenetic analysis (BLOSUM62 pairwise alignment) across all sequences")
    parser.add_argument("--mutate", action="store_true", help="Perform In Silico Mutagenesis on the selected sequence and visualize the effect")
    parser.add_argument("--mutate_phylogeny", action="store_true", help="Build an alignment-based phylogenetic tree of random point mutants of the selected sequence")
    parser.add_argument("--mutate_plddt", action="store_true", help="Query the ESMFold API for the mean pLDDT of the wild type and each random point mutant (requires internet, <=400aa)")
    parser.add_argument("--num_mutations", type=int, default=20, help="Number of random point mutations to generate (default: 20)")
    parser.add_argument("--hydrophobicity", action="store_true", help="Plot a Kyte-Doolittle hydrophobicity profile")
    parser.add_argument("--record_id", type=str, default=None, help="Record ID to use for hydrophobicity/mutation analysis (default: first record)")
    parser.add_argument("--report", action="store_true", help="Generate an HTML summary report (saved figures + tables) in the 'report' directory")
    parser.add_argument("--save", action="store_true", help="Save the output matrix to 'repeat_matrix.csv'")

    args = parser.parse_args()

    # --- Quick single-sequence mode (no FASTA file required) ---
    if args.seq:
        sequence = args.seq.replace("\n", "").strip().upper()
        print("\n--- Quick Sequence Analysis ---")
        print_quick_analysis(sequence)

        repeat_counts = getrep(sequence, args.min_len, args.max_len, args.min_reps)
        print("\nRepeated motifs (min_reps reached):")
        if repeat_counts:
            for motif, count in sorted(repeat_counts.items(), key=lambda x: -x[1]):
                print(f"  {motif}: {count}")
        else:
            print("  None found.")

        if args.hydrophobicity:
            plot_hydrophobicity(sequence, "Input Sequence")

        if args.mutate_phylogeny:
            mutagenesis_phylogeny(sequence, args.num_mutations)

        if args.mutate_plddt:
            mutagenesis_plddt(sequence, args.num_mutations)
        return

    if not args.fasta_file:
        parser.error("fasta_file is required unless --seq is provided.")

    report_dir = "report"
    figures = []

    print(f"Analyzing '{args.fasta_file}' using {args.workers} workers...")
    matrix = process_fasta_to_matrix(
        args.fasta_file,
        args.min_len,
        args.max_len,
        args.min_reps,
        args.workers,
        args.chunk_size
    )

    if matrix is not None:
        print("\nSequence matrix summary (Motifs & Properties):")
        print(matrix.head())

        # Print base motif heatmap excluding property data
        motif_cols = [c for c in matrix.columns if not c.startswith('Prop_')]
        if motif_cols:
            save_path = os.path.join(report_dir, "repeat_heatmap.png") if args.report else None
            visualize_matrix(matrix[motif_cols], save_path=save_path)
            if save_path:
                figures.append(("Exact Repeat Motif Frequency Heatmap", save_path))

        # Perform PCA and Clustering visualization if option is provided
        if args.pca:
            save_paths = (
                (os.path.join(report_dir, "pca_clusters.png"), os.path.join(report_dir, "clustermap.png"))
                if args.report else None
            )
            matrix = perform_pca_clustering(matrix, args.clusters, save_paths=save_paths)
            if save_paths:
                figures.append(("PCA & Clustering", save_paths[0]))
                figures.append(("Hierarchical Clustermap", save_paths[1]))

        if args.phylogeny:
            save_path = os.path.join(report_dir, "dendrogram.png") if args.report else None
            perform_phylogenetic_analysis(matrix, save_path=save_path)
            if save_path:
                figures.append(("Feature-based Phylogenetic Tree", save_path))

        if args.seq_phylogeny:
            records = list(SeqIO.parse(args.fasta_file, "fasta"))
            ids = [r.id for r in records]
            sequences = [str(r.seq) for r in records]
            save_path = os.path.join(report_dir, "seq_phylogeny.png") if args.report else None
            perform_sequence_phylogeny(ids, sequences, save_path=save_path)
            if save_path:
                figures.append(("Alignment-based Phylogenetic Tree", save_path))

        if args.hydrophobicity:
            target_record = find_record(args.fasta_file, args.record_id)
            if target_record is None:
                print(f"\nRecord '{args.record_id}' not found for hydrophobicity profile.")
            else:
                save_path = os.path.join(report_dir, "hydrophobicity.png") if args.report else None
                plot_hydrophobicity(str(target_record.seq), target_record.id, save_path=save_path)
                if save_path:
                    figures.append((f"Hydrophobicity Profile ({target_record.id})", save_path))

        if args.save:
            matrix.to_csv("repeat_matrix.csv")
            print("Saved successfully as 'repeat_matrix.csv'!")

        mut_matrix = None
        if args.mutate:
            target_record = find_record(args.fasta_file, args.record_id)
            save_path = os.path.join(report_dir, "mutagenesis_heatmap.png") if args.report else None
            mut_matrix = in_silico_mutagenesis(str(target_record.seq), args.min_len, args.max_len, args.min_reps, args.num_mutations, save_path=save_path)
            if save_path:
                figures.append((f"In Silico Mutagenesis Profile ({target_record.id})", save_path))
            if args.save:
                mut_matrix.to_csv("mutation_matrix.csv")
                print("Mutation analysis saved as 'mutation_matrix.csv'!")

        if args.mutate_phylogeny:
            target_record = find_record(args.fasta_file, args.record_id)
            save_path = os.path.join(report_dir, "mutagenesis_phylogeny.png") if args.report else None
            mutagenesis_phylogeny(str(target_record.seq), args.num_mutations, save_path=save_path)
            if save_path:
                figures.append((f"In Silico Mutagenesis Phylogenetic Tree ({target_record.id})", save_path))

        if args.mutate_plddt:
            target_record = find_record(args.fasta_file, args.record_id)
            save_path = os.path.join(report_dir, "mutagenesis_plddt.png") if args.report else None
            plddt_df = mutagenesis_plddt(str(target_record.seq), args.num_mutations, save_path=save_path)
            if plddt_df is not None and save_path:
                figures.append((f"ESMFold pLDDT Comparison ({target_record.id})", save_path))
            if plddt_df is not None and args.save:
                plddt_df.to_csv("mutation_plddt.csv")
                print("pLDDT comparison saved as 'mutation_plddt.csv'!")

        if args.report:
            generate_report(report_dir, args.fasta_file, matrix, figures, mut_matrix)

if __name__ == "__main__":
    main()
