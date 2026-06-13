import argparse
import os

from Bio import SeqIO

from aaseq.alignment import conservation_dataframe, consensus_sequence, reference_msa, write_alignment_fasta
from aaseq.annotations import annotate_functional_motifs
from aaseq.disorder import analyze_idr
from aaseq.enrichment import group_feature_enrichment
from aaseq.llps import analyze_llps
from aaseq.motifs import get_fuzzy_repeats, getrep
from aaseq.processing import process_fasta_to_matrix
from aaseq.qc import analyze_fasta_qc
from aaseq.report import find_record, generate_report, print_quick_analysis
from aaseq.statistics import motif_enrichment
from aaseq.tracks import motif_track_dataframe
from aaseq.visualization import (
    in_silico_mutagenesis,
    mutagenesis_phylogeny,
    mutagenesis_plddt,
    perform_pca_clustering,
    perform_phylogenetic_analysis,
    perform_sequence_phylogeny,
    plot_conservation,
    plot_hydrophobicity,
    plot_idr_profile,
    plot_llps_profile,
    plot_motif_tracks,
    plot_plddt_profile,
    visualize_matrix,
)


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
NON_MOTIF_COLUMNS = {
    "NGlyco_Sites",
    "MolecularWeight",
    "IsoelectricPoint",
    "InstabilityIndex",
    "Aromaticity",
    "GRAVY",
    "Cluster",
}


def get_motif_columns(matrix, min_len, max_len):
    motif_cols = []
    for col in matrix.columns:
        if col.startswith("Fuzzy_"):
            motif_cols.append(col)
        elif col.startswith("Prop_") or col.startswith("IDR_") or col.startswith("LLPS_") or col in NON_MOTIF_COLUMNS:
            continue
        elif min_len <= len(col) <= max_len and set(col).issubset(STANDARD_AA):
            motif_cols.append(col)
    return motif_cols


def print_table_preview(title, table, rows=10):
    print(f"\n{title}:")
    if table is None or table.empty:
        print("  None")
    else:
        print(table.head(rows))


def main():
    parser = argparse.ArgumentParser(description="Analyze Amino Acid Sequences for Exact Repeats and Properties.")
    parser.add_argument("fasta_file", nargs="?", help="Path to the FASTA file (.fasta). Not required when using --seq.")
    parser.add_argument("--seq", type=str, help="Directly analyze a single amino acid sequence string (no FASTA file needed)")
    parser.add_argument("--min_len", type=int, default=3, help="Minimum motif length (default: 3)")
    parser.add_argument("--max_len", type=int, default=10, help="Maximum motif length (default: 10)")
    parser.add_argument("--min_reps", type=int, default=2, help="Minimum repetition count (default: 2)")
    parser.add_argument("--workers", type=int, default=4, help="Number of CPU workers for parallel processing (default: 4)")
    parser.add_argument("--chunk_size", type=int, default=5000, help="Chunk size for parsing large FASTA files (default: 5000)")
    parser.add_argument("--fuzzy_motifs", action="store_true", help="Also find approximate repeated motif families")
    parser.add_argument("--fuzzy_mismatches", type=int, default=1, help="Allowed mismatches per fuzzy motif window (default: 1)")
    parser.add_argument("--motif_stats", action="store_true", help="Run shuffle-based statistical enrichment for repeats in the selected sequence")
    parser.add_argument("--num_shuffles", type=int, default=200, help="Number of composition-preserving shuffles for --motif_stats")
    parser.add_argument("--idr", action="store_true", help="Analyze intrinsic-disorder and low-complexity regions")
    parser.add_argument("--disorder_window", type=int, default=15, help="Window size for IDR smoothing (default: 15)")
    parser.add_argument("--complexity_window", type=int, default=12, help="Window size for low-complexity entropy (default: 12)")
    parser.add_argument("--disorder_threshold", type=float, default=0.2, help="Disorder score threshold (default: 0.2)")
    parser.add_argument("--low_complexity_threshold", type=float, default=0.55, help="Normalized entropy threshold for low complexity (default: 0.55)")
    parser.add_argument("--llps", action="store_true", help="Analyze IDR/LLPS-related sequence features and candidate regions")
    parser.add_argument("--llps_window", type=int, default=31, help="Sliding window size for LLPS local features (default: 31)")
    parser.add_argument("--llps_threshold", type=float, default=0.6, help="Local LLPS score threshold for candidate regions (default: 0.6)")
    parser.add_argument("--llps_min_region", type=int, default=20, help="Minimum LLPS candidate region length (default: 20)")
    parser.add_argument("--motif_tracks", action="store_true", help="Plot repeated motif positions along the selected sequence")
    parser.add_argument("--track_top_n", type=int, default=20, help="Number of repeated motifs to show in motif tracks")
    parser.add_argument("--metadata", type=str, help="CSV metadata file for group enrichment")
    parser.add_argument("--id_col", type=str, default="id", help="Metadata sequence ID column (default: id)")
    parser.add_argument("--group_col", type=str, default="group", help="Metadata group column (default: group)")
    parser.add_argument("--group_enrichment", action="store_true", help="Test feature enrichment for metadata groups")
    parser.add_argument("--msa", action="store_true", help="Build a simple reference-guided MSA and conservation profile")
    parser.add_argument("--msa_reference", type=str, default=None, help="Record ID to use as MSA reference (default: longest sequence)")
    parser.add_argument("--annotate", action="store_true", help="Annotate rule-based PTM, NLS, PEST, signal peptide, and transmembrane motifs")
    parser.add_argument("--qc", action="store_true", help="Run FASTA quality-control checks")
    parser.add_argument("--qc_min_length", type=int, default=30, help="Minimum expected sequence length for --qc")
    parser.add_argument("--pca", action="store_true", help="Perform PCA and Clustering on the analyzed data")
    parser.add_argument("--clusters", type=int, default=3, help="Number of clusters for KMeans (default: 3)")
    parser.add_argument("--phylogeny", action="store_true", help="Perform feature-based phylogenetic analysis (Dendrogram)")
    parser.add_argument("--seq_phylogeny", action="store_true", help="Perform alignment-based phylogenetic analysis (BLOSUM62 pairwise alignment) across all sequences")
    parser.add_argument("--mutate", action="store_true", help="Perform In Silico Mutagenesis on the selected sequence and visualize the effect")
    parser.add_argument("--mutate_phylogeny", action="store_true", help="Build an alignment-based phylogenetic tree of random point mutants of the selected sequence")
    parser.add_argument("--mutate_plddt", action="store_true", help="Query ESMFold for mean pLDDT of wild type and mutants; long sequences use chunked pLDDT profiles")
    parser.add_argument("--num_mutations", type=int, default=20, help="Number of random point mutations to generate (default: 20)")
    parser.add_argument("--mutation_mode", choices=["random", "systematic"], default="random", help="Mutation generation mode (default: random)")
    parser.add_argument("--mutation_positions", type=str, help="1-based mutation positions, e.g. '5,10-15'")
    parser.add_argument("--target_motif", type=str, help="Restrict mutations to exact occurrences of this motif")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducible mutation/statistical analyses")
    parser.add_argument("--hydrophobicity", action="store_true", help="Plot a Kyte-Doolittle hydrophobicity profile")
    parser.add_argument("--plddt_profile", action="store_true", help="Plot an ESMFold per-residue pLDDT profile; long sequences are chunked")
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

        if args.fuzzy_motifs:
            fuzzy_counts = get_fuzzy_repeats(
                sequence, args.min_len, args.max_len, args.min_reps, args.fuzzy_mismatches
            )
            print("\nFuzzy repeated motif families:")
            if fuzzy_counts:
                for motif, count in sorted(fuzzy_counts.items(), key=lambda x: -x[1]):
                    print(f"  {motif}: {count}")
            else:
                print("  None found.")

        if args.motif_stats:
            stats = motif_enrichment(
                sequence,
                args.min_len,
                args.max_len,
                args.min_reps,
                args.num_shuffles,
                seed=args.seed,
                observed_counts=repeat_counts,
            )
            print_table_preview("Shuffle-based motif enrichment", stats)

        if args.idr:
            scores, regions, summary = analyze_idr(
                sequence,
                args.disorder_window,
                args.complexity_window,
                args.disorder_threshold,
                args.low_complexity_threshold,
            )
            print("\nIDR / low-complexity summary:")
            for key, value in summary.items():
                print(f"  {key}: {value}")
            print_table_preview("IDR / low-complexity regions", regions)
            plot_idr_profile(scores, "Input Sequence", regions)

        if args.llps:
            llps_summary, llps_prof, llps_regions = analyze_llps(
                sequence,
                window=args.llps_window,
                threshold=args.llps_threshold,
                min_region=args.llps_min_region,
            )
            print("\nLLPS / condensate feature summary:")
            for key, value in llps_summary.items():
                print(f"  {key}: {value}")
            print_table_preview("LLPS candidate regions", llps_regions)
            plot_llps_profile(llps_prof, "Input Sequence", llps_regions)

        if args.motif_tracks:
            track_df = motif_track_dataframe(
                "Input Sequence",
                sequence,
                min_len=args.min_len,
                max_len=args.max_len,
                min_reps=args.min_reps,
                top_n=args.track_top_n,
            )
            print_table_preview("Motif position tracks", track_df)
            plot_motif_tracks(track_df, len(sequence), "Input Sequence")

        if args.annotate:
            annotations = annotate_functional_motifs(sequence)
            print_table_preview("Functional motif annotations", annotations)

        if args.hydrophobicity:
            plot_hydrophobicity(sequence, "Input Sequence")

        if args.plddt_profile:
            plot_plddt_profile(sequence, "Input Sequence")

        if args.mutate_phylogeny:
            mutagenesis_phylogeny(
                sequence,
                args.num_mutations,
                mutation_mode=args.mutation_mode,
                seed=args.seed,
                positions=args.mutation_positions,
                target_motif=args.target_motif,
            )

        if args.mutate_plddt:
            mutagenesis_plddt(
                sequence,
                args.num_mutations,
                mutation_mode=args.mutation_mode,
                seed=args.seed,
                positions=args.mutation_positions,
                target_motif=args.target_motif,
            )
        return

    if not args.fasta_file:
        parser.error("fasta_file is required unless --seq is provided.")

    report_dir = "report"
    figures = []
    tables = []

    if args.report:
        os.makedirs(report_dir, exist_ok=True)

    if args.qc:
        qc_summary, qc_issues = analyze_fasta_qc(args.fasta_file, min_length=args.qc_min_length)
        print("\nFASTA QC summary:")
        for key, value in qc_summary.items():
            print(f"  {key}: {value}")
        print_table_preview("FASTA QC issues", qc_issues)
        if args.save:
            qc_issues.to_csv("qc_issues.csv", index=False)
            print("QC issues saved as 'qc_issues.csv'!")
        if args.report:
            tables.append(("FASTA QC Issues", qc_issues))

    print(f"Analyzing '{args.fasta_file}' using {args.workers} workers...")
    matrix = process_fasta_to_matrix(
        args.fasta_file,
        args.min_len,
        args.max_len,
        args.min_reps,
        args.workers,
        args.chunk_size,
        include_fuzzy=args.fuzzy_motifs,
        fuzzy_mismatches=args.fuzzy_mismatches,
        include_idr=args.idr,
        include_llps=args.llps,
    )

    if matrix is not None:
        print("\nSequence matrix summary (Motifs & Properties):")
        print(matrix.head())

        # Print base motif heatmap excluding property data
        motif_cols = get_motif_columns(matrix, args.min_len, args.max_len)
        if motif_cols:
            save_path = os.path.join(report_dir, "repeat_heatmap.png") if args.report else None
            visualize_matrix(matrix[motif_cols], save_path=save_path)
            if save_path:
                figures.append(("Exact Repeat Motif Frequency Heatmap", save_path))

        target_record = None
        if any([
            args.hydrophobicity, args.mutate, args.mutate_phylogeny, args.mutate_plddt,
            args.idr, args.llps, args.motif_tracks, args.annotate, args.motif_stats, args.plddt_profile,
        ]):
            target_record = find_record(args.fasta_file, args.record_id)
            if target_record is None:
                print(f"\nRecord '{args.record_id}' not found. Selected-sequence analyses skipped.")

        if args.motif_stats and target_record is not None:
            repeat_counts = getrep(str(target_record.seq), args.min_len, args.max_len, args.min_reps)
            stats = motif_enrichment(
                str(target_record.seq),
                args.min_len,
                args.max_len,
                args.min_reps,
                args.num_shuffles,
                seed=args.seed,
                observed_counts=repeat_counts,
            )
            print_table_preview(f"Shuffle-based motif enrichment ({target_record.id})", stats)
            if args.save:
                stats.to_csv("motif_enrichment.csv", index=False)
                print("Motif enrichment saved as 'motif_enrichment.csv'!")
            if args.report:
                tables.append((f"Motif Enrichment ({target_record.id})", stats))

        if args.idr and target_record is not None:
            scores, regions, summary = analyze_idr(
                str(target_record.seq),
                args.disorder_window,
                args.complexity_window,
                args.disorder_threshold,
                args.low_complexity_threshold,
            )
            print("\nIDR / low-complexity summary:")
            for key, value in summary.items():
                print(f"  {key}: {value}")
            print_table_preview(f"IDR / low-complexity regions ({target_record.id})", regions)
            save_path = os.path.join(report_dir, "idr_profile.png") if args.report else None
            plot_idr_profile(scores, target_record.id, regions, save_path=save_path)
            if save_path:
                figures.append((f"IDR / Low-Complexity Profile ({target_record.id})", save_path))
            if args.save:
                scores.to_csv("idr_scores.csv", index=False)
                regions.to_csv("idr_regions.csv", index=False)
                print("IDR scores/regions saved as 'idr_scores.csv' and 'idr_regions.csv'!")
            if args.report:
                tables.append((f"IDR Regions ({target_record.id})", regions))

        if args.llps and target_record is not None:
            llps_summary, llps_prof, llps_regions = analyze_llps(
                str(target_record.seq),
                window=args.llps_window,
                threshold=args.llps_threshold,
                min_region=args.llps_min_region,
            )
            print("\nLLPS / condensate feature summary:")
            for key, value in llps_summary.items():
                print(f"  {key}: {value}")
            print_table_preview(f"LLPS candidate regions ({target_record.id})", llps_regions)
            save_path = os.path.join(report_dir, "llps_profile.png") if args.report else None
            plot_llps_profile(llps_prof, target_record.id, llps_regions, save_path=save_path)
            if save_path:
                figures.append((f"LLPS Candidate Profile ({target_record.id})", save_path))
            if args.save:
                llps_prof.to_csv("llps_profile.csv", index=False)
                llps_regions.to_csv("llps_regions.csv", index=False)
                print("LLPS profile/regions saved as 'llps_profile.csv' and 'llps_regions.csv'!")
            if args.report:
                tables.append((f"LLPS Candidate Regions ({target_record.id})", llps_regions))

        if args.motif_tracks and target_record is not None:
            track_df = motif_track_dataframe(
                target_record.id,
                str(target_record.seq),
                min_len=args.min_len,
                max_len=args.max_len,
                min_reps=args.min_reps,
                top_n=args.track_top_n,
            )
            print_table_preview(f"Motif position tracks ({target_record.id})", track_df)
            save_path = os.path.join(report_dir, "motif_tracks.png") if args.report else None
            plot_motif_tracks(track_df, len(target_record.seq), target_record.id, save_path=save_path)
            if save_path:
                figures.append((f"Motif Position Tracks ({target_record.id})", save_path))
            if args.save:
                track_df.to_csv("motif_tracks.csv", index=False)
                print("Motif tracks saved as 'motif_tracks.csv'!")
            if args.report:
                tables.append((f"Motif Tracks ({target_record.id})", track_df))

        if args.annotate and target_record is not None:
            annotations = annotate_functional_motifs(str(target_record.seq))
            print_table_preview(f"Functional motif annotations ({target_record.id})", annotations)
            if args.save:
                annotations.to_csv("annotations.csv", index=False)
                print("Annotations saved as 'annotations.csv'!")
            if args.report:
                tables.append((f"Functional Annotations ({target_record.id})", annotations))

        if args.group_enrichment:
            if not args.metadata:
                parser.error("--group_enrichment requires --metadata.")
            enrichment = group_feature_enrichment(matrix, args.metadata, args.id_col, args.group_col)
            print_table_preview("Group feature enrichment", enrichment)
            if args.save:
                enrichment.to_csv("group_enrichment.csv", index=False)
                print("Group enrichment saved as 'group_enrichment.csv'!")
            if args.report:
                tables.append(("Group Feature Enrichment", enrichment))

        if args.msa:
            records = list(SeqIO.parse(args.fasta_file, "fasta"))
            ids = [record.id for record in records]
            sequences = [str(record.seq) for record in records]
            alignment = reference_msa(ids, sequences, reference_id=args.msa_reference)
            conservation = conservation_dataframe(alignment)
            consensus = consensus_sequence(alignment)
            print(f"\nMSA built for {len(alignment)} sequences; alignment length {len(alignment[0][1]) if alignment else 0}.")
            print(f"Consensus preview: {consensus[:120]}{'...' if len(consensus) > 120 else ''}")
            print_table_preview("MSA conservation", conservation)
            if args.save or args.report:
                alignment_path = os.path.join(report_dir, "msa_alignment.fasta") if args.report else "msa_alignment.fasta"
                conservation_path = os.path.join(report_dir, "msa_conservation.csv") if args.report else "msa_conservation.csv"
                write_alignment_fasta(alignment, alignment_path)
                conservation.to_csv(conservation_path, index=False)
                print(f"MSA alignment saved as '{alignment_path}'!")
                print(f"MSA conservation saved as '{conservation_path}'!")
            save_path = os.path.join(report_dir, "msa_conservation.png") if args.report else None
            plot_conservation(conservation, save_path=save_path)
            if save_path:
                figures.append(("MSA Conservation", save_path))
            if args.report:
                tables.append(("MSA Conservation", conservation))

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
            if target_record is None:
                print(f"\nRecord '{args.record_id}' not found for hydrophobicity profile.")
            else:
                save_path = os.path.join(report_dir, "hydrophobicity.png") if args.report else None
                plot_hydrophobicity(str(target_record.seq), target_record.id, save_path=save_path)
                if save_path:
                    figures.append((f"Hydrophobicity Profile ({target_record.id})", save_path))

        if args.plddt_profile and target_record is not None:
            save_path = os.path.join(report_dir, "plddt_profile.png") if args.report else None
            plddt_profile = plot_plddt_profile(str(target_record.seq), target_record.id, save_path=save_path)
            if plddt_profile is not None:
                if args.save:
                    plddt_profile.to_csv("plddt_profile.csv", index=False)
                    print("pLDDT profile saved as 'plddt_profile.csv'!")
                if args.report:
                    figures.append((f"ESMFold pLDDT Profile ({target_record.id})", save_path))
                    tables.append((f"pLDDT Profile ({target_record.id})", plddt_profile))

        if args.save:
            matrix.to_csv("repeat_matrix.csv")
            print("Saved successfully as 'repeat_matrix.csv'!")

        mut_matrix = None
        if args.mutate:
            if target_record is None:
                print("\nSkipping mutagenesis because no target record was selected.")
            else:
                save_path = os.path.join(report_dir, "mutagenesis_heatmap.png") if args.report else None
                mut_matrix = in_silico_mutagenesis(
                    str(target_record.seq),
                    args.min_len,
                    args.max_len,
                    args.min_reps,
                    args.num_mutations,
                    save_path=save_path,
                    mutation_mode=args.mutation_mode,
                    seed=args.seed,
                    positions=args.mutation_positions,
                    target_motif=args.target_motif,
                )
                if save_path:
                    figures.append((f"In Silico Mutagenesis Profile ({target_record.id})", save_path))
                if args.save:
                    mut_matrix.to_csv("mutation_matrix.csv")
                    print("Mutation analysis saved as 'mutation_matrix.csv'!")

        if args.mutate_phylogeny:
            if target_record is None:
                print("\nSkipping mutagenesis phylogeny because no target record was selected.")
            else:
                save_path = os.path.join(report_dir, "mutagenesis_phylogeny.png") if args.report else None
                mutagenesis_phylogeny(
                    str(target_record.seq),
                    args.num_mutations,
                    save_path=save_path,
                    mutation_mode=args.mutation_mode,
                    seed=args.seed,
                    positions=args.mutation_positions,
                    target_motif=args.target_motif,
                )
                if save_path:
                    figures.append((f"In Silico Mutagenesis Phylogenetic Tree ({target_record.id})", save_path))

        if args.mutate_plddt:
            if target_record is None:
                print("\nSkipping mutagenesis pLDDT because no target record was selected.")
            else:
                save_path = os.path.join(report_dir, "mutagenesis_plddt.png") if args.report else None
                plddt_df = mutagenesis_plddt(
                    str(target_record.seq),
                    args.num_mutations,
                    save_path=save_path,
                    mutation_mode=args.mutation_mode,
                    seed=args.seed,
                    positions=args.mutation_positions,
                    target_motif=args.target_motif,
                )
                if plddt_df is not None and save_path:
                    figures.append((f"ESMFold pLDDT Comparison ({target_record.id})", save_path))
                if plddt_df is not None and args.save:
                    plddt_df.to_csv("mutation_plddt.csv")
                    print("pLDDT comparison saved as 'mutation_plddt.csv'!")

        if args.report:
            generate_report(report_dir, args.fasta_file, matrix, figures, mut_matrix, tables)

if __name__ == "__main__":
    main()
