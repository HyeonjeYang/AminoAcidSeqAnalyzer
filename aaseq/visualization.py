import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner, substitution_matrices
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.cluster import hierarchy
from scipy.spatial import distance

from .motifs import getrep
from .properties import (
    calculate_aa_properties,
    get_hydrophobicity_profile,
    generate_point_mutants,
    generate_systematic_single_mutants,
)
from .structure import predict_mean_plddt, predict_plddt_profile


def visualize_matrix(matrix, title="Exact Repeat Motif Frequency Heatmap", save_path=None):
    plt.figure(figsize=(10, 6))
    sns.heatmap(matrix, cmap="mako", linewidths=0.5)
    plt.xlabel("Motif (repeats)")
    plt.ylabel("Protein ID")
    plt.title(title)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=120)
    plt.show()


def plot_idr_profile(scores, record_id="Sequence", regions=None, save_path=None):
    """Plots disorder propensity and low-complexity entropy tracks."""
    plt.figure(figsize=(12, 5))
    plt.plot(scores["position"], scores["disorder_score"], label="Disorder score", color="darkorange")
    plt.plot(scores["position"], scores["normalized_entropy"], label="Normalized entropy", color="steelblue", alpha=0.8)
    plt.axhline(0.2, color="darkorange", linestyle="--", linewidth=0.8)
    plt.axhline(0.55, color="steelblue", linestyle="--", linewidth=0.8)

    if regions is not None and not regions.empty:
        for _, row in regions.iterrows():
            color = "orange" if row["type"] == "IDR" else "skyblue"
            plt.axvspan(row["start"], row["end"], color=color, alpha=0.18)

    plt.title(f"IDR and Low-Complexity Profile - {record_id}")
    plt.xlabel("Residue Position")
    plt.ylabel("Score")
    plt.ylim(-0.6, 1.05)
    plt.legend()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=120)
    plt.show()


def plot_motif_tracks(track_df, sequence_length, record_id="Sequence", save_path=None):
    """Plots repeated motif coordinates as horizontal tracks."""
    if track_df.empty:
        print("No motif occurrences available for track plotting.")
        return

    motifs = track_df["motif"].drop_duplicates().tolist()
    motif_to_y = {motif: idx for idx, motif in enumerate(motifs)}
    plt.figure(figsize=(12, max(3, 0.35 * len(motifs) + 1.5)))
    for _, row in track_df.iterrows():
        y = motif_to_y[row["motif"]]
        plt.plot([row["start"], row["end"]], [y, y], linewidth=6, solid_capstyle="butt")

    plt.yticks(range(len(motifs)), motifs)
    plt.xlim(1, max(sequence_length, 1))
    plt.xlabel("Residue Position")
    plt.ylabel("Motif")
    plt.title(f"Repeated Motif Position Tracks - {record_id}")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=120)
    plt.show()


def plot_conservation(conservation_df, title="MSA Conservation", save_path=None):
    """Plots MSA conservation and gap fraction by alignment column."""
    if conservation_df.empty:
        print("No conservation data available.")
        return

    plt.figure(figsize=(12, 4))
    plt.plot(conservation_df["alignment_position"], conservation_df["conservation"], label="Conservation", color="teal")
    plt.plot(conservation_df["alignment_position"], conservation_df["gap_fraction"], label="Gap fraction", color="gray", alpha=0.8)
    plt.ylim(0, 1.05)
    plt.xlabel("Alignment Position")
    plt.ylabel("Fraction")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=120)
    plt.show()


def plot_plddt_profile(sequence, record_id="Sequence", save_path=None):
    """Queries ESMFold and plots per-residue pLDDT, chunking long sequences."""
    profile = predict_plddt_profile(sequence)
    if profile is None:
        print("No pLDDT profile obtained.")
        return None

    positions = np.arange(1, len(profile) + 1)
    plt.figure(figsize=(12, 4))
    plt.plot(positions, profile, color="purple")
    plt.ylim(0, 100)
    plt.xlabel("Residue Position")
    plt.ylabel("pLDDT")
    plt.title(f"ESMFold pLDDT Profile - {record_id}")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=120)
    plt.show()

    return pd.DataFrame({"position": positions, "pLDDT": profile})


def perform_pca_clustering(matrix, n_clusters=3, save_paths=None):
    print("\nPerforming PCA and Clustering...")
    # Normalize motif frequencies and property ratios which have different units
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(matrix)

    # PCA dimensionality reduction (for 2D visualization)
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)

    # KMeans clustering (groups proteins with similar properties)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(scaled_data)
    matrix['Cluster'] = clusters

    # PCA scatter plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1], hue=clusters, palette="viridis", s=80, alpha=0.7)
    plt.title("PCA & Clustering of Protein Sequences")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.tight_layout()
    if save_paths:
        plt.savefig(save_paths[0], dpi=120)
    plt.show()

    # Hierarchical clustering-based heatmap visualization
    print("Generating Hierarchical Clustermap...")
    sns.clustermap(scaled_data, cmap="mako", figsize=(12, 10))
    plt.title("Hierarchical Clustermap")
    if save_paths:
        plt.savefig(save_paths[1], dpi=120)
    plt.show()

    return matrix


def perform_phylogenetic_analysis(matrix, save_path=None):
    print("\nPerforming Feature-based Phylogenetic Analysis (Dendrogram)...")
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(matrix)

    # Calculate distance matrix and linkage
    dist_matrix = distance.pdist(scaled_data, metric='euclidean')
    linkage_matrix = hierarchy.linkage(dist_matrix, method='ward')

    plt.figure(figsize=(10, 6))
    hierarchy.dendrogram(linkage_matrix, labels=matrix.index.tolist(), leaf_rotation=90)
    plt.title("Feature-based Phylogenetic Tree")
    plt.ylabel("Distance")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=120)
    plt.show()


def _generate_mutation_panel(
    sequence,
    num_mutations=20,
    mutation_mode="random",
    seed=None,
    positions=None,
    target_motif=None,
):
    if mutation_mode == "systematic":
        return generate_systematic_single_mutants(sequence, positions=positions, target_motif=target_motif)
    return generate_point_mutants(
        sequence,
        num_mutations=num_mutations,
        seed=seed,
        positions=positions,
        target_motif=target_motif,
    )


def in_silico_mutagenesis(
    sequence,
    min_len,
    max_len,
    min_reps,
    num_mutations=20,
    save_path=None,
    mutation_mode="random",
    seed=None,
    positions=None,
    target_motif=None,
):
    print("\nPerforming In Silico Mutagenesis (Random point mutations)...")
    mutant_ids, mutants = _generate_mutation_panel(
        sequence, num_mutations, mutation_mode, seed, positions, target_motif
    )

    all_data = []
    for mutant_id, seq in zip(mutant_ids, mutants):
        repeat_counts = getrep(seq, min_len, max_len, min_reps)
        properties = calculate_aa_properties(seq)
        result_dict = {**repeat_counts, **properties}
        all_data.append(pd.DataFrame([result_dict], index=[mutant_id]))

    matrix = pd.concat(all_data).fillna(0)

    plt.figure(figsize=(12, 8))
    sns.heatmap(matrix, cmap="vlag", linewidths=0.5)
    plt.title("In Silico Mutagenesis Profile Heatmap")
    plt.xlabel("Features (Motifs & Properties)")
    plt.ylabel("Mutants")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=120)
    plt.show()
    return matrix


def plot_hydrophobicity(sequence, record_id="Sequence", window=9, save_path=None):
    """Plots a Kyte-Doolittle hydrophobicity profile for a single sequence."""
    positions, scores = get_hydrophobicity_profile(sequence, window)
    plt.figure(figsize=(10, 4))
    plt.plot(positions, scores, color="teal")
    plt.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    plt.title(f"Hydrophobicity Profile (Kyte-Doolittle, window={window}) - {record_id}")
    plt.xlabel("Residue Position")
    plt.ylabel("Hydrophobicity Score")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=120)
    plt.show()


def _build_aligner():
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"
    return aligner


def compute_alignment_distance_matrix(sequences):
    """Computes a pairwise distance matrix (1 - normalized BLOSUM62 alignment score)."""
    aligner = _build_aligner()
    n = len(sequences)
    self_scores = [aligner.score(seq, seq) for seq in sequences]

    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            score = aligner.score(sequences[i], sequences[j])
            max_self = max(self_scores[i], self_scores[j])
            similarity = score / max_self if max_self else 0
            d = max(1 - similarity, 0)
            dist[i, j] = dist[j, i] = d
    return dist


def perform_sequence_phylogeny(ids, sequences, title="Sequence Alignment-based Phylogenetic Tree", save_path=None):
    """Builds a phylogenetic tree from real pairwise alignment distances (BLOSUM62)."""
    print(f"\nPerforming Alignment-based Phylogenetic Analysis ({len(sequences)} sequences)...")
    dist_matrix = compute_alignment_distance_matrix(sequences)
    condensed = distance.squareform(dist_matrix, checks=False)
    linkage_matrix = hierarchy.linkage(condensed, method='average')

    plt.figure(figsize=(10, 6))
    hierarchy.dendrogram(linkage_matrix, labels=ids, leaf_rotation=90)
    plt.title(title)
    plt.ylabel("Distance (1 - normalized alignment score)")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=120)
    plt.show()


def mutagenesis_phylogeny(
    sequence,
    num_mutations=20,
    save_path=None,
    mutation_mode="random",
    seed=None,
    positions=None,
    target_motif=None,
):
    """Generates random point mutants and builds an alignment-based phylogenetic tree,
    showing how individual mutations shift each mutant relative to the wild type."""
    mutant_ids, mutants = _generate_mutation_panel(
        sequence, num_mutations, mutation_mode, seed, positions, target_motif
    )
    perform_sequence_phylogeny(
        mutant_ids, mutants,
        title="In Silico Mutagenesis Phylogenetic Tree",
        save_path=save_path
    )
    return mutant_ids, mutants


def mutagenesis_plddt(
    sequence,
    num_mutations=20,
    save_path=None,
    mutation_mode="random",
    seed=None,
    positions=None,
    target_motif=None,
):
    """Generates random point mutants and queries the ESMFold API for each one's mean
    pLDDT, then plots a bar chart comparing structural confidence to the wild type.

    Requires an internet connection. Sequences longer than the ESMFold API limit
    (400 residues) are truncated by aaseq.structure.predict_structure.
    """
    mutant_ids, mutants = _generate_mutation_panel(
        sequence, num_mutations, mutation_mode, seed, positions, target_motif
    )

    print(f"\nQuerying ESMFold for {len(mutants)} sequences (this may take a while)...")
    results = []
    for mutant_id, seq in zip(mutant_ids, mutants):
        print(f"  Folding {mutant_id}...")
        mean_plddt = predict_mean_plddt(seq)
        if mean_plddt is not None:
            results.append((mutant_id, mean_plddt))

    if not results:
        print("  No pLDDT results obtained (ESMFold API unreachable). Skipping plot.")
        return None

    plddt_df = pd.DataFrame(results, columns=["Mutant", "Mean_pLDDT"]).set_index("Mutant")

    plt.figure(figsize=(10, 6))
    colors = ["steelblue" if idx != "Wild_Type" else "darkorange" for idx in plddt_df.index]
    plt.bar(plddt_df.index, plddt_df["Mean_pLDDT"], color=colors)
    plt.axhline(plddt_df.loc["Wild_Type", "Mean_pLDDT"], color="darkorange", linestyle="--", linewidth=0.8, label="Wild Type")
    plt.title("ESMFold Mean pLDDT: Wild Type vs Point Mutants")
    plt.xlabel("Sequence")
    plt.ylabel("Mean pLDDT")
    plt.xticks(rotation=90)
    plt.legend()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=120)
    plt.show()

    return plddt_df
