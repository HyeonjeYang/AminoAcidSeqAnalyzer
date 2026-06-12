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
from .properties import calculate_aa_properties, get_hydrophobicity_profile, generate_point_mutants


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


def in_silico_mutagenesis(sequence, min_len, max_len, min_reps, num_mutations=20, save_path=None):
    print("\nPerforming In Silico Mutagenesis (Random point mutations)...")
    mutant_ids, mutants = generate_point_mutants(sequence, num_mutations)

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


def mutagenesis_phylogeny(sequence, num_mutations=20, save_path=None):
    """Generates random point mutants and builds an alignment-based phylogenetic tree,
    showing how individual mutations shift each mutant relative to the wild type."""
    mutant_ids, mutants = generate_point_mutants(sequence, num_mutations)
    perform_sequence_phylogeny(
        mutant_ids, mutants,
        title="In Silico Mutagenesis Phylogenetic Tree",
        save_path=save_path
    )
    return mutant_ids, mutants
