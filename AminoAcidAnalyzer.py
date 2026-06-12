import argparse
import concurrent.futures
import random
from collections import Counter
from itertools import islice

from tqdm import tqdm
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.cluster import hierarchy
from scipy.spatial import distance

def getrep(sequence, min_len=3, max_len=10, min_reps=2):
    """
    sequence: Amino acid sequence to analyze (string)
    min_len: Minimum motif length
    max_len: Maximum motif length
    min_reps: Minimum repetition count
    """
    sequence = sequence.replace("\n", "").strip()
    count_dict = {}
    
    for length in range(min_len, max_len + 1):
        motifs_gen = (sequence[i:i + length] for i in range(0, len(sequence) - length + 1))
        motif_counts = Counter(motifs_gen)
        for motif, c in motif_counts.items():
            if c >= min_reps:
                count_dict[motif] = c
    return count_dict

def calculate_aa_properties(sequence):
    """Calculates the composition ratio (%) of chemical properties of amino acids in the sequence."""
    groups = {
        "Non-polar": set("AILMFVPG"),
        "Polar": set("STCNQY"),
        "Aromatic": set("FYW"),
        "Positive": set("KRH"),
        "Negative": set("DE")
    }
    seq_len = len(sequence) if len(sequence) > 0 else 1
    props = {}
    for prop_name, aa_set in groups.items():
        count = sum(1 for aa in sequence if aa in aa_set)
        props[f"Prop_{prop_name}"] = round((count / seq_len) * 100, 2)
    return props

def process_record(record_id, seq, min_len, max_len, min_reps):
    """Worker function for single sequence analysis (for parallel processing)"""
    repeat_counts = getrep(seq, min_len, max_len, min_reps)
    properties = calculate_aa_properties(seq)
    
    # Merge motif search results and chemical property data
    result_dict = {**repeat_counts, **properties}
    return pd.DataFrame([result_dict], index=[record_id])

def process_fasta_to_matrix(infile, min_len=3, max_len=10, min_reps=2, workers=4, chunk_size=5000):
    all_data = []
    record_iter = SeqIO.parse(infile, "fasta")
    
    # Implement parallel processing to maximize CPU core utilization using ProcessPoolExecutor
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        chunk_idx = 1
        while True:
            # Load only chunk_size items into memory using islice to prevent overflow of large files
            batch = list(islice(record_iter, chunk_size))
            if not batch:
                break
                
            futures = []
            for record in batch:
                # Pass record.id and string sequence instead of SeqRecord object for serialization speed
                futures.append(
                    executor.submit(process_record, record.id, str(record.seq), min_len, max_len, min_reps)
                )
            
            # Apply tqdm to track parallel processing progress by chunk
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc=f"Analyzing chunk {chunk_idx}"):
                res = future.result()
                if res is not None:
                    all_data.append(res)
            chunk_idx += 1

    if not all_data:
        print("No sequences processed.")
        return None

    # Do not cast to int as property data contains floats
    matrix = pd.concat(all_data).fillna(0)
    return matrix

def visualize_matrix(matrix, title="Exact Repeat Motif Frequency Heatmap"):
    plt.figure(figsize=(10, 6))
    sns.heatmap(matrix, cmap="mako", linewidths=0.5)
    plt.xlabel("Motif (repeats)")
    plt.ylabel("Protein ID")
    plt.title(title)
    plt.tight_layout()
    plt.show()

def perform_pca_clustering(matrix, n_clusters=3):
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
    plt.show()
    
    # Hierarchical clustering-based heatmap visualization
    print("Generating Hierarchical Clustermap...")
    sns.clustermap(scaled_data, cmap="mako", figsize=(12, 10))
    plt.title("Hierarchical Clustermap")
    plt.show()
    
    return matrix

def perform_phylogenetic_analysis(matrix):
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
    plt.show()

def in_silico_mutagenesis(sequence, min_len, max_len, min_reps, num_mutations=20):
    print("\nPerforming In Silico Mutagenesis (Random point mutations)...")
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    mutants = [sequence]
    mutant_ids = ["Wild_Type"]
    
    seq_list = list(sequence)
    for i in range(num_mutations):
        idx = random.randint(0, len(seq_list) - 1)
        original_aa = seq_list[idx]
        new_aa = random.choice(amino_acids.replace(original_aa, ""))
        
        mutated_seq = seq_list.copy()
        mutated_seq[idx] = new_aa
        mutants.append("".join(mutated_seq))
        mutant_ids.append(f"Mut_{original_aa}{idx+1}{new_aa}")
        
    all_data = []
    for i, seq in enumerate(mutants):
        repeat_counts = getrep(seq, min_len, max_len, min_reps)
        properties = calculate_aa_properties(seq)
        result_dict = {**repeat_counts, **properties}
        all_data.append(pd.DataFrame([result_dict], index=[mutant_ids[i]]))
        
    matrix = pd.concat(all_data).fillna(0)
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(matrix, cmap="vlag", linewidths=0.5)
    plt.title("In Silico Mutagenesis Profile Heatmap")
    plt.xlabel("Features (Motifs & Properties)")
    plt.ylabel("Mutants")
    plt.tight_layout()
    plt.show()
    return matrix

def main():
    parser = argparse.ArgumentParser(description="Analyze Amino Acid Sequences for Exact Repeats and Properties.")
    parser.add_argument("fasta_file", help="Path to the FASTA file (.fasta)")
    parser.add_argument("--min_len", type=int, default=3, help="Minimum motif length (default: 3)")
    parser.add_argument("--max_len", type=int, default=10, help="Maximum motif length (default: 10)")
    parser.add_argument("--min_reps", type=int, default=2, help="Minimum repetition count (default: 2)")
    parser.add_argument("--workers", type=int, default=4, help="Number of CPU workers for parallel processing (default: 4)")
    parser.add_argument("--chunk_size", type=int, default=5000, help="Chunk size for parsing large FASTA files (default: 5000)")
    parser.add_argument("--pca", action="store_true", help="Perform PCA and Clustering on the analyzed data")
    parser.add_argument("--clusters", type=int, default=3, help="Number of clusters for KMeans (default: 3)")
    parser.add_argument("--phylogeny", action="store_true", help="Perform feature-based phylogenetic analysis (Dendrogram)")
    parser.add_argument("--mutate", action="store_true", help="Perform In Silico Mutagenesis on the first sequence and visualize the effect")
    parser.add_argument("--save", action="store_true", help="Save the output matrix to 'repeat_matrix.csv'")
    
    args = parser.parse_args()
    
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
            visualize_matrix(matrix[motif_cols])
            
        # Perform PCA and Clustering visualization if option is provided
        if args.pca:
            matrix = perform_pca_clustering(matrix, args.clusters)
            
        if args.phylogeny:
            perform_phylogenetic_analysis(matrix)
            
        if args.save:
            matrix.to_csv("repeat_matrix.csv")
            print("Saved successfully as 'repeat_matrix.csv'!")

        if args.mutate:
            first_record = next(SeqIO.parse(args.fasta_file, "fasta"))
            mut_matrix = in_silico_mutagenesis(str(first_record.seq), args.min_len, args.max_len, args.min_reps)
            if args.save:
                mut_matrix.to_csv("mutation_matrix.csv")
                print("Mutation analysis saved as 'mutation_matrix.csv'!")

if __name__ == "__main__":
    main()