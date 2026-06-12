import argparse
import concurrent.futures
from collections import Counter

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO

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

def process_record(record_id, seq, min_len, max_len, min_reps):
    """단일 서열 분석을 담당하는 워커 함수 (병렬 처리용)"""
    repeat_counts = getrep(seq, min_len, max_len, min_reps)
    if repeat_counts:
        return pd.DataFrame([repeat_counts], index=[record_id])
    return None

def process_fasta_to_matrix(infile, min_len=3, max_len=10, min_reps=2, workers=4):
    all_data = []
    
    # ProcessPoolExecutor를 이용해 CPU 코어를 최대로 활용하는 병렬 처리 구현
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        futures = []
        for record in SeqIO.parse(infile, "fasta"):
            # 직렬화 속도를 위해 객체 대신 record.id 와 str 서열만 전달
            futures.append(
                executor.submit(process_record, record.id, str(record.seq), min_len, max_len, min_reps)
            )
        
        for future in concurrent.futures.as_completed(futures):
            res = future.result()
            if res is not None:
                all_data.append(res)

    if not all_data:
        print("No repeated sequences found. Please relax the conditions.")
        return None

    matrix = pd.concat(all_data).fillna(0).astype(int)
    return matrix

def visualize_matrix(matrix, title="Exact Repeat Motif Frequency Heatmap"):
    plt.figure(figsize=(10, 6))
    sns.heatmap(matrix, cmap="mako", linewidths=0.5)
    plt.xlabel("Motif (repeats)")
    plt.ylabel("Protein ID")
    plt.title(title)
    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Analyze Amino Acid Sequences for Exact Repeats.")
    parser.add_argument("fasta_file", help="Path to the FASTA file (.fasta)")
    parser.add_argument("--min_len", type=int, default=3, help="Minimum motif length (default: 3)")
    parser.add_argument("--max_len", type=int, default=10, help="Maximum motif length (default: 10)")
    parser.add_argument("--min_reps", type=int, default=2, help="Minimum repetition count (default: 2)")
    parser.add_argument("--workers", type=int, default=4, help="Number of CPU workers for parallel processing (default: 4)")
    parser.add_argument("--save", action="store_true", help="Save the output matrix to 'repeat_matrix.csv'")
    
    args = parser.parse_args()
    
    print(f"Analyzing '{args.fasta_file}' using {args.workers} workers...")
    matrix = process_fasta_to_matrix(
        args.fasta_file, 
        args.min_len, 
        args.max_len, 
        args.min_reps, 
        args.workers
    )
    
    if matrix is not None:
        print("\nRepeated sequence matrix summary:")
        print(matrix.head())
        visualize_matrix(matrix)
        if args.save:
            matrix.to_csv("repeat_matrix.csv")
            print("Saved successfully as 'repeat_matrix.csv'!")

if __name__ == "__main__":
    main()