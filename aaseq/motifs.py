from collections import Counter


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
