from collections import Counter


def getrep(sequence, min_len=3, max_len=10, min_reps=2):
    """
    sequence: Amino acid sequence to analyze (string)
    min_len: Minimum motif length
    max_len: Maximum motif length
    min_reps: Minimum repetition count
    """
    sequence = sequence.replace("\n", "").strip().upper()
    count_dict = {}

    for length in range(min_len, max_len + 1):
        motifs_gen = (sequence[i:i + length] for i in range(0, len(sequence) - length + 1))
        motif_counts = Counter(motifs_gen)
        for motif, c in motif_counts.items():
            if c >= min_reps:
                count_dict[motif] = c
    return count_dict


def hamming_distance(seq_a, seq_b):
    """Returns the Hamming distance between equal-length strings."""
    if len(seq_a) != len(seq_b):
        raise ValueError("Hamming distance requires equal-length sequences.")
    return sum(a != b for a, b in zip(seq_a, seq_b))


def get_fuzzy_repeats(sequence, min_len=3, max_len=10, min_reps=2, max_mismatches=1):
    """Finds repeated motif families allowing a small number of mismatches.

    The exact repeat finder treats every window as its own motif. This helper
    clusters same-length windows greedily around the most frequent exact windows,
    then reports a family if the total number of windows within max_mismatches
    reaches min_reps.
    """
    sequence = sequence.replace("\n", "").strip().upper()
    fuzzy_counts = {}

    for length in range(min_len, max_len + 1):
        if length > len(sequence):
            continue

        windows = [sequence[i:i + length] for i in range(0, len(sequence) - length + 1)]
        window_counts = Counter(windows)
        unassigned = set(window_counts)
        representatives = sorted(window_counts, key=lambda motif: (-window_counts[motif], motif))

        for representative in representatives:
            if representative not in unassigned:
                continue

            members = [
                motif for motif in list(unassigned)
                if hamming_distance(representative, motif) <= max_mismatches
            ]
            total = sum(window_counts[motif] for motif in members)
            for motif in members:
                unassigned.remove(motif)

            if total >= min_reps:
                feature = f"Fuzzy_{representative}_len{length}_mm{max_mismatches}"
                fuzzy_counts[feature] = total

    return fuzzy_counts


def get_repeat_positions(sequence, motifs):
    """Returns overlapping exact motif occurrences as 1-based coordinates."""
    sequence = sequence.replace("\n", "").strip().upper()
    rows = []
    for motif in motifs:
        motif = motif.upper()
        start = 0
        while True:
            idx = sequence.find(motif, start)
            if idx == -1:
                break
            rows.append({
                "motif": motif,
                "start": idx + 1,
                "end": idx + len(motif),
                "length": len(motif),
            })
            start = idx + 1
    return rows
