import math
from collections import Counter

import pandas as pd
from Bio.Align import substitution_matrices


def _load_blosum62():
    try:
        return substitution_matrices.load("BLOSUM62")
    except Exception:
        return None


_BLOSUM62 = _load_blosum62()


def _substitution_score(a, b):
    if a == "-" or b == "-":
        return -5
    if _BLOSUM62 is not None:
        try:
            return _BLOSUM62[a, b]
        except Exception:
            pass
    return 1 if a == b else -1


def needleman_wunsch(seq_a, seq_b, gap_penalty=-5):
    """Global pairwise alignment returning two aligned strings."""
    seq_a = seq_a.replace("\n", "").strip().upper()
    seq_b = seq_b.replace("\n", "").strip().upper()
    n = len(seq_a)
    m = len(seq_b)
    scores = [[0.0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        scores[i][0] = scores[i - 1][0] + gap_penalty
        trace[i][0] = "up"
    for j in range(1, m + 1):
        scores[0][j] = scores[0][j - 1] + gap_penalty
        trace[0][j] = "left"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = scores[i - 1][j - 1] + _substitution_score(seq_a[i - 1], seq_b[j - 1])
            up = scores[i - 1][j] + gap_penalty
            left = scores[i][j - 1] + gap_penalty
            best = max(diagonal, up, left)
            scores[i][j] = best
            if best == diagonal:
                trace[i][j] = "diag"
            elif best == up:
                trace[i][j] = "up"
            else:
                trace[i][j] = "left"

    aligned_a = []
    aligned_b = []
    i, j = n, m
    while i > 0 or j > 0:
        direction = trace[i][j]
        if direction == "diag":
            aligned_a.append(seq_a[i - 1])
            aligned_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif direction == "up":
            aligned_a.append(seq_a[i - 1])
            aligned_b.append("-")
            i -= 1
        else:
            aligned_a.append("-")
            aligned_b.append(seq_b[j - 1])
            j -= 1

    return "".join(reversed(aligned_a)), "".join(reversed(aligned_b))


def _pairwise_to_reference_map(reference, aligned_reference, aligned_sequence):
    insertions = ["" for _ in range(len(reference) + 1)]
    aligned_to_reference = ["-" for _ in reference]
    ref_pos = 0

    for ref_char, seq_char in zip(aligned_reference, aligned_sequence):
        if ref_char == "-":
            insertions[ref_pos] += seq_char
        else:
            if ref_pos < len(reference):
                aligned_to_reference[ref_pos] = seq_char
            ref_pos += 1

    return insertions, aligned_to_reference


def reference_msa(ids, sequences, reference_id=None):
    """Builds a simple star MSA around a selected reference sequence.

    This avoids external aligner dependencies and is intended for moderate
    exploratory datasets. For publication-grade MSAs, export FASTA and run a
    dedicated aligner such as MAFFT or Clustal Omega.
    """
    if len(ids) != len(sequences):
        raise ValueError("ids and sequences must have the same length.")
    if not ids:
        return []

    clean_sequences = [seq.replace("\n", "").strip().upper() for seq in sequences]
    if reference_id is None:
        reference_index = max(range(len(clean_sequences)), key=lambda idx: len(clean_sequences[idx]))
    else:
        try:
            reference_index = ids.index(reference_id)
        except ValueError as exc:
            raise ValueError(f"Reference ID '{reference_id}' not found.") from exc

    reference = clean_sequences[reference_index]
    mapped = []
    for seq in clean_sequences:
        aligned_ref, aligned_seq = needleman_wunsch(reference, seq)
        mapped.append(_pairwise_to_reference_map(reference, aligned_ref, aligned_seq))

    max_insertions = [
        max(len(insertions[slot]) for insertions, _ in mapped)
        for slot in range(len(reference) + 1)
    ]

    alignment = []
    for record_id, (insertions, aligned_to_reference) in zip(ids, mapped):
        pieces = []
        for slot in range(len(reference) + 1):
            insertion = insertions[slot]
            pieces.append(insertion + "-" * (max_insertions[slot] - len(insertion)))
            if slot < len(reference):
                pieces.append(aligned_to_reference[slot])
        alignment.append((record_id, "".join(pieces)))

    return alignment


def write_alignment_fasta(alignment, output_path):
    with open(output_path, "w", encoding="utf-8") as handle:
        for record_id, aligned_sequence in alignment:
            handle.write(f">{record_id}\n")
            for i in range(0, len(aligned_sequence), 80):
                handle.write(aligned_sequence[i:i + 80] + "\n")


def conservation_dataframe(alignment):
    """Calculates per-column consensus, conservation, gap fraction, and entropy."""
    if not alignment:
        return pd.DataFrame(columns=[
            "alignment_position", "reference_position", "consensus",
            "conservation", "gap_fraction", "entropy",
        ])

    aligned_sequences = [seq for _, seq in alignment]
    width = len(aligned_sequences[0])
    if any(len(seq) != width for seq in aligned_sequences):
        raise ValueError("All aligned sequences must have the same length.")

    reference_position = 0
    rows = []
    for col_idx in range(width):
        column = [seq[col_idx] for seq in aligned_sequences]
        residues = [aa for aa in column if aa != "-"]
        if alignment[0][1][col_idx] != "-":
            reference_position += 1

        if residues:
            counts = Counter(residues)
            consensus, consensus_count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
            conservation = consensus_count / len(residues)
            entropy = -sum((count / len(residues)) * math.log2(count / len(residues)) for count in counts.values())
        else:
            consensus = "-"
            conservation = 0.0
            entropy = 0.0

        rows.append({
            "alignment_position": col_idx + 1,
            "reference_position": reference_position if alignment[0][1][col_idx] != "-" else None,
            "consensus": consensus,
            "conservation": conservation,
            "gap_fraction": column.count("-") / len(column),
            "entropy": entropy,
        })

    return pd.DataFrame(rows)


def consensus_sequence(alignment, min_conservation=0.5):
    conservation = conservation_dataframe(alignment)
    chars = []
    for _, row in conservation.iterrows():
        if row["consensus"] != "-" and row["conservation"] >= min_conservation:
            chars.append(row["consensus"])
        elif row["gap_fraction"] < 1.0:
            chars.append("X")
    return "".join(chars)
