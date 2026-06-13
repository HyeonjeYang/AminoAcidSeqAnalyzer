import math
from collections import Counter

import pandas as pd


# TOP-IDP amino-acid disorder propensity scale.
TOP_IDP = {
    "A": 0.06,
    "R": 0.18,
    "N": 0.007,
    "D": 0.192,
    "C": -0.02,
    "Q": 0.318,
    "E": 0.736,
    "G": 0.166,
    "H": 0.303,
    "I": -0.486,
    "L": -0.326,
    "K": 0.586,
    "M": -0.397,
    "F": -0.697,
    "P": 0.987,
    "S": 0.341,
    "T": 0.059,
    "W": -0.884,
    "Y": -0.51,
    "V": -0.121,
}

DISORDER_PROMOTING = set("PEKQSDRG")


def _window_bounds(index, length, window):
    half = window // 2
    start = max(0, index - half)
    end = min(length, index + half + 1)
    return start, end


def _normalized_entropy(residues):
    if not residues:
        return 0.0
    counts = Counter(residues)
    total = len(residues)
    entropy = -sum((count / total) * math.log2(count / total) for count in counts.values())
    return entropy / math.log2(20)


def _segments_from_mask(mask, min_region):
    segments = []
    start = None
    for idx, flag in enumerate(mask, start=1):
        if flag and start is None:
            start = idx
        elif not flag and start is not None:
            end = idx - 1
            if end - start + 1 >= min_region:
                segments.append((start, end))
            start = None

    if start is not None:
        end = len(mask)
        if end - start + 1 >= min_region:
            segments.append((start, end))
    return segments


def analyze_idr(
    sequence,
    disorder_window=15,
    complexity_window=12,
    disorder_threshold=0.2,
    low_complexity_threshold=0.55,
    min_region=10,
):
    """Scores intrinsic disorder tendency and low-complexity regions.

    This is a fast heuristic, not a replacement for dedicated predictors such as
    IUPred2A. It is useful for prioritizing repeat-rich candidate IDRs locally.
    """
    sequence = sequence.replace("\n", "").strip().upper()
    length = len(sequence)
    rows = []

    for idx, aa in enumerate(sequence):
        disorder_start, disorder_end = _window_bounds(idx, length, disorder_window)
        complexity_start, complexity_end = _window_bounds(idx, length, complexity_window)
        disorder_window_seq = sequence[disorder_start:disorder_end]
        complexity_window_seq = sequence[complexity_start:complexity_end]

        raw_score = TOP_IDP.get(aa, 0.0)
        smoothed = sum(TOP_IDP.get(residue, 0.0) for residue in disorder_window_seq) / len(disorder_window_seq)
        entropy = _normalized_entropy(complexity_window_seq)
        rows.append({
            "position": idx + 1,
            "aa": aa,
            "raw_disorder": raw_score,
            "disorder_score": smoothed,
            "normalized_entropy": entropy,
            "disordered": smoothed >= disorder_threshold,
            "low_complexity": entropy <= low_complexity_threshold,
            "disorder_promoting": aa in DISORDER_PROMOTING,
        })

    scores = pd.DataFrame(rows)
    disorder_regions = _segments_from_mask(scores["disordered"].tolist(), min_region)
    low_complexity_regions = _segments_from_mask(scores["low_complexity"].tolist(), min_region)

    disorder_df = pd.DataFrame([
        {
            "type": "IDR",
            "start": start,
            "end": end,
            "length": end - start + 1,
            "mean_disorder": scores.loc[start - 1:end - 1, "disorder_score"].mean(),
        }
        for start, end in disorder_regions
    ])
    low_complexity_df = pd.DataFrame([
        {
            "type": "LowComplexity",
            "start": start,
            "end": end,
            "length": end - start + 1,
            "mean_entropy": scores.loc[start - 1:end - 1, "normalized_entropy"].mean(),
        }
        for start, end in low_complexity_regions
    ])

    region_tables = [df for df in [disorder_df, low_complexity_df] if not df.empty]
    regions = pd.concat(region_tables, ignore_index=True) if region_tables else pd.DataFrame(
        columns=["type", "start", "end", "length"]
    )

    summary = {
        "Length": length,
        "DisorderFraction": round(float(scores["disordered"].mean()) if length else 0.0, 4),
        "LowComplexityFraction": round(float(scores["low_complexity"].mean()) if length else 0.0, 4),
        "MeanDisorderScore": round(float(scores["disorder_score"].mean()) if length else 0.0, 4),
        "MeanEntropy": round(float(scores["normalized_entropy"].mean()) if length else 0.0, 4),
        "IDR_Count": len(disorder_df),
        "LowComplexity_Count": len(low_complexity_df),
        "Longest_IDR": int(disorder_df["length"].max()) if not disorder_df.empty else 0,
        "Longest_LowComplexity": int(low_complexity_df["length"].max()) if not low_complexity_df.empty else 0,
    }
    return scores, regions, summary
