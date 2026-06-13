import math
import random

import pandas as pd

from .motifs import getrep


def benjamini_hochberg(p_values):
    """Returns Benjamini-Hochberg FDR-adjusted p-values."""
    n = len(p_values)
    if n == 0:
        return []

    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [1.0] * n
    running_min = 1.0

    for rank, (idx, p_value) in enumerate(reversed(indexed), start=1):
        original_rank = n - rank + 1
        bh_value = p_value * n / original_rank
        running_min = min(running_min, bh_value)
        adjusted[idx] = min(running_min, 1.0)

    return adjusted


def motif_enrichment(
    sequence,
    min_len=3,
    max_len=10,
    min_reps=2,
    num_shuffles=200,
    seed=None,
    observed_counts=None,
):
    """Tests whether exact repeated motifs exceed composition-preserving shuffles.

    The null model shuffles the same residues, preserving amino acid composition
    while removing local order. P-values are one-sided with a +1 correction.
    """
    sequence = sequence.replace("\n", "").strip().upper()
    observed = observed_counts or getrep(sequence, min_len, max_len, min_reps)
    if not observed:
        return pd.DataFrame(columns=[
            "motif", "length", "observed", "expected_mean", "expected_sd",
            "z_score", "p_value", "fdr",
        ])

    rng = random.Random(seed)
    motifs = sorted(observed)
    null_counts = {motif: [] for motif in motifs}
    residues = list(sequence)

    for _ in range(num_shuffles):
        shuffled = residues.copy()
        rng.shuffle(shuffled)
        shuffled_counts = getrep("".join(shuffled), min_len, max_len, min_reps=1)
        for motif in motifs:
            null_counts[motif].append(shuffled_counts.get(motif, 0))

    rows = []
    for motif in motifs:
        values = null_counts[motif]
        expected_mean = sum(values) / len(values) if values else 0.0
        variance = sum((value - expected_mean) ** 2 for value in values) / len(values) if values else 0.0
        expected_sd = math.sqrt(variance)
        observed_value = observed[motif]
        exceedances = sum(value >= observed_value for value in values)
        p_value = (exceedances + 1) / (len(values) + 1)
        if expected_sd == 0:
            z_score = math.inf if observed_value > expected_mean else 0.0
        else:
            z_score = (observed_value - expected_mean) / expected_sd

        rows.append({
            "motif": motif,
            "length": len(motif),
            "observed": observed_value,
            "expected_mean": round(expected_mean, 4),
            "expected_sd": round(expected_sd, 4),
            "z_score": z_score,
            "p_value": p_value,
        })

    fdr_values = benjamini_hochberg([row["p_value"] for row in rows])
    for row, fdr in zip(rows, fdr_values):
        row["fdr"] = fdr

    return pd.DataFrame(rows).sort_values(["fdr", "p_value", "motif"]).reset_index(drop=True)
