import math
import re
from collections import Counter

import pandas as pd

from .disorder import analyze_idr


POSITIVE = set("KR")
NEGATIVE = set("DE")
CHARGED = POSITIVE | NEGATIVE
AROMATIC = set("FYW")
STICKERS = set("FYWR")
PRION_LIKE = set("QNGSY")
POLAR_SPACERS = set("QNSTG")
RGG_PATTERN = re.compile(r"(?=(RGG|GRG|GGR))")


def _fraction(sequence, residues):
    if not sequence:
        return 0.0
    return sum(aa in residues for aa in sequence) / len(sequence)


def _longest_run(sequence, residues):
    longest = 0
    current = 0
    for aa in sequence:
        if aa in residues:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def _interaction_pair_density(sequence, residues_a, residues_b, max_separation=6):
    if not sequence:
        return 0.0
    count = 0
    for idx, aa in enumerate(sequence):
        if aa not in residues_a:
            continue
        end = min(len(sequence), idx + max_separation + 1)
        for other in sequence[idx + 1:end]:
            if other in residues_b:
                count += 1
    return count / len(sequence)


def _mean_spacing(positions):
    if len(positions) < 2:
        return 0.0
    distances = [right - left for left, right in zip(positions, positions[1:])]
    return sum(distances) / len(distances)


def _charge_value(aa):
    if aa in POSITIVE:
        return 1
    if aa in NEGATIVE:
        return -1
    return 0


def sequence_charge_decoration(sequence):
    """Approximate SCD-style charge-pattern score.

    Negative values indicate nearby opposite-charge patterning; positive values
    indicate like-charge blockiness. This is a compact descriptor, not a full
    polymer-physics model.
    """
    charges = [_charge_value(aa) for aa in sequence]
    n = len(charges)
    if n < 2:
        return 0.0

    score = 0.0
    for i in range(n):
        if charges[i] == 0:
            continue
        for j in range(i + 1, n):
            if charges[j] == 0:
                continue
            score += charges[i] * charges[j] * math.sqrt(j - i)
    return score / n


def charge_blockiness_proxy(sequence, window=10):
    """Returns a local-NCPR variance proxy for charge segregation."""
    if len(sequence) < 2:
        return 0.0
    window = max(2, min(window, len(sequence)))
    local = []
    for start in range(0, len(sequence) - window + 1):
        segment = sequence[start:start + window]
        local.append(sum(_charge_value(aa) for aa in segment) / window)
    if not local:
        return 0.0
    mean = sum(local) / len(local)
    return sum((value - mean) ** 2 for value in local) / len(local)


def rgg_motif_count(sequence):
    return sum(1 for _ in RGG_PATTERN.finditer(sequence))


def llps_summary_features(sequence):
    """Calculates global sequence features associated with IDR-mediated LLPS."""
    sequence = sequence.replace("\n", "").strip().upper()
    length = len(sequence)
    if length == 0:
        return {
            "LLPS_Length": 0,
            "LLPS_HeuristicScore": 0.0,
        }

    counts = Counter(sequence)
    positive = sum(counts[aa] for aa in POSITIVE)
    negative = sum(counts[aa] for aa in NEGATIVE)
    sticker_positions = [idx for idx, aa in enumerate(sequence, start=1) if aa in STICKERS]
    _, _, idr_summary = analyze_idr(sequence)

    aromatic_fraction = _fraction(sequence, AROMATIC)
    sticker_fraction = _fraction(sequence, STICKERS)
    prion_like_fraction = _fraction(sequence, PRION_LIKE)
    polar_spacer_fraction = _fraction(sequence, POLAR_SPACERS)
    rgg_count = rgg_motif_count(sequence)
    cation_pi_density = _interaction_pair_density(sequence, POSITIVE, AROMATIC)
    aromatic_pair_density = _interaction_pair_density(sequence, AROMATIC, AROMATIC)
    fcr = (positive + negative) / length
    ncpr = (positive - negative) / length
    scd = sequence_charge_decoration(sequence)
    blockiness = charge_blockiness_proxy(sequence)

    score = (
        0.28 * idr_summary["DisorderFraction"]
        + 0.18 * idr_summary["LowComplexityFraction"]
        + 0.16 * min(prion_like_fraction / 0.45, 1.0)
        + 0.13 * min(sticker_fraction / 0.18, 1.0)
        + 0.10 * min(cation_pi_density / 0.20, 1.0)
        + 0.08 * min(rgg_count / max(length / 80, 1.0), 1.0)
        + 0.07 * min(fcr / 0.35, 1.0)
    )

    return {
        "LLPS_Length": length,
        "LLPS_HeuristicScore": round(score, 4),
        "LLPS_DisorderFraction": idr_summary["DisorderFraction"],
        "LLPS_LowComplexityFraction": idr_summary["LowComplexityFraction"],
        "LLPS_PrionLikeFraction": round(prion_like_fraction, 4),
        "LLPS_PolarSpacerFraction": round(polar_spacer_fraction, 4),
        "LLPS_StickerFraction": round(sticker_fraction, 4),
        "LLPS_AromaticFraction": round(aromatic_fraction, 4),
        "LLPS_CationPiDensity": round(cation_pi_density, 4),
        "LLPS_AromaticPairDensity": round(aromatic_pair_density, 4),
        "LLPS_FCR": round(fcr, 4),
        "LLPS_NCPR": round(ncpr, 4),
        "LLPS_AbsNCPR": round(abs(ncpr), 4),
        "LLPS_SCD": round(scd, 4),
        "LLPS_ChargeBlockinessProxy": round(blockiness, 4),
        "LLPS_RGGMotifCount": rgg_count,
        "LLPS_RGGDensity": round(rgg_count / length, 4),
        "LLPS_StickerMeanSpacing": round(_mean_spacing(sticker_positions), 4),
        "LLPS_LongestPolyQ": _longest_run(sequence, {"Q"}),
        "LLPS_LongestPolyG": _longest_run(sequence, {"G"}),
        "LLPS_LongestPrionLikeRun": _longest_run(sequence, PRION_LIKE),
    }


def llps_profile(sequence, window=31):
    """Returns sliding-window LLPS-related features."""
    sequence = sequence.replace("\n", "").strip().upper()
    if not sequence:
        return pd.DataFrame(columns=[
            "position", "aa", "local_score", "prion_like_fraction",
            "sticker_fraction", "aromatic_fraction", "fcr", "ncpr",
            "rgg_count",
        ])

    window = max(3, min(window, len(sequence)))
    half = window // 2
    rows = []
    for idx, aa in enumerate(sequence):
        start = max(0, idx - half)
        end = min(len(sequence), idx + half + 1)
        segment = sequence[start:end]
        positive = sum(residue in POSITIVE for residue in segment)
        negative = sum(residue in NEGATIVE for residue in segment)
        prion_like = _fraction(segment, PRION_LIKE)
        stickers = _fraction(segment, STICKERS)
        aromatic = _fraction(segment, AROMATIC)
        fcr = (positive + negative) / len(segment)
        ncpr = (positive - negative) / len(segment)
        rgg = rgg_motif_count(segment)
        local_score = (
            0.35 * min(prion_like / 0.45, 1.0)
            + 0.25 * min(stickers / 0.18, 1.0)
            + 0.15 * min(aromatic / 0.08, 1.0)
            + 0.15 * min(fcr / 0.35, 1.0)
            + 0.10 * min(rgg / max(len(segment) / 80, 1.0), 1.0)
        )
        rows.append({
            "position": idx + 1,
            "aa": aa,
            "local_score": local_score,
            "prion_like_fraction": prion_like,
            "sticker_fraction": stickers,
            "aromatic_fraction": aromatic,
            "fcr": fcr,
            "ncpr": ncpr,
            "rgg_count": rgg,
        })
    return pd.DataFrame(rows)


def _segments_from_scores(profile, threshold, min_region):
    segments = []
    start = None
    for _, row in profile.iterrows():
        flag = row["local_score"] >= threshold
        position = int(row["position"])
        if flag and start is None:
            start = position
        elif not flag and start is not None:
            end = position - 1
            if end - start + 1 >= min_region:
                segments.append((start, end))
            start = None

    if start is not None:
        end = int(profile["position"].iloc[-1])
        if end - start + 1 >= min_region:
            segments.append((start, end))
    return segments


def llps_candidate_regions(profile, threshold=0.6, min_region=20):
    """Calls candidate LLPS-prone windows from a local heuristic score."""
    if profile.empty:
        return pd.DataFrame(columns=["start", "end", "length", "mean_local_score"])

    rows = []
    for start, end in _segments_from_scores(profile, threshold, min_region):
        region = profile.loc[(profile["position"] >= start) & (profile["position"] <= end)]
        rows.append({
            "start": start,
            "end": end,
            "length": end - start + 1,
            "mean_local_score": float(region["local_score"].mean()),
            "mean_prion_like_fraction": float(region["prion_like_fraction"].mean()),
            "mean_sticker_fraction": float(region["sticker_fraction"].mean()),
            "mean_fcr": float(region["fcr"].mean()),
        })
    return pd.DataFrame(rows)


def analyze_llps(sequence, window=31, threshold=0.6, min_region=20):
    """Returns global LLPS features, local profile, and candidate regions."""
    summary = llps_summary_features(sequence)
    profile = llps_profile(sequence, window=window)
    regions = llps_candidate_regions(profile, threshold=threshold, min_region=min_region)
    return summary, profile, regions
