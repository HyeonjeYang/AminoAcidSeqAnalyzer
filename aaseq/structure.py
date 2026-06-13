import re

import numpy as np
import requests

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_MAX_LEN = 400

_PLDDT_PATTERN = re.compile(r"^ATOM\s+\d+\s+CA\s")


def predict_structure(sequence, timeout=120):
    """Submits a sequence to the ESMFold API and returns the predicted PDB text.

    Raises requests.RequestException on network/HTTP errors.
    Sequences longer than ESMFOLD_MAX_LEN are truncated (with a printed warning),
    since the public ESMFold API rejects longer inputs.
    """
    sequence = sequence.replace("\n", "").strip().upper()
    if len(sequence) > ESMFOLD_MAX_LEN:
        print(f"  Warning: sequence length {len(sequence)} exceeds ESMFold API limit "
              f"({ESMFOLD_MAX_LEN}aa). Truncating to the first {ESMFOLD_MAX_LEN} residues.")
        sequence = sequence[:ESMFOLD_MAX_LEN]

    response = requests.post(ESMFOLD_URL, data=sequence, timeout=timeout)
    response.raise_for_status()
    return response.text


def parse_plddt(pdb_text):
    """Extracts per-residue pLDDT scores on a 0-100 scale from an ESMFold PDB."""
    scores = []
    for line in pdb_text.splitlines():
        if _PLDDT_PATTERN.match(line):
            scores.append(float(line[60:66]))
    if scores and max(scores) <= 1.5:
        scores = [score * 100 for score in scores]
    return scores


def predict_mean_plddt(sequence, timeout=120):
    """Returns the mean pLDDT score for a sequence, or None if the request fails.

    For sequences longer than the public API limit, this uses overlapping chunks
    and averages residue-level pLDDT over chunk coverage. This bypasses the mean
    pLDDT length limit, but it does not produce a full-length global structure.
    """
    scores = predict_plddt_profile(sequence, timeout=timeout)
    if scores is None or len(scores) == 0:
        return None
    return float(np.nanmean(scores))


def chunk_sequence(sequence, max_len=ESMFOLD_MAX_LEN, overlap=50):
    """Returns (start, chunk) pairs with 0-based starts."""
    sequence = sequence.replace("\n", "").strip().upper()
    if len(sequence) <= max_len:
        return [(0, sequence)]
    if overlap >= max_len:
        raise ValueError("overlap must be smaller than max_len.")

    chunks = []
    step = max_len - overlap
    start = 0
    while start < len(sequence):
        end = min(start + max_len, len(sequence))
        chunks.append((start, sequence[start:end]))
        if end == len(sequence):
            break
        start += step
    return chunks


def _chunk_weights(length, edge_fraction=0.15):
    """Center-weight chunk residue scores to reduce boundary artifacts."""
    if length <= 2:
        return np.ones(length, dtype=float)

    edge = max(1, int(length * edge_fraction))
    weights = np.ones(length, dtype=float)
    for idx in range(edge):
        weight = (idx + 1) / edge
        weights[idx] = min(weights[idx], weight)
        weights[length - idx - 1] = min(weights[length - idx - 1], weight)
    return weights


def predict_plddt_profile(sequence, timeout=120, max_len=ESMFOLD_MAX_LEN, overlap=50, center_weight=True):
    """Returns per-residue pLDDT scores, chunking long sequences when needed.

    Chunking is an approximation for local pLDDT only. It should not be treated
    as a full-length structure prediction because inter-domain context is lost.
    Overlap plus center-weighted merging reduces edge artifacts.
    """
    sequence = sequence.replace("\n", "").strip().upper()
    if not sequence:
        return None

    if len(sequence) <= max_len:
        try:
            pdb_text = predict_structure(sequence, timeout=timeout)
        except requests.RequestException as exc:
            print(f"  ESMFold request failed: {exc}")
            return None
        scores = parse_plddt(pdb_text)
        return np.array(scores, dtype=float) if scores else None

    print(
        f"  Sequence length {len(sequence)} exceeds {max_len}aa. "
        f"Using overlapping {max_len}aa chunks for pLDDT profile."
    )
    totals = np.zeros(len(sequence), dtype=float)
    counts = np.zeros(len(sequence), dtype=float)

    for start, chunk in chunk_sequence(sequence, max_len=max_len, overlap=overlap):
        try:
            pdb_text = predict_structure(chunk, timeout=timeout)
        except requests.RequestException as exc:
            print(f"  ESMFold chunk {start + 1}-{start + len(chunk)} failed: {exc}")
            continue
        scores = parse_plddt(pdb_text)
        usable = min(len(scores), len(chunk))
        if usable == 0:
            continue
        chunk_scores = np.array(scores[:usable], dtype=float)
        weights = _chunk_weights(usable) if center_weight else np.ones(usable, dtype=float)
        totals[start:start + usable] += chunk_scores * weights
        counts[start:start + usable] += weights

    if not np.any(counts):
        return None
    profile = np.full(len(sequence), np.nan, dtype=float)
    covered = counts > 0
    profile[covered] = totals[covered] / counts[covered]
    return profile
