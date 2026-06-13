import re

import pandas as pd
from Bio.SeqUtils.ProtParamData import kd


ANNOTATION_COLUMNS = ["type", "start", "end", "motif", "description", "score"]


def _regex_rows(sequence, pattern, annotation_type, description):
    rows = []
    regex = re.compile(f"(?=({pattern}))")
    for match in regex.finditer(sequence):
        motif = match.group(1)
        rows.append({
            "type": annotation_type,
            "start": match.start(1) + 1,
            "end": match.start(1) + len(motif),
            "motif": motif,
            "description": description,
            "score": None,
        })
    return rows


def _merge_mask(mask, min_length):
    regions = []
    start = None
    for idx, flag in enumerate(mask, start=1):
        if flag and start is None:
            start = idx
        elif not flag and start is not None:
            end = idx - 1
            if end - start + 1 >= min_length:
                regions.append((start, end))
            start = None
    if start is not None:
        end = len(mask)
        if end - start + 1 >= min_length:
            regions.append((start, end))
    return regions


def _hydrophobic_average(segment):
    return sum(kd.get(aa, 0.0) for aa in segment) / len(segment)


def annotate_functional_motifs(
    sequence,
    tm_window=19,
    tm_threshold=1.6,
    pest_window=12,
    pest_threshold=0.45,
):
    """Annotates common rule-based protein sequence motifs."""
    sequence = sequence.replace("\n", "").strip().upper()
    rows = []

    rows.extend(_regex_rows(
        sequence, r"N[^P][ST]", "N-glycosylation", "N-X-S/T sequon where X is not proline"
    ))
    rows.extend(_regex_rows(
        sequence, r"[ST]P", "Phosphorylation", "Proline-directed Ser/Thr kinase motif"
    ))
    rows.extend(_regex_rows(
        sequence, r"[ST]..[DE]", "Phosphorylation", "CK2-like acidic Ser/Thr motif"
    ))
    rows.extend(_regex_rows(
        sequence, r"[KR][KR].?[ST]", "Phosphorylation", "Basic kinase Ser/Thr motif"
    ))
    rows.extend(_regex_rows(
        sequence, r"[KR]{4,}", "NLS", "Monopartite basic nuclear localization signal"
    ))
    rows.extend(_regex_rows(
        sequence, r"[KR]{2}.{10,12}[KR]{3,}", "NLS", "Bipartite nuclear localization signal"
    ))

    pest_mask = [False] * len(sequence)
    pest_letters = set("PESTD")
    for start in range(0, max(len(sequence) - pest_window + 1, 0)):
        segment = sequence[start:start + pest_window]
        fraction = sum(aa in pest_letters for aa in segment) / len(segment)
        if fraction >= pest_threshold and "P" in segment and any(aa in segment for aa in "ED"):
            for idx in range(start, start + pest_window):
                pest_mask[idx] = True
    for start, end in _merge_mask(pest_mask, pest_window):
        segment = sequence[start - 1:end]
        score = sum(aa in pest_letters for aa in segment) / len(segment)
        rows.append({
            "type": "PEST",
            "start": start,
            "end": end,
            "motif": segment,
            "description": "PEST-like degradation-prone region",
            "score": score,
        })

    tm_mask = [False] * len(sequence)
    for start in range(0, max(len(sequence) - tm_window + 1, 0)):
        segment = sequence[start:start + tm_window]
        avg = _hydrophobic_average(segment)
        if avg >= tm_threshold:
            for idx in range(start, start + tm_window):
                tm_mask[idx] = True
    for start, end in _merge_mask(tm_mask, tm_window):
        segment = sequence[start - 1:end]
        rows.append({
            "type": "Transmembrane",
            "start": start,
            "end": end,
            "motif": segment,
            "description": "Hydrophobic helix candidate",
            "score": _hydrophobic_average(segment),
        })

    n_term = sequence[:30]
    if len(n_term) >= 15:
        basic_n_region = sum(aa in "KR" for aa in n_term[:8])
        best_core = max(
            (_hydrophobic_average(n_term[i:i + 8]) for i in range(0, len(n_term) - 7)),
            default=0.0,
        )
        if basic_n_region >= 1 and best_core >= 1.6:
            rows.append({
                "type": "SignalPeptideLike",
                "start": 1,
                "end": min(30, len(sequence)),
                "motif": n_term,
                "description": "N-terminal signal peptide-like hydrophobic core",
                "score": best_core,
            })

    if not rows:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)
    return pd.DataFrame(rows, columns=ANNOTATION_COLUMNS).sort_values(["start", "end", "type"]).reset_index(drop=True)
