import re

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
    """Extracts per-residue pLDDT scores (0-100 scale) from an ESMFold PDB
    (stored as a 0-1 fraction in the B-factor column)."""
    scores = []
    for line in pdb_text.splitlines():
        if _PLDDT_PATTERN.match(line):
            scores.append(float(line[60:66]) * 100)
    return scores


def predict_mean_plddt(sequence, timeout=120):
    """Returns the mean pLDDT score for a sequence, or None if the request fails."""
    try:
        pdb_text = predict_structure(sequence, timeout=timeout)
    except requests.RequestException as exc:
        print(f"  ESMFold request failed: {exc}")
        return None

    scores = parse_plddt(pdb_text)
    if not scores:
        return None
    return sum(scores) / len(scores)
