from collections import Counter, defaultdict

import pandas as pd
from Bio import SeqIO


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = set("XBZUO")


def analyze_fasta_qc(fasta_file, min_length=30):
    """Checks FASTA records for common sequence-quality issues."""
    with open(fasta_file, "r", encoding="utf-8") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    ids = [record.id for record in records]
    sequences = [str(record.seq).replace("\n", "").strip().upper() for record in records]
    lengths = [len(seq) for seq in sequences]

    issues = []
    id_counts = Counter(ids)
    sequence_to_ids = defaultdict(list)

    for record_id, seq in zip(ids, sequences):
        sequence_to_ids[seq].append(record_id)
        invalid = sorted(set(seq) - STANDARD_AA - AMBIGUOUS_AA - set("*-"))
        ambiguous = sorted(set(seq) & AMBIGUOUS_AA)

        if id_counts[record_id] > 1:
            issues.append({"record_id": record_id, "issue": "DuplicateID", "detail": f"Seen {id_counts[record_id]} times"})
        if len(seq) < min_length:
            issues.append({"record_id": record_id, "issue": "ShortSequence", "detail": f"Length {len(seq)} < {min_length}"})
        if "*" in seq:
            issues.append({"record_id": record_id, "issue": "StopSymbol", "detail": "Contains '*'"})
        if "-" in seq:
            issues.append({"record_id": record_id, "issue": "GapSymbol", "detail": "Contains '-'"})
        if invalid:
            issues.append({"record_id": record_id, "issue": "InvalidResidues", "detail": ",".join(invalid)})
        if ambiguous:
            issues.append({"record_id": record_id, "issue": "AmbiguousResidues", "detail": ",".join(ambiguous)})

    for seq, duplicate_ids in sequence_to_ids.items():
        if len(duplicate_ids) > 1:
            for record_id in duplicate_ids:
                issues.append({
                    "record_id": record_id,
                    "issue": "DuplicateSequence",
                    "detail": ",".join(duplicate_ids),
                })

    summary = {
        "Records": len(records),
        "UniqueIDs": len(set(ids)),
        "DuplicateIDCount": sum(1 for count in id_counts.values() if count > 1),
        "DuplicateSequenceCount": sum(1 for record_ids in sequence_to_ids.values() if len(record_ids) > 1),
        "MinLength": min(lengths) if lengths else 0,
        "MedianLength": float(pd.Series(lengths).median()) if lengths else 0.0,
        "MaxLength": max(lengths) if lengths else 0,
        "MeanLength": float(pd.Series(lengths).mean()) if lengths else 0.0,
        "IssueCount": len(issues),
    }
    issues_df = pd.DataFrame(issues, columns=["record_id", "issue", "detail"])
    return summary, issues_df
