import pandas as pd

from .motifs import get_repeat_positions, getrep


def select_top_repeated_motifs(sequence, min_len=3, max_len=10, min_reps=2, top_n=20):
    """Returns the most frequent exact repeated motifs for track plotting."""
    counts = getrep(sequence, min_len=min_len, max_len=max_len, min_reps=min_reps)
    return [
        motif for motif, _ in sorted(counts.items(), key=lambda item: (-item[1], len(item[0]), item[0]))[:top_n]
    ]


def motif_track_dataframe(
    record_id,
    sequence,
    motifs=None,
    min_len=3,
    max_len=10,
    min_reps=2,
    top_n=20,
):
    """Builds a table of motif positions for one sequence."""
    if motifs is None:
        motifs = select_top_repeated_motifs(sequence, min_len, max_len, min_reps, top_n)

    rows = get_repeat_positions(sequence, motifs)
    for row in rows:
        row["record_id"] = record_id
    columns = ["record_id", "motif", "start", "end", "length"]
    return pd.DataFrame(rows, columns=columns)
