import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

from .statistics import benjamini_hochberg


def load_metadata(metadata):
    if isinstance(metadata, pd.DataFrame):
        return metadata.copy()
    return pd.read_csv(metadata)


def group_feature_enrichment(
    matrix,
    metadata,
    id_col="id",
    group_col="group",
    feature_cols=None,
    min_present=1,
):
    """Tests feature enrichment for each group against all other sequences.

    Binary presence is evaluated with Fisher's exact test. Count/intensity shifts
    are also summarized with a Mann-Whitney U p-value when possible.
    """
    meta = load_metadata(metadata)
    if id_col not in meta.columns:
        raise ValueError(f"Metadata is missing id column '{id_col}'.")
    if group_col not in meta.columns:
        raise ValueError(f"Metadata is missing group column '{group_col}'.")

    df = matrix.copy()
    df.index = df.index.astype(str)
    meta[id_col] = meta[id_col].astype(str)
    merged = df.merge(meta[[id_col, group_col]], left_index=True, right_on=id_col, how="inner")
    if merged.empty:
        raise ValueError("No matrix rows matched metadata IDs.")

    if feature_cols is None:
        feature_cols = [
            col for col in matrix.columns
            if col != "Cluster" and pd.api.types.is_numeric_dtype(matrix[col])
        ]

    rows = []
    for group_name in sorted(merged[group_col].dropna().unique()):
        in_group = merged[group_col] == group_name
        out_group = ~in_group
        for feature in feature_cols:
            in_values = merged.loc[in_group, feature]
            out_values = merged.loc[out_group, feature]
            if in_values.empty or out_values.empty:
                continue

            a = int((in_values > 0).sum())
            b = int((in_values <= 0).sum())
            c = int((out_values > 0).sum())
            d = int((out_values <= 0).sum())
            if a + c < min_present:
                continue

            odds_ratio, fisher_p = fisher_exact([[a, b], [c, d]], alternative="greater")
            try:
                _, mannwhitney_p = mannwhitneyu(in_values, out_values, alternative="two-sided")
            except ValueError:
                mannwhitney_p = 1.0

            rows.append({
                "group": group_name,
                "feature": feature,
                "group_present": a,
                "group_absent": b,
                "other_present": c,
                "other_absent": d,
                "group_mean": float(in_values.mean()),
                "other_mean": float(out_values.mean()),
                "odds_ratio": float(odds_ratio),
                "fisher_p": float(fisher_p),
                "mannwhitney_p": float(mannwhitney_p),
            })

    if not rows:
        return pd.DataFrame(columns=[
            "group", "feature", "group_present", "group_absent", "other_present",
            "other_absent", "group_mean", "other_mean", "odds_ratio",
            "fisher_p", "fisher_fdr", "mannwhitney_p",
        ])

    result = pd.DataFrame(rows)
    result["fisher_fdr"] = benjamini_hochberg(result["fisher_p"].tolist())
    return result.sort_values(["fisher_fdr", "fisher_p", "group", "feature"]).reset_index(drop=True)
