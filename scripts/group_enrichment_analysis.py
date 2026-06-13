"""Standalone: group-wise motif/property enrichment from a feature matrix and metadata CSV."""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.enrichment import group_feature_enrichment


def main():
    parser = argparse.ArgumentParser(description="Group-wise feature enrichment.")
    parser.add_argument("matrix_csv", help="Feature matrix CSV, indexed by sequence ID")
    parser.add_argument("metadata_csv", help="Metadata CSV with sequence IDs and groups")
    parser.add_argument("--id_col", default="id")
    parser.add_argument("--group_col", default="group")
    parser.add_argument("--min_present", type=int, default=1)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix_csv, index_col=0)
    result = group_feature_enrichment(matrix, args.metadata_csv, args.id_col, args.group_col, min_present=args.min_present)
    print(result.head(30))
    if args.save:
        result.to_csv("group_enrichment.csv", index=False)


if __name__ == "__main__":
    main()
