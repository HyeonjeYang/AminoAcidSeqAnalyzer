"""Standalone: optional Foldseek structure-search wrapper.

Requires a separately installed `foldseek` executable. No Foldseek source code
is bundled in this repository.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.external_tools import ExternalToolError, run_foldseek_easy_search


def main():
    parser = argparse.ArgumentParser(description="Search structures with an installed Foldseek binary.")
    parser.add_argument("query_structures", help="Query structure file or directory")
    parser.add_argument("target_database", help="Foldseek target database or structure directory")
    parser.add_argument("--output_tsv", default="foldseek_results.tsv")
    parser.add_argument("--tmp_dir", default="foldseek_tmp")
    parser.add_argument("--foldseek_bin", default="foldseek")
    args = parser.parse_args()

    try:
        output = run_foldseek_easy_search(
            args.query_structures,
            args.target_database,
            args.output_tsv,
            args.tmp_dir,
            foldseek_bin=args.foldseek_bin,
        )
    except ExternalToolError as exc:
        parser.error(str(exc))

    print(f"Foldseek results: {output}")


if __name__ == "__main__":
    main()
