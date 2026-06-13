"""Standalone: optional MMseqs2 sequence clustering wrapper.

Requires a separately installed `mmseqs` executable. No MMseqs2 source code is
bundled in this repository.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.external_tools import ExternalToolError, run_mmseqs_easy_cluster


def main():
    parser = argparse.ArgumentParser(description="Cluster protein sequences with an installed MMseqs2 binary.")
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("--output_prefix", default="mmseqs_clusters")
    parser.add_argument("--tmp_dir", default="mmseqs_tmp")
    parser.add_argument("--min_seq_id", type=float, default=None)
    parser.add_argument("--coverage", type=float, default=None)
    parser.add_argument("--sensitivity", type=float, default=None)
    parser.add_argument("--mmseqs_bin", default="mmseqs")
    args = parser.parse_args()

    try:
        outputs = run_mmseqs_easy_cluster(
            args.fasta_file,
            args.output_prefix,
            args.tmp_dir,
            min_seq_id=args.min_seq_id,
            coverage=args.coverage,
            sensitivity=args.sensitivity,
            mmseqs_bin=args.mmseqs_bin,
        )
    except ExternalToolError as exc:
        parser.error(str(exc))

    print("MMseqs2 clustering outputs:")
    for label, path in outputs.items():
        print(f"  {label}: {path}")


if __name__ == "__main__":
    main()
