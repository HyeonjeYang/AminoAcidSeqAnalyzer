"""Standalone: FASTA quality-control checks."""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.qc import analyze_fasta_qc


def main():
    parser = argparse.ArgumentParser(description="FASTA quality-control checks.")
    parser.add_argument("fasta_file", help="Path to FASTA file")
    parser.add_argument("--min_length", type=int, default=30)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()

    summary, issues = analyze_fasta_qc(args.fasta_file, min_length=args.min_length)
    print("\nFASTA QC summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    print("\nIssues:")
    print(issues if not issues.empty else "  None")
    if args.save:
        issues.to_csv("qc_issues.csv", index=False)


if __name__ == "__main__":
    main()
