"""Standalone: plot exact repeated motif positions along one sequence."""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from aaseq.report import find_record
from aaseq.tracks import motif_track_dataframe
from aaseq.visualization import plot_motif_tracks


def main():
    parser = argparse.ArgumentParser(description="Repeated motif position track analysis.")
    parser.add_argument("fasta_file", nargs="?", help="Path to FASTA file. Not required with --seq.")
    parser.add_argument("--seq", type=str, help="Directly analyze one amino acid sequence")
    parser.add_argument("--record_id", type=str, default=None)
    parser.add_argument("--min_len", type=int, default=3)
    parser.add_argument("--max_len", type=int, default=10)
    parser.add_argument("--min_reps", type=int, default=2)
    parser.add_argument("--top_n", type=int, default=20)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()

    if args.seq:
        sequence = args.seq.replace("\n", "").strip().upper()
        label = "Input_Sequence"
    elif args.fasta_file:
        record = find_record(args.fasta_file, args.record_id)
        if record is None:
            parser.error(f"Record '{args.record_id}' not found in '{args.fasta_file}'.")
        sequence = str(record.seq)
        label = record.id
    else:
        parser.error("Either fasta_file or --seq must be provided.")

    tracks = motif_track_dataframe(label, sequence, min_len=args.min_len, max_len=args.max_len, min_reps=args.min_reps, top_n=args.top_n)
    print(tracks if not tracks.empty else "No repeated motif tracks found.")
    if args.save:
        tracks.to_csv(f"{label}_motif_tracks.csv", index=False)
    plot_motif_tracks(tracks, len(sequence), label)


if __name__ == "__main__":
    main()
