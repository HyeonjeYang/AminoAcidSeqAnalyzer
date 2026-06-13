import os
import tempfile
import unittest

import pandas as pd

from aaseq.alignment import conservation_dataframe, reference_msa
from aaseq.annotations import annotate_functional_motifs
from aaseq.disorder import analyze_idr
from aaseq.enrichment import group_feature_enrichment
from aaseq.motifs import get_fuzzy_repeats, get_repeat_positions
from aaseq.properties import generate_point_mutants, generate_systematic_single_mutants
from aaseq.qc import analyze_fasta_qc
from aaseq.statistics import motif_enrichment
from aaseq.tracks import motif_track_dataframe


class AdvancedFeatureTests(unittest.TestCase):
    def test_fuzzy_repeats_and_positions(self):
        sequence = "AAAATAAA"
        fuzzy = get_fuzzy_repeats(sequence, min_len=4, max_len=4, min_reps=2, max_mismatches=1)
        self.assertTrue(any(key.startswith("Fuzzy_AAAA") for key in fuzzy))

        positions = get_repeat_positions(sequence, ["AAA"])
        self.assertEqual(positions[0]["start"], 1)
        self.assertGreaterEqual(len(positions), 2)

    def test_motif_enrichment(self):
        stats = motif_enrichment("AAAAAKKKKAAAAA", min_len=3, max_len=3, min_reps=2, num_shuffles=20, seed=1)
        self.assertIn("fdr", stats.columns)
        self.assertFalse(stats.empty)

    def test_idr_and_tracks(self):
        sequence = "M" + "P" * 20 + "EEEEKKKKSSSS" + "L" * 10
        scores, regions, summary = analyze_idr(sequence, min_region=5)
        self.assertEqual(summary["Length"], len(sequence))
        self.assertGreater(summary["IDR_Count"], 0)
        self.assertEqual(len(scores), len(sequence))

        tracks = motif_track_dataframe("seq1", "AAAAAKAAAAA", min_len=3, max_len=4, min_reps=2)
        self.assertFalse(tracks.empty)

    def test_annotations(self):
        annotations = annotate_functional_motifs("MKKKKAAAAAANSTPESPESTEDLLLLLLLLLLLLLLLLLLLL")
        self.assertIn("N-glycosylation", set(annotations["type"]))
        self.assertIn("NLS", set(annotations["type"]))

    def test_mutant_generation(self):
        ids_a, mutants_a = generate_point_mutants("ACDE", num_mutations=3, seed=7)
        ids_b, mutants_b = generate_point_mutants("ACDE", num_mutations=3, seed=7)
        self.assertEqual(ids_a, ids_b)
        self.assertEqual(mutants_a, mutants_b)

        ids, mutants = generate_systematic_single_mutants("ACDE", positions="2")
        self.assertEqual(len(ids), 20)
        self.assertEqual(len(mutants), 20)

    def test_reference_msa_conservation(self):
        alignment = reference_msa(["a", "b", "c"], ["MKTAA", "MKTA", "MRTAA"])
        conservation = conservation_dataframe(alignment)
        self.assertEqual(len(alignment), 3)
        self.assertFalse(conservation.empty)
        self.assertLessEqual(conservation["conservation"].max(), 1.0)

    def test_group_enrichment(self):
        matrix = pd.DataFrame(
            {"AAA": [3, 2, 0, 0], "Prop_Polar": [10.0, 11.0, 30.0, 31.0]},
            index=["s1", "s2", "s3", "s4"],
        )
        metadata = pd.DataFrame({"id": ["s1", "s2", "s3", "s4"], "group": ["A", "A", "B", "B"]})
        result = group_feature_enrichment(matrix, metadata)
        self.assertFalse(result.empty)
        self.assertIn("fisher_fdr", result.columns)

    def test_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, "w", encoding="utf-8") as handle:
                handle.write(">seq1\nACDEFGHIK\n>seq1\nACDZX*\n>seq3\nACDEFGHIK\n")
            summary, issues = analyze_fasta_qc(fasta_path, min_length=10)
        self.assertEqual(summary["Records"], 3)
        self.assertGreater(summary["IssueCount"], 0)
        self.assertIn("DuplicateID", set(issues["issue"]))


if __name__ == "__main__":
    unittest.main()
