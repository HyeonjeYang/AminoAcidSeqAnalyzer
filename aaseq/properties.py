import random
import re

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd

# N-glycosylation sequon: N-X-[S/T], where X is any residue except Proline
NGLYCO_PATTERN = re.compile(r"N[^P][ST]")


def calculate_aa_properties(sequence):
    """Calculates the composition ratio (%) of chemical properties of amino acids in the sequence,
    plus ProtParam-based physicochemical metrics and N-glycosylation sequon count."""
    groups = {
        "Non-polar": set("AILMFVPG"),
        "Polar": set("STCNQY"),
        "Aromatic": set("FYW"),
        "Positive": set("KRH"),
        "Negative": set("DE")
    }
    seq_len = len(sequence) if len(sequence) > 0 else 1
    props = {}
    for prop_name, aa_set in groups.items():
        count = sum(1 for aa in sequence if aa in aa_set)
        props[f"Prop_{prop_name}"] = round((count / seq_len) * 100, 2)

    # N-glycosylation sequon (N-X-S/T) count
    props["NGlyco_Sites"] = len(NGLYCO_PATTERN.findall(sequence))

    # ProtParam-based physicochemical properties (skipped if sequence contains
    # non-standard residues that ProteinAnalysis cannot handle)
    try:
        analysis = ProteinAnalysis(sequence)
        props["MolecularWeight"] = round(analysis.molecular_weight(), 2)
        props["IsoelectricPoint"] = round(analysis.isoelectric_point(), 2)
        props["InstabilityIndex"] = round(analysis.instability_index(), 2)
        props["Aromaticity"] = round(analysis.aromaticity(), 4)
        props["GRAVY"] = round(analysis.gravy(), 4)
    except (ValueError, KeyError, ZeroDivisionError):
        pass
    return props


def get_hydrophobicity_profile(sequence, window=9):
    """Returns (positions, scores) for a Kyte-Doolittle hydrophobicity profile."""
    sequence = sequence.replace("\n", "").strip()
    analysis = ProteinAnalysis(sequence)
    scores = analysis.protein_scale(kd, window=window, edge=1.0)
    offset = window // 2
    positions = list(range(offset + 1, offset + 1 + len(scores)))
    return positions, scores


def generate_point_mutants(sequence, num_mutations=20):
    """Generates random single-point amino acid mutants of a sequence.

    Returns (mutant_ids, mutant_sequences), where the first entry is always
    the unmutated 'Wild_Type' sequence.
    """
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    seq_list = list(sequence)
    mutants = [sequence]
    mutant_ids = ["Wild_Type"]

    for _ in range(num_mutations):
        idx = random.randint(0, len(seq_list) - 1)
        original_aa = seq_list[idx]
        new_aa = random.choice(amino_acids.replace(original_aa, ""))

        mutated_seq = seq_list.copy()
        mutated_seq[idx] = new_aa
        mutants.append("".join(mutated_seq))
        mutant_ids.append(f"Mut_{original_aa}{idx+1}{new_aa}")

    return mutant_ids, mutants
