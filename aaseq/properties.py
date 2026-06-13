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


def parse_position_spec(position_spec, sequence_length):
    """Parses 1-based positions such as '3,10-15' into 0-based indexes."""
    if not position_spec:
        return None

    positions = set()
    for token in str(position_spec).split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            start_text, end_text = token.split("-", 1)
            start = int(start_text)
            end = int(end_text)
            for pos in range(start, end + 1):
                if 1 <= pos <= sequence_length:
                    positions.add(pos - 1)
        else:
            pos = int(token)
            if 1 <= pos <= sequence_length:
                positions.add(pos - 1)

    return sorted(positions)


def motif_positions(sequence, motif):
    """Returns 0-based positions covered by all exact motif matches."""
    sequence = sequence.replace("\n", "").strip().upper()
    motif = motif.replace("\n", "").strip().upper()
    positions = set()
    start = 0
    while motif:
        idx = sequence.find(motif, start)
        if idx == -1:
            break
        positions.update(range(idx, idx + len(motif)))
        start = idx + 1
    return sorted(positions)


def generate_point_mutants(sequence, num_mutations=20, seed=None, positions=None, target_motif=None):
    """Generates random single-point amino acid mutants of a sequence.

    Returns (mutant_ids, mutant_sequences), where the first entry is always
    the unmutated 'Wild_Type' sequence.
    """
    sequence = sequence.replace("\n", "").strip().upper()
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    seq_list = list(sequence)
    rng = random.Random(seed)
    mutants = [sequence]
    mutant_ids = ["Wild_Type"]
    if target_motif:
        positions = motif_positions(sequence, target_motif)
    elif isinstance(positions, str):
        positions = parse_position_spec(positions, len(sequence))
    elif positions is not None:
        positions = sorted({int(pos) for pos in positions if 0 <= int(pos) < len(sequence)})
    else:
        positions = list(range(len(sequence)))

    if not positions:
        raise ValueError("No valid mutation positions were provided.")

    for _ in range(num_mutations):
        idx = rng.choice(positions)
        original_aa = seq_list[idx]
        if original_aa not in amino_acids:
            choices = amino_acids
        else:
            choices = amino_acids.replace(original_aa, "")
        new_aa = rng.choice(choices)

        mutated_seq = seq_list.copy()
        mutated_seq[idx] = new_aa
        mutants.append("".join(mutated_seq))
        mutant_ids.append(f"Mut_{original_aa}{idx+1}{new_aa}")

    return mutant_ids, mutants


def generate_systematic_single_mutants(sequence, positions=None, target_motif=None):
    """Generates all single amino-acid substitutions at selected positions."""
    sequence = sequence.replace("\n", "").strip().upper()
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    if target_motif:
        positions = motif_positions(sequence, target_motif)
    elif isinstance(positions, str):
        positions = parse_position_spec(positions, len(sequence))
    elif positions is not None:
        positions = sorted({int(pos) for pos in positions if 0 <= int(pos) < len(sequence)})
    else:
        positions = list(range(len(sequence)))

    if not positions:
        raise ValueError("No valid mutation positions were provided.")

    mutants = [sequence]
    mutant_ids = ["Wild_Type"]
    for idx in positions:
        original_aa = sequence[idx]
        choices = amino_acids.replace(original_aa, "") if original_aa in amino_acids else amino_acids
        for new_aa in choices:
            mutated_seq = list(sequence)
            mutated_seq[idx] = new_aa
            mutants.append("".join(mutated_seq))
            mutant_ids.append(f"Mut_{original_aa}{idx+1}{new_aa}")
    return mutant_ids, mutants
