# AminoAcidSeqAnalyzer

Personal toolkit for analyzing amino acid sequences: find meaningful repeated motifs, compute physicochemical properties, and explore patterns across multiple sequences.

## Requirements

- A `.fasta` file (or a raw sequence string for quick mode)
- Python packages: `biopython`, `pandas`, `numpy`, `scikit-learn`, `scipy`, `seaborn`, `matplotlib`, `tqdm`

Install with:

```bash
pip install -r requirements.txt
```

## Usage

```bash
python AminoAcidAnalyzer.py <fasta_file> [options]
```

### Core options

- `--min_len`, `--max_len`, `--min_reps`: control repeat-motif search range and minimum repetition count
- `--workers`, `--chunk_size`: tune parallel processing for large FASTA files
- `--save`: save the resulting feature matrix to `repeat_matrix.csv`

### Analysis features

- **Repeat motif heatmap**: visualizes exact repeated motifs across sequences
- **Fuzzy repeat motif families** (`--fuzzy_motifs`): clusters same-length windows that differ by a small number of residues (`--fuzzy_mismatches`)
- **Motif significance** (`--motif_stats`): composition-preserving sequence shuffles with z-scores, empirical p-values, and FDR
- **IDR / low-complexity analysis** (`--idr`): TOP-IDP-style disorder propensity, normalized sequence entropy, region calls, and plots
- **LLPS / condensate feature analysis** (`--llps`): heuristic IDR-mediated phase-separation features, including prion-like/polar composition, aromatic/cationic stickers, RGG motifs, and charge-pattern descriptors
- **Motif position tracks** (`--motif_tracks`): maps repeated motifs along the selected sequence
- **Physicochemical properties**: per-sequence composition ratios (non-polar, polar, aromatic, positive, negative), N-glycosylation sequon count (`N-X-[S/T]`), plus ProtParam metrics (molecular weight, isoelectric point, instability index, aromaticity, GRAVY)
- **Functional annotation** (`--annotate`): rule-based N-glycosylation, phosphorylation-like, NLS, PEST-like, transmembrane, and signal-peptide-like motif calls
- **FASTA QC** (`--qc`): duplicate IDs/sequences, invalid/ambiguous residues, stop/gap symbols, and length summaries
- **Group enrichment** (`--group_enrichment --metadata metadata.csv`): tests motif/property enrichment between metadata groups
- **Reference-guided MSA** (`--msa`): dependency-free star MSA around a reference/longest sequence plus conservation profile
- `--pca`: PCA + KMeans clustering and a hierarchical clustermap (`--clusters` sets the number of clusters)
- `--phylogeny`: feature-based phylogenetic dendrogram (clusters sequences by motif/property similarity)
- `--seq_phylogeny`: alignment-based phylogenetic tree across all sequences in the FASTA file, using real pairwise BLOSUM62 alignment distances
- `--mutate`: in silico point-mutagenesis profile heatmap for a chosen sequence (`--num_mutations`, `--seed`, `--mutation_mode`, `--mutation_positions`, and `--target_motif` control the panel)
- `--mutate_phylogeny`: builds an alignment-based phylogenetic tree of the random point mutants, showing how each mutation shifts the sequence relative to the wild type
- `--mutate_plddt`: queries the [ESMFold](https://esmatlas.com/about) API for mean pLDDT (structural confidence, 0-100) of the wild type and mutants. **Requires an internet connection**. For sequences over 400 residues, the tool estimates pLDDT from overlapping <=400aa chunks; this avoids truncating the pLDDT profile but does not produce a full-length global structure.
- `--plddt_profile`: plots per-residue ESMFold pLDDT for a selected sequence, using overlapping chunks for long sequences
- `--hydrophobicity`: Kyte-Doolittle hydrophobicity profile plot for a chosen sequence
- `--record_id`: pick which FASTA record to use for selected-sequence analyses (defaults to the first record)
- `--report`: generate an HTML summary report (`report/summary.html`) bundling all generated tables and figures

### Quick single-sequence mode

Analyze a sequence directly without a FASTA file:

```bash
python AminoAcidAnalyzer.py --seq "MKTAYIAKQRQ..." --hydrophobicity --mutate_phylogeny
```

This prints repeat motifs and physicochemical properties, and can also plot the hydrophobicity profile and a mutagenesis phylogenetic tree.

Examples:

```bash
python AminoAcidAnalyzer.py proteins.fasta --qc --fuzzy_motifs --idr --motif_tracks --annotate --save
python AminoAcidAnalyzer.py proteins.fasta --msa --seq_phylogeny --report
python AminoAcidAnalyzer.py proteins.fasta --metadata metadata.csv --group_enrichment --save
python AminoAcidAnalyzer.py proteins.fasta --record_id NP_000000.1 --motif_stats --num_shuffles 1000 --idr --llps
python AminoAcidAnalyzer.py proteins.fasta --record_id NP_000000.1 --mutate --mutation_mode systematic --target_motif PEST --save
```

## Standalone feature scripts

Every major feature can also be run on its own via the scripts in `scripts/`, without going through the combined CLI:

- `python scripts/motif_analysis.py <fasta_file> [--save]` - repeat-motif/property matrix + heatmap
- `python scripts/fuzzy_motif_analysis.py <fasta_file> [--record_id ID] [--fuzzy] [--stats] [--save]` - exact/fuzzy motifs and shuffle enrichment
- `python scripts/idr_analysis.py <fasta_file> [--record_id ID] [--save]` - IDR and low-complexity analysis
- `python scripts/llps_analysis.py <fasta_file> [--record_id ID] [--save]` - IDR/LLPS-related sequence features and candidate regions
- `python scripts/motif_track_analysis.py <fasta_file> [--record_id ID] [--save]` - repeated motif coordinate tracks
- `python scripts/annotation_analysis.py <fasta_file> [--record_id ID] [--save]` - rule-based functional motif annotation
- `python scripts/qc_analysis.py <fasta_file> [--save]` - FASTA quality checks
- `python scripts/msa_conservation_analysis.py <fasta_file> [--save]` - reference-guided MSA and conservation profile
- `python scripts/group_enrichment_analysis.py repeat_matrix.csv metadata.csv [--save]` - group-wise feature enrichment
- `python scripts/plddt_profile_analysis.py <fasta_file> [--record_id ID] [--save]` - per-residue ESMFold pLDDT profile
- `python scripts/cluster_analysis.py <fasta_file> [--clusters N] [--save]` - PCA + KMeans clustering and clustermap
- `python scripts/phylogeny_analysis.py <fasta_file> [--feature] [--sequence]` - feature-based and/or alignment-based phylogenetic tree
- `python scripts/mutagenesis_analysis.py <fasta_file> [--record_id ID] [--heatmap] [--phylogeny] [--plddt]` - in silico mutagenesis (heatmap / phylogeny / ESMFold pLDDT), also accepts `--seq` instead of a FASTA file
- `python scripts/hydrophobicity_analysis.py <fasta_file> [--record_id ID] [--window N]` - Kyte-Doolittle hydrophobicity profile, also accepts `--seq`
- `python scripts/quick_seq_analysis.py "MKT..."` - repeat motifs + physicochemical properties for a single sequence

### Optional external tools

These wrappers require separately installed command-line tools. No source code from these projects is copied into this repository.

- `python scripts/mmseqs_cluster_analysis.py proteins.fasta --min_seq_id 0.5 --coverage 0.8` - clusters sequences with an installed [MMseqs2](https://github.com/soedinglab/MMseqs2) binary
- `python scripts/foldseek_search_analysis.py query_structures target_database --output_tsv foldseek_results.tsv` - searches structures with an installed [Foldseek](https://github.com/steineggerlab/foldseek) binary

## Notebooks

The original `AminoAcidAnalyzer.ipynb` is left unchanged. Additional lightweight notebooks are in `notebooks/`:

- `notebooks/IDR_Motif_Workflow.ipynb` - selected-sequence IDR, repeat motif, motif-track, and annotation workflow
- `notebooks/Batch_MSA_Group_Workflow.ipynb` - FASTA QC, matrix creation, MSA conservation, and optional metadata group enrichment

## Long-sequence pLDDT chunking

The public ESMFold endpoint can reject long sequences. For `--plddt_profile`, this project can split long proteins into overlapping <=400aa chunks, run ESMFold per chunk, and merge residue-level pLDDT with center-weighted averaging. This is reasonable for exploratory local confidence profiling, especially around IDRs and low-complexity regions, because pLDDT is a per-residue confidence estimate and overlap reduces boundary artifacts.

It is not a replacement for a full-length structure prediction. Cutting a protein removes long-range domain/domain and terminal context, so chunked pLDDT should not be used to infer global folds, interfaces, allostery, domain packing, or precise effects of mutations outside the chunk context.

## Third-party acknowledgements

No third-party GitHub source code is vendored or copied into this repository. The analyzer uses external Python libraries and public web services through their normal APIs. Those projects retain their own copyright and license terms:

- [Biopython](https://github.com/biopython/biopython/blob/master/LICENSE.rst) - Biopython License Agreement / BSD 3-Clause; used for FASTA parsing, ProtParam metrics, hydrophobicity scale data, and alignment utilities.
- [NumPy](https://github.com/numpy/numpy/blob/main/LICENSE.txt), [pandas](https://github.com/pandas-dev/pandas/blob/main/LICENSE), [scikit-learn](https://github.com/scikit-learn/scikit-learn/blob/main/COPYING), [SciPy](https://github.com/scipy/scipy/blob/main/LICENSE.txt), and [seaborn](https://github.com/mwaskom/seaborn/blob/master/LICENSE.md) - BSD-style licenses; used for numerical analysis, tables, clustering, statistics, and plotting.
- [Matplotlib](https://github.com/matplotlib/matplotlib/blob/main/LICENSE/LICENSE) - Matplotlib license; used for plots.
- [tqdm](https://github.com/tqdm/tqdm/blob/master/LICENCE) - MIT license with noted MPL-2.0 portions; used for progress bars.
- [requests](https://github.com/psf/requests/blob/main/LICENSE) - Apache License 2.0; used for HTTP requests.
- [ESMFold / ESM Atlas](https://esmatlas.com/about) is accessed as a remote service for structure/pLDDT prediction. No Meta ESM source code is bundled here; the upstream [ESM repository](https://github.com/facebookresearch/esm/blob/main/LICENSE) is MIT licensed.
- [MMseqs2](https://github.com/soedinglab/MMseqs2/blob/master/LICENSE.md), associated with the Söding/Steinegger ecosystem, is MIT licensed. This project does not copy MMseqs2 code; it only provides an optional wrapper for a user-installed binary.
- [Foldseek](https://github.com/steineggerlab/foldseek/blob/master/LICENSE.md) and [HH-suite](https://github.com/soedinglab/hh-suite/blob/master/LICENSE) are GPL-3.0 licensed upstream. To avoid mixing GPL source into this MIT project, no Foldseek or HH-suite source code is copied here; optional use should remain through separately installed command-line tools.
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/) were used to download example protein sequences. Follow NCBI's usage guidelines, Disclaimer, and Copyright notice when using or redistributing data obtained from NCBI.

Scientific methods and scales implemented or accessed here build on prior work, including ProtParam, Kyte-Doolittle hydrophobicity, BLOSUM62 alignment scoring, TOP-IDP intrinsic-disorder propensity, charge-pattern descriptors for IDRs, sticker-spacer thinking for condensates, prion-like composition, and pi/pi or cation/pi interaction concepts used in LLPS studies. The LLPS module is an original heuristic feature calculator, not a copied implementation of PScore, PLAAC, FuzDrop, PSPredictor, ParSe, or another trained predictor. Please cite the relevant original scientific publications when using those analyses in publications or reports.

## Project layout

- `AminoAcidAnalyzer.py` - combined CLI entry point and orchestration
- `scripts/` - standalone CLIs for each individual feature (see above)
- `aaseq/motifs.py` - repeat motif detection
- `aaseq/statistics.py` - shuffle-based motif enrichment
- `aaseq/disorder.py` - IDR and low-complexity heuristics
- `aaseq/llps.py` - IDR/LLPS-related sequence feature heuristics
- `aaseq/tracks.py` - motif coordinate tables
- `aaseq/annotations.py` - rule-based functional motif annotations
- `aaseq/alignment.py` - reference-guided MSA and conservation
- `aaseq/enrichment.py` - metadata group enrichment
- `aaseq/qc.py` - FASTA quality control
- `aaseq/external_tools.py` - optional wrappers for separately installed external tools
- `aaseq/properties.py` - physicochemical properties, hydrophobicity profile, mutant generation
- `aaseq/processing.py` - FASTA parsing and parallel feature-matrix construction
- `aaseq/visualization.py` - plots: heatmaps, PCA/clustering, dendrograms, hydrophobicity, ESMFold pLDDT comparison
- `aaseq/structure.py` - ESMFold API client (structure/pLDDT prediction)
- `aaseq/report.py` - HTML report generation and helper lookups
