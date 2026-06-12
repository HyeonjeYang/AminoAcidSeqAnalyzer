# AminoAcidSeqAnalyzer

Personal toolkit for analyzing amino acid sequences: find meaningful repeated motifs, compute physicochemical properties, and explore patterns across multiple sequences.

## Requirements

- A `.fasta` file (or a raw sequence string for quick mode)
- Python packages: `biopython`, `pandas`, `numpy`, `scikit-learn`, `scipy`, `seaborn`, `matplotlib`, `tqdm`

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
- **Physicochemical properties**: per-sequence composition ratios (non-polar, polar, aromatic, positive, negative), N-glycosylation sequon count (`N-X-[S/T]`), plus ProtParam metrics (molecular weight, isoelectric point, instability index, aromaticity, GRAVY)
- `--pca`: PCA + KMeans clustering and a hierarchical clustermap (`--clusters` sets the number of clusters)
- `--phylogeny`: feature-based phylogenetic dendrogram (clusters sequences by motif/property similarity)
- `--seq_phylogeny`: alignment-based phylogenetic tree across all sequences in the FASTA file, using real pairwise BLOSUM62 alignment distances
- `--mutate`: in silico random point-mutagenesis profile heatmap for a chosen sequence (`--num_mutations` sets the mutation count)
- `--mutate_phylogeny`: builds an alignment-based phylogenetic tree of the random point mutants, showing how each mutation shifts the sequence relative to the wild type
- `--hydrophobicity`: Kyte-Doolittle hydrophobicity profile plot for a chosen sequence
- `--record_id`: pick which FASTA record to use for `--mutate` / `--mutate_phylogeny` / `--hydrophobicity` (defaults to the first record)
- `--report`: generate an HTML summary report (`report/summary.html`) bundling all generated tables and figures

### Quick single-sequence mode

Analyze a sequence directly without a FASTA file:

```bash
python AminoAcidAnalyzer.py --seq "MKTAYIAKQRQ..." --hydrophobicity --mutate_phylogeny
```

This prints repeat motifs and physicochemical properties, and can also plot the hydrophobicity profile and a mutagenesis phylogenetic tree.

## Project layout

- `AminoAcidAnalyzer.py` - CLI entry point and orchestration
- `aaseq/motifs.py` - repeat motif detection
- `aaseq/properties.py` - physicochemical properties, hydrophobicity profile, mutant generation
- `aaseq/processing.py` - FASTA parsing and parallel feature-matrix construction
- `aaseq/visualization.py` - plots: heatmaps, PCA/clustering, dendrograms, hydrophobicity
- `aaseq/report.py` - HTML report generation and helper lookups
