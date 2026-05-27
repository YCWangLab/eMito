# eMito probe-design scripts

This directory contains the reproducible probe-design workflow for the eMito mitochondrial capture probe set.

The scripts are intentionally written as standalone command-line programs so that the repository remains easy to audit and easy to run on an HPC server without installing a custom Python package.

## Design logic

The design workflow follows the eMito manuscript logic:

1. prepare a curated mitochondrial reference panel;
2. generate 52 bp k-mers from all reference mitogenomes using a 1 bp step;
3. generate 52 bp probe tiles from representative mitogenomes using a 5 bp step;
4. identify species- and genus-specific k-mers by cross-comparison across the panel;
5. retain candidate probes containing informative k-mers and filter them by GC content, sequence complexity, and self-dimerization potential;
6. reduce redundancy by grouping overlapping probes and keeping the central probe from each overlap cluster;
7. add supplementary unfiltered tiling probes at family and haplogroup levels;
8. merge and deduplicate all probe components.

## Scripts

```text
01_prepare_reference_panel.py
02_generate_52bp_kmers.py
03_generate_52bp_probe_tiles.py
04_identify_species_genus_specific_kmers.py
05_filter_probe_gc_complexity_dimer.py
06_reduce_overlapping_probes.py
07_generate_haplogroup_family_tiling.py
08_merge_deduplicate_probes.py
```

## Recommended input metadata

The core input is a tab-separated metadata file describing the mitochondrial reference panel.

Recommended columns:

```text
sequence_id
accession
fasta_path
family
genus
species
group
haplogroup
population
is_representative
used_for_kmer_screening
used_for_probe_tiling
used_for_family_tiling
used_for_haplogroup_tiling
```

Minimal columns:

```text
sequence_id
fasta_path
family
genus
species
```

Notes:

- `sequence_id` should be unique.
- `fasta_path` can be an absolute path or a path relative to `--fasta-dir`.
- `used_for_probe_tiling` should mark one high-quality representative mitogenome per target taxon.
- `used_for_family_tiling` can mark one representative species/genome per family.
- `used_for_haplogroup_tiling` can mark one mitochondrial genome per Homo/domestic haplogroup or lineage.
- If representative flags are missing, the tiling scripts can fall back to using the first sequence per family or haplogroup.

## Example workflow

The example below assumes this project structure:

```text
metadata/reference_panel_raw.tsv
data/reference_mitogenomes/
results/design/
probes/
```

Create output directories:

```bash
mkdir -p results/design/{01_reference_panel,02_kmers,03_probe_tiles,04_informative_kmers,05_filtered_probes,06_overlap_reduced,07_supplementary_tiling,08_final}
mkdir -p probes
```

### 1. Prepare reference panel

```bash
python scripts/design/01_prepare_reference_panel.py \
  --metadata metadata/reference_panel_raw.tsv \
  --fasta-dir data/reference_mitogenomes \
  --out-fasta results/design/01_reference_panel/reference_panel.clean.fasta \
  --out-metadata results/design/01_reference_panel/reference_panel.clean.tsv
```

### 2. Generate 52 bp k-mers

```bash
python scripts/design/02_generate_52bp_kmers.py \
  --fasta results/design/01_reference_panel/reference_panel.clean.fasta \
  --metadata results/design/01_reference_panel/reference_panel.clean.tsv \
  --k 52 \
  --step 1 \
  --out results/design/02_kmers/reference_panel.52bp_kmers.tsv.gz
```

### 3. Generate 52 bp probe tiles

```bash
python scripts/design/03_generate_52bp_probe_tiles.py \
  --fasta results/design/01_reference_panel/reference_panel.clean.fasta \
  --metadata results/design/01_reference_panel/reference_panel.clean.tsv \
  --representative-column used_for_probe_tiling \
  --probe-len 52 \
  --step 5 \
  --out-fasta results/design/03_probe_tiles/representative_tiles.fasta \
  --out-metadata results/design/03_probe_tiles/representative_tiles.metadata.tsv
```

### 4. Identify species/genus-specific k-mers

```bash
python scripts/design/04_identify_species_genus_specific_kmers.py \
  --kmers results/design/02_kmers/reference_panel.52bp_kmers.tsv.gz \
  --out results/design/04_informative_kmers/species_genus_specific_kmers.tsv.gz \
  --summary results/design/04_informative_kmers/species_genus_specific_kmers.summary.tsv
```

### 5. Filter candidate probes

Default manuscript-like thresholds:

- GC content: 40–60%
- DUST-style complexity score: < 2
- self-dimer score: <= 8

```bash
python scripts/design/05_filter_probe_gc_complexity_dimer.py \
  --probes-fasta results/design/03_probe_tiles/representative_tiles.fasta \
  --probes-metadata results/design/03_probe_tiles/representative_tiles.metadata.tsv \
  --informative-kmers results/design/04_informative_kmers/species_genus_specific_kmers.tsv.gz \
  --gc-min 40 \
  --gc-max 60 \
  --dust-max 2 \
  --dimer-max 8 \
  --out-fasta results/design/05_filtered_probes/filtered_informative_probes.fasta \
  --out-metadata results/design/05_filtered_probes/filtered_informative_probes.metadata.tsv \
  --rejected results/design/05_filtered_probes/rejected_probes.tsv
```

### 6. Reduce overlapping probes

```bash
python scripts/design/06_reduce_overlapping_probes.py \
  --probes-fasta results/design/05_filtered_probes/filtered_informative_probes.fasta \
  --metadata results/design/05_filtered_probes/filtered_informative_probes.metadata.tsv \
  --out-fasta results/design/06_overlap_reduced/informative_overlap_reduced.fasta \
  --out-metadata results/design/06_overlap_reduced/informative_overlap_reduced.metadata.tsv \
  --cluster-summary results/design/06_overlap_reduced/overlap_clusters.tsv
```

### 7. Generate family/haplogroup supplementary tiling probes

```bash
python scripts/design/07_generate_haplogroup_family_tiling.py \
  --fasta results/design/01_reference_panel/reference_panel.clean.fasta \
  --metadata results/design/01_reference_panel/reference_panel.clean.tsv \
  --levels family,haplogroup \
  --probe-len 52 \
  --step 5 \
  --out-fasta results/design/07_supplementary_tiling/family_haplogroup_tiling.fasta \
  --out-metadata results/design/07_supplementary_tiling/family_haplogroup_tiling.metadata.tsv
```

### 8. Merge and deduplicate final probe set

```bash
python scripts/design/08_merge_deduplicate_probes.py \
  --fasta \
    results/design/06_overlap_reduced/informative_overlap_reduced.fasta \
    results/design/07_supplementary_tiling/family_haplogroup_tiling.fasta \
  --metadata \
    results/design/06_overlap_reduced/informative_overlap_reduced.metadata.tsv \
    results/design/07_supplementary_tiling/family_haplogroup_tiling.metadata.tsv \
  --dedup-mode exact \
  --out-fasta probes/emito_probes.fasta.gz \
  --out-metadata probes/emito_probes.metadata.tsv.gz \
  --summary results/design/08_final/final_probe_summary.tsv
```

## Important implementation notes

### DUST score

`05_filter_probe_gc_complexity_dimer.py` implements a lightweight DUST-style triplet complexity score. If the final manuscript used a specific external e-probe or DUST implementation, the README and methods should report the exact implementation. The function in this repository is designed to make the workflow reproducible and transparent.

### Self-dimer score

The self-dimer score is an approximate maximum contiguous complementary run between a probe and its reverse complement. It is intended as a simple, auditable filter for obvious high-risk probes.

### Reverse complements

By default, probe generation emits the forward tile sequence only. Use `--include-rc` in `03_generate_52bp_probe_tiles.py` if reverse-complement tiles are explicitly required.

### Deduplication

`08_merge_deduplicate_probes.py` supports two modes:

```text
exact         collapse exact duplicate sequences only
canonical_rc  collapse exact duplicates and reverse-complement pairs
```

For most capture-probe manufacturing tables, `exact` is usually the safer default.

## Suggested output files for GitHub release

```text
probes/emito_probes.fasta.gz
probes/emito_probes.metadata.tsv.gz
results/design/08_final/final_probe_summary.tsv
```

## Minimal software requirements

The scripts use only the Python standard library.

Recommended:

```text
Python >= 3.8
```
