# eMito

**eMito** designs mitochondrial capture probes using taxonomy-aware sequence comparisons and reference tiling. Six independent modules let users choose how to prepare genomes, generate probes, filter them, reduce overlap and combine sets.

| Module | Function |
| --- | --- |
| `eMito-prepare` | Validate inputs, filter genomes by length, normalise circular sequences, align and sample probe/k-mer windows |
| `eMito-taxa-generate` | Generate species-, genus- and subgroup-specific candidates |
| `eMito-node-generate` | Tile a selected reference for each requested taxonomic node |
| `eMito-access` | Filter each input set by GC content, DUST complexity and complementary-sequence score |
| `eMito-collapse` | Retain non-overlapping probes from each source genome within an input set |
| `eMito-merge` | Combine sets and remove exact duplicate sequences |

Each module runs only the requested operation. Assessment, collapse and merging are optional and can be selected and ordered to suit an experiment.

## Installation

Requires Python 3.9 or newer. Install [MAFFT](https://mafft.cbrc.jp/alignment/software/) separately for the default alignment-enabled preparation.

```bash
git clone https://github.com/YCWangLab/eMito.git
cd eMito
python3 -m pip install .
emito --help
emito prepare --help
```

The installed commands `emito prepare` and `eMito-prepare` are equivalent; the same convention applies to the other five modules. Source-checkout wrappers also work without installing the package:

```bash
python3 scripts/eMito-prepare.py --help
```

The Python implementation uses the standard library. MAFFT is not required when preparation uses `--no-align` or when processing existing probe sets.

## Inputs

Provide an accession FASTA directory, NCBI `nodes.dmp` and `names.dmp` from the same taxonomy release, and a tab-separated metadata file with these four columns:

```text
accession_id	species_name	taxid	subgroup_label
NC_012920.1	Homo sapiens	9606	
```

Use one ungapped mitochondrial genome per FASTA, named `<accession_id>.fasta`, including the accession version. `species_name` must match the scientific name of the species ancestor of `taxid`. Leave `subgroup_label` empty for ordinary species records; use it to distinguish subspecies, populations or other intended within-species groups. See [examples/metadata.tsv](examples/metadata.tsv). Input accessions must be unique.

## Standard panel

This example assesses species/subgroup candidates, supplements them with family-reference tiling and merges the results. It does not collapse either branch.

```bash
emito prepare --metadata metadata.tsv --fasta-dir fasta \
  --names-dmp taxdump/names.dmp --nodes-dmp taxdump/nodes.dmp \
  --min-genome-length 10000 --threads 4 --output results/01_prepare

emito taxa-generate --prepared results/01_prepare --output results/02_taxa_generate
emito access --inputs results/02_taxa_generate --output results/03_taxa_access
emito node-generate --prepared results/01_prepare --node-rank family \
  --output results/04_node_generate
emito merge --inputs results/03_taxa_access results/04_node_generate \
  --output results/05_final_merge
```

The final FASTA is `results/05_final_merge/probe_sets/merged.fasta`. Counts are recorded in `merge_summary.json`. Every module requires a **new output directory**; completed outputs are not overwritten or resumed. A successful stage writes `COMPLETE`.

For this same combination, a configurable launcher is provided:

```bash
export FASTA_DIR=/absolute/path/fasta
export NAMES_DMP=/absolute/path/taxdump/names.dmp
export NODES_DMP=/absolute/path/taxdump/nodes.dmp
MIN_GENOME_LENGTH=8000 bash run_standard_panel.sh \
  /absolute/path/metadata.tsv /absolute/path/new_results
```

`MIN_GENOME_LENGTH` defaults to **10000 bp**, and accepts any positive integer. The CLI equivalent is `--min-genome-length 8000`. `PIPELINE_PYTHON`, `MAFFT_BIN` and `THREADS` can override executable selection and worker count. The launcher requires absolute input/output paths.

For SLURM, submit from the repository root after exporting the input variables:

```bash
sbatch submit_standard_panel.sh /absolute/path/metadata.tsv /absolute/path/new_results
```

Use `EMITO_CODE_DIR=/absolute/path/eMito` if submitting from another directory. Set partition, memory and time through your site's `sbatch` options. The example requests four CPUs and 100 GB; adjust these to the dataset and cluster.

## Preparation and specificity

`prepare` filters genomes below the requested minimum length. Within each species pool or subgroup, it selects a reference by preferring `NC_` accessions, then the lowest proportion of non-ATCG bases, with accession order breaking ties. It normalises strand and circular origin to that reference and aligns pools containing multiple genomes with MAFFT `--auto`.

Candidate probes are sampled at fixed intervals: **52 bp long, every 5 reference bases** by default. Reference positions are mapped through the alignment to each genome. A gap at the start skips that window; otherwise the program extracts consecutive ungapped bases. Windows crossing the end continue from the beginning of the circular genome. Windows containing non-ATCG bases are removed. Comparison k-mers use the same length with a default **1 bp** step and include reverse complements; candidate probes do not receive an additional reverse-complement copy. `--window-length`, `--probe-step` and `--kmer-step` are adjustable. `--no-align` samples the raw circular genomes and skips both alignment and strand/origin normalisation.

`taxa-generate` builds representative k-mer sets from sequences present in all genomes of a target with up to three genomes, or at least `ceil(0.75 × n)` genomes for larger targets. Configure this using `--small-group-all-max` and `--representative-fraction`. Species-specific sequences are absent from the representative sets of other species. Genus-restricted sequences can be shared within a genus but are absent from other genera. Each species' candidates are intersected with its eligible sequences and deduplicated within the target.

Species analysis preferentially uses records without a subgroup assignment. If none are available and only one taxonomic subgroup represents a species, that subgroup supplies the species analysis. Multiple subgroups are compared separately; explicit populations defined at the species TaxID remain subgroup targets. When direct species records and one taxonomic subgroup coexist, the subgroup is checked against both other subgroup representatives and other species' representatives. These decisions are saved in the preparation output. Subgroup generation is part of `taxa-generate`; there is no separate group-generation command.

`node-generate` uses all genomes that passed preparation QC, independently of their species/subgroup role. It selects one reference per family by default, using the same reference preference, and tiles 52 bp probes at 5 bp intervals. Choose another rank with `--node-rank` or supply `--selection selection.tsv` with `node_taxid` and `accession_id` columns. Every node must have a manual selection unless `--allow-auto-unlisted` is supplied. Node windows use the sequence coordinate system saved during preparation.

## Filtering, collapse and merging

`access`, `collapse` and `merge` accept completed output bundles or individual probe FASTAs through `--inputs`.

- **Assessment:** default GC and DUST ranges are 35–65% and 0–2. Dimer scores count reverse-complement 11-mer matches in the pool remaining after these filters, exclude self-contribution and normalise by window count and pool size. The default `--dimer 0.15` uses the score at zero-based index `floor(0.15 × N)` as its cutoff, retaining scores at or below it, including ties. Values at least 1 specify an absolute cutoff; values at most 0 disable this filter. Exact duplicates are removed within each assessed input set.
- **Collapse:** probes are processed by genomic position and retained only when they do not overlap an earlier retained probe from the same accession, including across the circular origin. Different genomes' coordinates are not collapsed together. Input headers must retain position information; conflicting coordinate frames cause an error.
- **Merge:** exact sequence deduplication is enabled by default. `--no-dedup` disables deduplication at this step. Identical strings retain the first source header; reverse-complement strings are not considered identical.

For example, create a lower-density alternative without changing the assessed input:

```bash
emito collapse --inputs results/03_taxa_access --output results/03_taxa_collapsed
emito merge --inputs results/03_taxa_collapsed results/04_node_generate \
  --output results/06_collapsed_panel
```

The chosen order matters: merging before assessment changes the Dimer scoring pool, while deduplication before collapse discards alternative source coordinates.

## Matched comparisons with and without collapse

The comparison runner creates species-only panels with both branches: `access → merge` and `access → collapse → merge`. It excludes subgroup and node supplements from the comparison.

```bash
python3 scripts/run_matched_comparison.py \
  --metadata reference_metadata.tsv --fasta-dir reference_fasta \
  --names-dmp taxdump/names.dmp --nodes-dmp taxdump/nodes.dmp \
  --min-genome-length 10000 --threads 4 --output comparisons/reference

python3 scripts/run_matched_comparison.py \
  --existing-standard results --output comparisons/emito
```

For an already selected one-reference-per-species dataset, add `--no-align --require-single-reference`. Use an explicit smaller `--min-genome-length` if short references should be included. Reusing a standard panel preserves its original preparation and length filter. The runner writes `comparison_summary.tsv`, separate FASTAs in `04_no_collapse_merge` and `06_collapse_merge`, and `COMPLETE` after both branches finish. `submit_matched_panel.sh` accepts the same arguments for SLURM. `scripts/audit_prepare_refseq.py --help` describes the optional offline reference audit.

## Outputs and migration

Preparation records input/QC decisions, normalized genomes, alignments, windows and `prepared.json`. Probe bundles contain `manifest.tsv`, `bundle.json`, FASTAs and `COMPLETE`. The manifest maps targets to FASTAs and counts. Filtering and collapse also write `processing_summary.tsv`; merge writes `merge_inputs.tsv` and `merge_summary.json`. Keep stage directories together: prepared data and some comparison views refer to upstream files.

The modular commands replace the previous `emito run` workflow. Run old commands explicitly as `emito legacy run`, `emito legacy summarize`, `emito legacy stage-summary`, `emito legacy taxonomy-summary` or `emito legacy access-summary` when reproducing or inspecting the old output layout. Existing reporting wrappers and publication exporters target that legacy layout; they do not summarize modular output bundles. Use the new manifests and summaries for current runs. The earlier implementation and history remain available for reproducibility.

## Testing and citation

```bash
python3 -m pip install -e '.[test]'
python3 -m unittest discover -s tests -v
```

Tests cover circular windows and collapse, ambiguous-base handling, routing, subgroup backgrounds, independent modules and both matched-comparison branches. Synthetic alignment tests mock MAFFT; they do not replace validation on a production dataset. A GitHub Actions template is available in `examples/github-actions-tests.yml`.

See [CITATION.cff](CITATION.cff) for citation metadata and [LICENSE](LICENSE) for the MIT licence.
