# Publication source tables

This directory contains compact, tab-separated source tables for the eMito
manuscript figures and aggregate workflow summaries. Tables are generated from
the final production run rather than maintained manually, so obsolete results
cannot silently remain in the repository.

Expected generated files:

| File | Contents |
|---|---|
| `EXPORT_MANIFEST.tsv` | File-level inventory and description of the exported tables |
| `pipeline_parameters.tsv` | Probe-generation, optional-mode, and filtering parameters |
| `mode_routing.tsv` | Per-species taxa/group routing decisions and selected accessions after input QC |
| `input_taxonomy_summary.tsv` | Genome, species, genus, and family counts before/after input QC |
| `mode_processing_summary.tsv` | Generation/access/collapse counts for each mode |
| `final_merge_summary.tsv` | Final merge and exact-sequence deduplication counts |
| `probe_stage_summary.tsv` | Aggregate counts across probe-processing stages |
| `probe_terminal_mode_contributions.tsv` | Unique/overlapping taxa, group, and node contributions |
| `taxa_access_filter_summary.tsv` | GC, complexity, dimer, and within-target deduplication losses |
| `panel_A_probe_counts.tsv` | Source data for the normalized probe-filtering panel |
| `panel_B_focal_probe_counts.tsv` | Source data for the focal-taxon bar plot |
| `genus_probe_counts.tsv` | Source data for the circular genus tree |

Generate all tables after running both manuscript plotting scripts:

```bash
python scripts/export_publication_metadata.py \
  --pipeline-root /path/to/completed_emito_run \
  --plot-dir /path/to/generated_plot_tables \
  --output-dir metadata
```

Replace the two `/path/to/...` placeholders with locations on your own system.
The exporter does not modify the pipeline result. Paths recorded in exported
plot tables are stored relative to `--pipeline-root`, so the publication tables
remain portable and do not expose machine-specific directory layouts.
