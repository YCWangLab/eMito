# Scripts

Run these from a source checkout with Python 3.9 or newer; installation is optional.

- `eMito-prepare.py`, `eMito-taxa-generate.py`, `eMito-node-generate.py`, `eMito-access.py`, `eMito-collapse.py`, `eMito-merge.py`: independent current modules; each accepts `--help`.
- `run_emito_pipeline.py`: equivalent to `emito`; pass a module name followed by its options.
- `run_matched_comparison.py`: species-only comparison with and without collapse, using fresh inputs or a completed standard panel.
- `audit_prepare_refseq.py`: offline audit/reference preparation; see `--help` for required inputs.

The `count_*` scripts and `export_publication_metadata.py` support the earlier integrated pipeline output layout only. They require an installed package and should be used for legacy outputs. Current modules write their own manifests and summary tables, as described in the root README.
