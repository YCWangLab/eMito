# Script entry points

Install eMito first with `python -m pip install -e .`, then use the `emito`
command documented in the repository README. These thin wrappers are provided
for clusters and legacy workflows that prefer an explicit Python script.

| Script | Equivalent command |
|---|---|
| `run_emito_pipeline.py` | `emito run` |
| `count_emito_mode_merge.py` | `emito summarize` |
| `count_emito_probe_stages.py` | `emito stage-summary` |
| `count_input_taxonomy.py` | `emito taxonomy-summary` |

Every argument is forwarded unchanged, for example:

```bash
python scripts/run_emito_pipeline.py --help
```
