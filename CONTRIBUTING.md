# Contributing

Bug reports and pull requests are welcome. Please open an issue before a large
change so that input semantics and probe-design behavior can be discussed.

For local development:

```bash
git clone https://github.com/YCWangLab/eMito.git
cd eMito
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
pytest
```

Changes to taxonomy grouping, representative k-mer thresholds, specificity
backgrounds, coordinate handling, or probe filtering must include a focused
test describing the intended biological behavior.
