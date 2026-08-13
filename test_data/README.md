# Test data

The automated test suite builds a small synthetic taxonomy, metadata table,
mitochondrial FASTA collection, and a deterministic MAFFT stand-in inside a
temporary directory. This keeps the repository small while exercising the
complete pipeline, including:

- species-level and subgroup-level alignment;
- circular-origin and reverse-complement normalization;
- short-genome exclusion and non-ATCG window removal;
- species/genus and subgroup specificity;
- node selection with `NC_` preference;
- access, collapse, and optional final sequence deduplication; and
- all three summary commands.

Run it with:

```bash
python -m unittest discover -s tests -v
```
