# PBMC3k atlas slice: senescence-adjacent module scores across immune cell types

One public single-cell dataset — 10x Genomics PBMC3k (2,700 peripheral blood
cells) — full Scanpy QC and Leiden clustering, rule-based cell-type
annotation from canonical markers, and two senescence-adjacent gene-module
scores (`p21_like`, `sasp_core`) computed per cell type.

**This is a method demonstration on a healthy-donor reference atlas.**
PBMC3k contains no senescence induction experiment. Module-score differences
across cell types do not identify senescent cells; the repo exists to show a
reproducible, provenance-tracked way to slice one public dataset with one
clear question — the pattern to copy for real senescence/organoid datasets.

## Outputs

- `data/results/cluster_module_scores.tsv` — per-cluster cell type, size,
  mean mito %, mean module scores
- `figures/module_scores_by_cell_type.png` — score distributions by cell type
- `METHODS.md` — dataset URL + sha256, software versions, all QC parameters
- `data/provenance/*_run.json` — full run receipt, including which module
  genes were present vs missing in the data

## Usage

```bash
pip install -r requirements.txt
python scripts/analyze.py --email you@university.edu
python -m unittest discover -s tests -v
```

The dataset (~7 MB) downloads on first run; everything runs on a laptop in
a couple of minutes.

## Scientific honesty rules

- Module genes absent from the dataset are recorded in the provenance
  receipt, not silently dropped.
- The `sasp_core` module is *expected* to be near-zero in healthy PBMC;
  the report says so rather than over-interpreting.
- QC thresholds, normalization, HVG, PCA, neighbors, Leiden resolution,
  and the random seed are all recorded in METHODS.md.
- Swap the dataset + modules for a real senescence study (e.g., irradiated
  fibroblast scRNA-seq) to use this as an actual analysis template.

## License

MIT. The PBMC3k dataset belongs to 10x Genomics and is redistributed by
this project only via download from their servers, under their terms.
