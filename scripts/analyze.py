#!/usr/bin/env python3
"""Public atlas slice: senescence-adjacent module scores across PBMC3k cell types.

One public dataset (10x Genomics PBMC3k, 2,700 cells), full Scanpy QC and
clustering, rule-based cell-type annotation from canonical markers, and two
senescence-adjacent gene-module scores per cell type.

This is a METHOD DEMONSTRATION on a healthy-donor reference dataset. PBMC3k
contains no senescence induction experiment; module-score differences across
cell types do NOT indicate senescent cells and no such claim is made.

Usage:
  python scripts/analyze.py --email you@university.edu
"""
from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

DATASET_URL = (
    "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/"
    "pbmc3k_filtered_gene_bc_matrices.tar.gz"
)

MARKER_SETS: dict[str, list[str]] = {
    "B cell": ["MS4A1", "CD79A"],
    "T cell": ["CD3E", "CD3D"],
    "NK cell": ["NKG7", "GNLY"],
    "CD14+ monocyte": ["CD14", "LST1"],
    "FCGR3A+ monocyte": ["FCGR3A", "MS4A7"],
    "Dendritic cell": ["FCER1A", "CST3"],
    "Platelet": ["PPBP"],
}

MODULES: dict[str, list[str]] = {
    # p53/p21-style cell-cycle entry module (senescence-adjacent, curated)
    "p21_like": ["CDKN1A", "GADD45A", "MDM2", "BTG2", "CDKN1B"],
    # core inflammatory SASP genes (mostly expected to be near-zero in healthy PBMC)
    "sasp_core": ["IL6", "CXCL8", "MMP1", "SERPINE1", "CXCL1", "IL1A"],
}

QC_MIN_GENES = 200
QC_MAX_GENES = 2500
QC_MAX_MITO_FRACTION = 0.05
QC_MIN_CELLS_PER_GENE = 3


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def download_file(url: str, destination: Path, user_agent: str) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": user_agent})
    with urlopen(request, timeout=300) as response, destination.open("wb") as handle:
        digest = hashlib.sha256()
        while True:
            chunk = response.read(1 << 16)
            if not chunk:
                break
            digest.update(chunk)
            handle.write(chunk)
    return digest.hexdigest()


def extract_tarball(tarball: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tarball) as archive:
        archive.extractall(destination)


def filter_module(module: list[str], available_genes: set[str]) -> tuple[list[str], list[str]]:
    """Split a module into genes present vs missing from the dataset."""
    present = [gene for gene in module if gene in available_genes]
    missing = [gene for gene in module if gene not in available_genes]
    return present, missing


def pick_cell_type(marker_means: dict[str, float]) -> str | None:
    """Assign the cell type with the highest mean marker expression."""
    if not marker_means:
        return None
    return max(marker_means, key=lambda key: marker_means[key])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="analyze")
    parser.add_argument("--email", default="pbmc3k-atlas-slice@local", help="polite-pool contact")
    parser.add_argument("--out", default=".", help="project root")
    parser.add_argument("--seed", type=int, default=0, help="random seed for reproducibility")
    args = parser.parse_args(argv)

    try:
        import anndata
        import matplotlib

        matplotlib.use("Agg")
        import numpy as np
        import scanpy as sc
    except ImportError as exc:
        raise RuntimeError(
            "scanpy/anndata/matplotlib are required: pip install -r requirements.txt"
        ) from exc

    root = Path(args.out).resolve()
    user_agent = f"pbmc3k-atlas-slice/1.0 (mailto:{args.email})"
    receipt: dict = {
        "tool": "pbmc3k-atlas-slice",
        "started": now_utc(),
        "dataset_url": DATASET_URL,
        "versions": {
            "scanpy": sc.__version__,
            "anndata": anndata.__version__,
            "numpy": np.__version__,
        },
        "qc": {
            "min_genes": QC_MIN_GENES,
            "max_genes": QC_MAX_GENES,
            "max_mito_fraction": QC_MAX_MITO_FRACTION,
            "min_cells_per_gene": QC_MIN_CELLS_PER_GENE,
        },
        "modules": MODULES,
        "seed": args.seed,
    }

    # 1. Download + extract the public dataset.
    tarball = root / "data" / "raw" / "pbmc3k_filtered_gene_bc_matrices.tar.gz"
    sha256 = download_file(DATASET_URL, tarball, user_agent)
    receipt["dataset_sha256"] = sha256
    extract_dir = root / "data" / "raw" / "extracted"
    extract_tarball(tarball, extract_dir)
    matrix_candidates = list(extract_dir.rglob("matrix.mtx.gz"))
    if not matrix_candidates:
        matrix_candidates = list(extract_dir.rglob("matrix.mtx"))
    if not matrix_candidates:
        raise FileNotFoundError(f"no 10x matrix found under {extract_dir}")
    matrix_dir = matrix_candidates[0].parent
    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=False)

    # 2. QC.
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[
        (adata.obs["n_genes_by_counts"] >= QC_MIN_GENES)
        & (adata.obs["n_genes_by_counts"] <= QC_MAX_GENES)
        & (adata.obs["pct_counts_mt"] <= QC_MAX_MITO_FRACTION * 100)
    ].copy()
    sc.pp.filter_genes(adata, min_cells=QC_MIN_CELLS_PER_GENE)
    receipt["cells_after_qc"] = int(adata.n_obs)
    receipt["genes_after_qc"] = int(adata.n_vars)

    # 3. Standard preprocessing + clustering.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=40, random_state=args.seed)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=args.seed)
    sc.tl.leiden(adata, resolution=1.0, random_state=args.seed, key_added="cluster")

    # 4. Rule-based cell-type annotation from canonical markers.
    available = set(adata.var_names)
    cluster_types: dict[str, str] = {}
    scaled = adata.to_df()
    for cluster, indices in adata.obs.groupby("cluster").indices.items():
        cluster_scaled = scaled.iloc[indices]
        means = {}
        for cell_type, markers in MARKER_SETS.items():
            present, _ = filter_module(markers, available)
            means[cell_type] = float(cluster_scaled[present].mean(axis=1).mean()) if present else 0.0
        cluster_types[str(cluster)] = pick_cell_type(means)
    adata.obs["cell_type"] = adata.obs["cluster"].astype(str).map(cluster_types)
    receipt["cluster_types"] = cluster_types
    return finish(root, receipt, adata, args.seed)


def finish(root: Path, receipt: dict, adata, seed: int) -> int:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    available = set(adata.var_names)
    module_results: dict[str, str] = {}
    for module_name, genes in MODULES.items():
        present, missing = filter_module(genes, available)
        receipt.setdefault("module_genes", {})[module_name] = {
            "present": present,
            "missing": missing,
        }
        if len(present) < 2:
            print(f"module {module_name}: only {len(present)} genes present; skipped")
            continue
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{module_name}")
        module_results[module_name] = f"score_{module_name}"

    rows = ["cluster\tcell_type\tn_cells\tmean_mito_pct"
            + "".join(f"\tmean_{name}" for name in module_results)]
    for cluster, indices in adata.obs.groupby("cluster").indices.items():
        sub = adata.obs.iloc[indices]
        line = f"{cluster}\t{sub['cell_type'].iloc[0]}\t{len(sub)}\t{round(float(sub['pct_counts_mt'].mean()), 2)}"
        for name in module_results:
            line += f"\t{round(float(sub[module_results[name]].mean()), 4)}"
        rows.append(line)
    results_path = root / "data" / "results" / "cluster_module_scores.tsv"
    results_path.parent.mkdir(parents=True, exist_ok=True)
    results_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    figure_path = root / "figures" / "module_scores_by_cell_type.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    cell_types = sorted(adata.obs["cell_type"].unique())
    n_modules = len(module_results)
    fig, axes = plt.subplots(1, max(n_modules, 1), figsize=(6 * max(n_modules, 1), 5))
    axes_array = [axes] if n_modules <= 1 else list(axes)
    for ax, (name, score_col) in zip(axes_array, module_results.items()):
        data = [
            adata.obs.loc[adata.obs["cell_type"] == cell_type, score_col].values
            for cell_type in cell_types
        ]
        ax.boxplot(data, tick_labels=cell_types, showfliers=False)
        ax.set_title(f"{name} module score by cell type")
        ax.set_ylabel("score_genes (additive)")
        ax.tick_params(axis="x", rotation=45)
    if n_modules == 0:
        axes_array[0].text(0.5, 0.5, "no modules scoreable", ha="center")
    fig.suptitle("PBMC3k reference atlas: senescence-adjacent module scores (method demo)")
    fig.tight_layout()
    fig.savefig(figure_path, dpi=150)
    plt.close(fig)

    receipt["finished"] = now_utc()
    receipt_path = root / "data" / "provenance" / f"{time.time_ns()}_run.json"
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2), encoding="utf-8")

    (root / "METHODS.md").write_text(render_methods(receipt, seed), encoding="utf-8")
    print(
        f"cells: {receipt['cells_after_qc']}, genes: {receipt['genes_after_qc']}, "
        f"clusters: {len(receipt['cluster_types'])}; results: {results_path.name}; "
        f"figure: {figure_path.name}; provenance: {receipt_path.name}"
    )
    return 0


def render_methods(receipt: dict, seed: int) -> str:
    return (
        "# Methods\n\n"
        f"Generated by `pbmc3k-atlas-slice` on {receipt['started']}.\n\n"
        f"- Dataset: 10x Genomics PBMC3k (`{DATASET_URL}`), "
        f"sha256 `{receipt['dataset_sha256']}`.\n"
        f"- Software: scanpy {receipt['versions']['scanpy']}, "
        f"anndata {receipt['versions']['anndata']}, numpy {receipt['versions']['numpy']}.\n"
        f"- QC: {receipt['qc']['min_genes']}-{receipt['qc']['max_genes']} genes/cell, "
        f"mito fraction <= {receipt['qc']['max_mito_fraction']}, genes in >= "
        f"{receipt['qc']['min_cells_per_gene']} cells.\n"
        "- Normalization: per-cell totals to 1e4, log1p; top-2,000 HVGs; scaled to "
        f"max 10; 40-component PCA; 10-NN graph; Leiden (resolution 1.0, seed {seed}).\n"
        "- Cell types: rule-based assignment by maximum mean scaled expression of "
        "canonical marker pairs (see scripts/analyze.py MARKER_SETS).\n"
        "- Module scores: `sc.tl.score_genes` over the genes present in the data "
        "(present/missing lists in data/provenance receipts).\n\n"
        "**Interpretation limits:** PBMC3k is a healthy-donor reference with no\n"
        "senescence induction; module-score variation across cell types is shown as a\n"
        "method demonstration and does not identify senescent cells.\n\n"
        f"Reproduce: `python scripts/analyze.py --email <contact> --seed {seed}`\n"
    )


if __name__ == "__main__":
    raise SystemExit(main())
