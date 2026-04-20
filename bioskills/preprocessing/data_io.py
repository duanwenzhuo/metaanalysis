"""
L0: DataIOSkill — scRNA-seq / Bulk RNA-seq 数据加载

bioSkills: bio-sc-data-io → bioskills L0
支持格式:
- 10X CellRanger (mtx + barcodes + features)
- 10X H5
- H5AD (AnnData)
- Loom
- CSV / TSV
- Geomtx (gzip compressed)

bioSkills Method Comparison:
| 格式         | 工具          | 适用场景                     |
|-------------|--------------|------------------------------|
| 10X MTX     | scanpy       | CellRanger 输出              |
| 10X H5      | scanpy       | CellRanger 输出(压缩)        |
| H5AD        | anndata      | AnnData 标准                 |
| Loom        | loompy       | 旧格式兼容                   |
| CSV/TSV     | pandas       | 简单矩阵                     |
"""

import os
from typing import Optional
from bioskills.core.base import AbstractSkill


class DataIOSkill(AbstractSkill):
    name = "data_io"
    stage = "preprocessing"
    input_contract = ["input_path"]
    output_contract = ["adata", "metadata", "raw_counts"]

    tunables = {
        "format": "auto",            # auto | 10x_mtx | 10x_h5 | h5ad | loom | csv | tsv
        "genome": None,              # for 10X H5 (e.g. "GRCh38")
        "sheet": None,               # for Excel files
        "sep": None,                 # for CSV/TSV override
        "sparse": True,              # load as sparse if possible
        "cache": True,               # cache loaded data
        "backup_url": None,          # download if not found locally
        "min_genes": 0,              # basic cell filter on load
        "min_cells": 0,              # basic gene filter on load
        "obs_columns": None,         # specific metadata columns to load
        "var_columns": None,         # specific gene metadata columns
    }

    FORMAT_DETECTORS = {
        ".h5ad": "h5ad",
        ".h5": "10x_h5",
        ".loom": "loom",
        ".csv": "csv",
        ".tsv": "tsv",
        ".mtx": "10x_mtx",
        ".mtx.gz": "10x_mtx",
    }

    def _detect_format(self, path: str) -> str:
        if path.endswith(".mtx.gz") or path.endswith(".mtx"):
            return "10x_mtx"
        for ext, fmt in self.FORMAT_DETECTORS.items():
            if path.endswith(ext):
                return fmt
        # Check if directory (10X format)
        if os.path.isdir(path):
            if os.path.exists(os.path.join(path, "matrix.mtx")) or \
               os.path.exists(os.path.join(path, "matrix.mtx.gz")):
                return "10x_mtx"
            if any(f.endswith(".h5") for f in os.listdir(path)):
                return "10x_h5"
        return "csv"  # fallback

    def _run(self, state: dict) -> dict:
        params = {**self.tunables, **state.get("params", {})}
        input_path = state.get("input_path")

        if input_path is None:
            return {"error": "input_path required in state"}

        fmt = params.get("format", "auto")
        if fmt == "auto":
            fmt = self._detect_format(input_path)

        try:
            import scanpy as sc
            import anndata as ad
            import numpy as np
            from scipy import sparse
        except ImportError:
            return {"error": "scanpy/anndata required: pip install scanpy anndata"}

        try:
            adata = None
            loader_report = {"format": fmt, "path": input_path, "status": "success"}

            if fmt == "10x_mtx":
                genome = params.get("genome")
                adata = sc.read_10x_mtx(input_path, var_names="gene_symbols",
                                        cache=params.get("cache", True),
                                        gex_only=True)
                loader_report["genome"] = genome or "auto"

            elif fmt == "10x_h5":
                genome = params.get("genome", None)
                kwargs = {"genome": genome} if genome else {}
                adata = sc.read_10x_h5(input_path, **kwargs)
                # Convert var names to gene symbols if available
                if "gene_ids" in adata.var.columns:
                    adata.var_names_make_unique()

            elif fmt == "h5ad":
                adata = sc.read_h5ad(input_path)

            elif fmt == "loom":
                adata = sc.read_loom(input_path)

            elif fmt in ("csv", "tsv"):
                import pandas as pd
                sep = params.get("sep") or ("," if fmt == "csv" else "\t")
                df = pd.read_csv(input_path, sep=sep, index_col=0)
                # Rows = genes, Cols = cells (standard for scRNA-seq)
                adata = ad.AnnData(X=sparse.csr_matrix(df.values) if params.get("sparse", True) else df.values)
                adata.var_names = df.index.astype(str)
                adata.obs_names = df.columns.astype(str)

            else:
                return {"error": f"Unsupported format: {fmt}"}

            if adata is None:
                return {"error": f"Failed to load data from {input_path}"}

            # Basic filtering on load
            min_genes = params.get("min_genes", 0)
            min_cells = params.get("min_cells", 0)
            if min_genes > 0:
                sc.pp.filter_cells(adata, min_genes=min_genes)
            if min_cells > 0:
                sc.pp.filter_genes(adata, min_cells=min_cells)

            # Make var_names unique
            adata.var_names_make_unique()

            # Report
            loader_report.update({
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "sparse": sparse.issparse(adata.X),
                "layers": list(adata.layers.keys()) if adata.layers else [],
                "obs_keys": list(adata.obs.columns),
                "var_keys": list(adata.var.columns),
            })

            return {
                "status": "success",
                "state_updates": {
                    "adata": adata,
                    "metadata": adata.obs.to_dict(orient="records") if len(adata.obs.columns) > 0 else [],
                    "raw_counts": adata.X,
                    "data_io_report": loader_report,
                }
            }

        except Exception as e:
            return {
                "error": f"DataIO failed: {str(e)}",
                "state_updates": {
                    "data_io_report": {
                        "format": fmt,
                        "path": input_path,
                        "status": "error",
                        "error": str(e),
                    }
                }
            }
