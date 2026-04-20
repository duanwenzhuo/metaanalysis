"""
L0: QC Skill — 质量控制

输入: adata
输出: adata (filtered), n_cells_before, n_cells_after, n_genes_after,
      pct_filtered, pct_mito

契约:
  Input:  ["adata"]
  Output: ["adata", "qc_report"]
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
import numpy as np

@register
class QCSkill(AbstractSkill):
    """单细胞 RNA-seq 质量控制"""
    
    name = "qc"
    description = "Quality control: filter cells/genes by min_genes, min_cells, max_mito_pct"
    input_contract = ["adata"]
    output_contract = ["adata", "qc_report"]
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.GENERIC]
    
    tunable_parameters = {
        "min_genes": {"type": "int", "default": 200, "min": 50, "max": 1000},
        "min_cells": {"type": "int", "default": 3, "min": 1, "max": 50},
        "max_mito_pct": {"type": "float", "default": 0.20, "min": 0.05, "max": 0.50},
        "mito_genes_regex": {"type": "str", "default": "^MT-"},
        "max_genes": {"type": "int", "default": 6000, "min": 500, "max": 50000},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        
        adata = state["adata"]
        params = state.get("params", {})
        
        p = {
            "min_genes": params.get("min_genes", 200),
            "min_cells": params.get("min_cells", 3),
            "max_mito_pct": params.get("max_mito_pct", 0.20),
            "mito_genes_regex": params.get("mito_genes_regex", "^MT-"),
            "max_genes": params.get("max_genes", 6000),
        }
        
        adata = adata.copy()
        n_before = adata.n_obs
        g_before = adata.n_vars
        
        # 线粒体基因比例
        mito_genes = adata.var_names.str.match(p["mito_genes_regex"], na=False)
        if mito_genes.any():
            adata.obs["pct_mito"] = np.asarray(
                adata[:, mito_genes].X.sum(axis=1)
            ).flatten() / np.asarray(adata.X.sum(axis=1)).clip(min=1e-9)
        else:
            adata.obs["pct_mito"] = 0.0
        
        # 过滤
        sc.pp.filter_cells(adata, min_genes=p["min_genes"])
        sc.pp.filter_genes(adata, min_cells=p["min_cells"])
        adata = adata[adata.obs["pct_mito"] < p["max_mito_pct"]].copy()
        adata = adata[adata.obs.n_genes < p["max_genes"]].copy()
        
        n_after = adata.n_obs
        g_after = adata.n_vars
        
        qc_report = {
            "n_cells_before": int(n_before),
            "n_cells_after": int(n_after),
            "n_genes_before": int(g_before),
            "n_genes_after": int(g_after),
            "cells_filtered_pct": round((n_before - n_after) / n_before, 4) if n_before else 0,
            "genes_filtered_pct": round((g_before - g_after) / g_before, 4) if g_before else 0,
            "mean_mito_pct": round(float(adata.obs["pct_mito"].mean()), 4),
            "max_mito_pct_threshold": p["max_mito_pct"],
            "params": p,
        }
        
        self._log(state, f"  QC: {n_before}→{n_after} cells, {g_before}→{g_after} genes")
        
        return {"adata": adata, "qc_report": qc_report}
    
    def _log(self, state: State, msg: str):
        print(msg)
