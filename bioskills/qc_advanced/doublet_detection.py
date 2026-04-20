"""
L0: DoubletDetection Skill — 双峰检测（Scrublet / scDblFinder）

bioSkills 原版: bio-single-cell-doublet-detection
翻译重构: 适配 bioskills State Dict 契约系统

功能：
- Scrublet（Python）：计算 doublet scores，预测 doublets
- scDblFinder（R）：通过 R 接口调用（仅当adata有R伴侣时）

输出:
  - adata.obs['doublet_score']: float array
  - adata.obs['is_doublet']: bool array
  - doublet_report: 报告
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State


@register
class DoubletDetectionSkill(AbstractSkill):
    """
    Detect and remove doublets from scRNA-seq data.
    
    工具选择策略（bioSkills method comparison）：
    | Method       | Speed | Accuracy | Language |
    |-------------|-------|----------|----------|
    | Scrublet    | Fast  | Good     | Python   |
    | scDblFinder | Fast  | Excellent| R        |
    | DoubletFinder| Slow | Good     | R        |
    
    预期 doublet rates（bioSkills）：
    | Cells Loaded | Expected Rate |
    |-------------|--------------|
    | 1,000       | ~0.8%        |
    | 5,000       | ~4.0%        |
    | 10,000      | ~8.0%        |
    
    公式: rate ≈ cells_loaded / 1000 * 0.008
    """
    
    name = "doublet_detection"
    description = (
        "Detect doublets using Scrublet (Python) or scDblFinder (R). "
        "Remove droplets containing multiple cells before clustering. "
        "bioSkills: bio-single-cell-doublet-detection"
    )
    input_contract = ["adata"]
    output_contract = ["adata", "doublet_report"]
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.SPATIAL]
    
    tunable_parameters = {
        "method": {"type": "str", "default": "scrublet"},
        "expected_doublet_rate": {"type": "float", "default": 0.06},
        "threshold": {"type": "float", "default": None},  # None = auto
        "n_prin_comps": {"type": "int", "default": 30},
        "min_counts": {"type": "int", "default": 2},
        "min_cells": {"type": "int", "default": 3},
        "min_gene_variability_pctl": {"type": "int", "default": 85},
        "filter_doublets": {"type": "bool", "default": False},
    }
    
    def _run(self, state: State) -> dict:
        import numpy as np
        
        adata = state["adata"]
        params = state.get("params", {})
        method = params.get("method", "scrublet")
        
        if method == "scrublet":
            return self._run_scrublet(state)
        elif method == "scdblfinder":
            return self._run_scdblfinder(state)
        else:
            return {
                "adata": adata,
                "doublet_report": {
                    "status": "failed",
                    "error": f"Unknown method: {method}. Use 'scrublet' or 'scdblfinder'.",
                }
            }
    
    def _run_scrublet(self, state: State) -> dict:
        """Scrublet: Fast Python doublet detection"""
        import scanpy as sc
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        # ── Step 1: 检查是否有 scrublet ──────────
        try:
            import scrublet as scr
        except ImportError:
            print("  [DoubletDetection] ⚠️  scrublet not installed. Installing...")
            import subprocess
            subprocess.run(["pip3", "install", "scrublet", "-q"], check=True)
            import scrublet as scr
        
        # ── Step 2: 计算 doublet scores ──────────
        expected_rate = params.get("expected_doublet_rate", 0.06)
        n_cells = adata.n_obs
        
        # bioSkills 建议：按细胞数自动计算 expected_doublet_rate
        auto_rate = n_cells / 1000 * 0.008
        if params.get("expected_doublet_rate") is None:
            expected_rate = auto_rate
            print(f"  [Scrublet] Auto doublet rate: {expected_rate:.3f} for {n_cells} cells")
        
        scrub = scr.Scrublet(
            adata.X,
            expected_doublet_rate=expected_rate,
        )
        
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=params.get("min_counts", 2),
            min_cells=params.get("min_cells", 3),
            min_gene_variability_pctl=params.get("min_gene_variability_pctl", 85),
            n_prin_comps=min(params.get("n_prin_comps", 30), min(adata.n_obs, adata.n_vars)-1),
        )
        
        adata.obs["doublet_score"] = doublet_scores
        adata.obs["is_doublet"] = predicted_doublets
        
        n_doublets = int(predicted_doublets.sum())
        pct_doublets = n_doublets / adata.n_obs
        
        print(f"  [Scrublet] ✅ Detected {n_doublets}/{adata.n_obs} doublets "
              f"({pct_doublets:.1%})")
        
        # ── Step 3: 可视化 ──────────────────────
        try:
            # 简单直方图
            threshold = scrub.threshold_ if hasattr(scrub, 'threshold_') else None
            if threshold is None:
                threshold = params.get("threshold")
            
            # 设置默认 threshold
            if threshold is None:
                threshold = np.median(doublet_scores[doublet_scores > np.percentile(doublet_scores, 50)])
            
            adata.obs["doublet_threshold"] = threshold
            
        except Exception as e:
            print(f"  [Scrublet] ⚠️  Visualization skipped: {e}")
        
        # ── Step 4: 可选过滤 ───────────────────
        filter_doublets = params.get("filter_doublets", False)
        if filter_doublets:
            adata_filt = adata[~adata.obs["is_doublet"]].copy()
            print(f"  [Scrublet] ✅ Filtered → {adata_filt.n_obs} cells remain")
            adata = adata_filt
        else:
            print(f"  [Scrublet] ℹ️  Doublets marked but NOT filtered "
                  f"(set filter_doublets=True to remove)")
        
        report = {
            "method": "scrublet",
            "n_cells_input": int(state["adata"].n_obs),
            "n_cells_output": int(adata.n_obs),
            "n_doublets": n_doublets,
            "pct_doublets": round(pct_doublets, 4),
            "expected_doublet_rate": expected_rate,
            "doublet_score_key": "doublet_score",
            "is_doublet_key": "is_doublet",
            "status": "success",
        }
        
        return {"adata": adata, "doublet_report": report}
    
    def _run_scdblfinder(self, state: State) -> dict:
        """scDblFinder: R/Bioconductor gradient-boosted classifier"""
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {})
        n_cells = adata.n_obs
        
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            from rpy2.robjects.packages import importr
            pandas2ri.activate()
        except ImportError:
            return {
                "adata": adata,
                "doublet_report": {
                    "status": "skipped",
                    "reason": "rpy2 not installed. Use scrublet (Python) instead.",
                }
            }
        
        try:
            scDblFinder = importr("scDblFinder")
            SCE = importr("SingleCellExperiment")
        except Exception:
            return {
                "adata": adata,
                "doublet_report": {
                    "status": "failed",
                    "error": "scDblFinder R package not installed. "
                              "Install with: R -e 'BiocManager::install(\"scDblFinder\")'",
                }
            }
        
        # Convert to R
        ro.globalenv["counts_matrix"] = adata.X.T.tocsr() if hasattr(adata.X, 'tocsr') else adata.X
        ro.globalenv["n_cells"] = n_cells
        
        ro.r("""
        library(SingleCellExperiment)
        counts <- as.matrix(counts_matrix)
        rownames(counts) <- NULL
        sce <- SingleCellExperiment(list(counts=counts))
        sce <- scDblFinder::scDblFinder(sce, dbr=0.06, verbose=FALSE)
        """ + f"""
        doublet_class <- as.character(sce$scDblFinder.class)
        doublet_score <- as.numeric(sce$scDblFinder.score)
        """)
        
        doublet_class = ro.globalenv["doublet_class"]
        doublet_scores_r = ro.globalenv["doublet_score"]
        
        adata.obs["is_doublet"] = doublet_class == "doublet"
        adata.obs["doublet_score"] = np.array(doublet_scores_r)
        
        n_doublets = int((adata.obs["is_doublet"]).sum())
        print(f"  [scDblFinder] ✅ Detected {n_doublets}/{n_cells} doublets")
        
        filter_doublets = params.get("filter_doublets", False)
        if filter_doublets:
            adata = adata[~adata.obs["is_doublet"]].copy()
        
        return {
            "adata": adata,
            "doublet_report": {
                "method": "scDblFinder",
                "n_cells_input": int(state["adata"].n_obs),
                "n_cells_output": int(adata.n_obs),
                "n_doublets": n_doublets,
                "pct_doublets": round(n_doublets / n_cells, 4),
                "status": "success",
            }
        }
