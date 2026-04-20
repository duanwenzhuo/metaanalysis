"""
L0: MarkerScore Skill — 通用细胞注释（替换硬编码 T cell markers）

功能：
1. 从 state['params']['markers'] 或 state['params']['cell_type'] 获取 marker 基因
2. GeneSetDB 动态查询（不再硬编码）
3. 对任意细胞类型通用
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
import numpy as np

@register
class MarkerScoreSkill(AbstractSkill):
    """
    通用 marker gene 细胞注释。
    
    动态 marker 来源（优先级）：
    1. state['params']['markers']（显式传入）
    2. state['params']['cell_type'] + GeneSetDB 查询
    3. state['params']['gene_sets']（Phase 2 产出）
    
    输出:
      - is_<cell_type>: bool array in adata.obs
      - <cell_type>_score: float array in adata.obs
      - annotation_report: 报告
    """
    
    name = "marker_score"
    description = (
        "Score cells by marker gene expression. "
        "Automatically resolves markers from GeneSetDB if not provided."
    )
    input_contract = ["adata"]
    output_contract = ["adata", "annotation_report"]
    stage = Stage.ANNOTATION
    modality = [Modality.SCRNA, Modality.SPATIAL, Modality.GENERIC]
    
    tunable_parameters = {
        "markers": {"type": "list", "default": []},
        "cell_type": {"type": "str", "default": "Unknown"},
        "process": {"type": "str", "default": ""},
        "quantile": {"type": "float", "default": 0.75},
        "min_markers_required": {"type": "int", "default": 3},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        import scanpy as sc_utils  # alias to avoid confusion
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        # ── Step 1: 动态获取 markers ─────────────
        markers = params.get("markers", [])
        cell_type = params.get("cell_type", "Unknown")
        process = params.get("process", "")
        
        if not markers:
            # 尝试从 GeneSetDB 动态查询
            from bioskills.knowledge.gene_set_db import resolve_gene_set
            hyp_text = params.get("hypothesis_text", state.get("hypothesis_text", ""))
            
            gene_sets = resolve_gene_set(entity=cell_type,
                relation=process,
                hypothesis_text=hyp_text,
            )
            
            # 取第一个基因集作为主基因集
            if gene_sets:
                main_gs_name = next(iter(gene_sets))
                markers = gene_sets[main_gs_name]
                print(f"  [MarkerScore] Resolved {len(markers)} markers from GeneSetDB "
                      f"(cell_type={cell_type}, process={process})")
            else:
                markers = []
        
        # ── Step 2: 过滤可用基因 ─────────────────
        available = [g for g in markers if g in adata.var_names]
        pct_available = len(available) / len(markers) if markers else 0
        min_required = params.get("min_markers_required", 3)
        
        score_name = f"{cell_type.lower().replace(' ', '_')}_score"
        is_positive_key = f"is_{cell_type.lower().replace(' ', '_')}"
        
        if len(available) < min_required:
            print(f"  [MarkerScore] ⚠️  Only {len(available)}/{len(markers)} markers "
                  f"found, below minimum {min_required} — skipping scoring")
            return {
                "adata": adata,
                "annotation_report": {
                    "cell_type": cell_type,
                    "n_markers_requested": len(markers),
                    "n_markers_available": len(available),
                    "pct_available": pct_available,
                    "status": "skipped",
                    "reason": f"Only {len(available)} markers available",
                }
            }
        
        # ── Step 3: 计算 Marker Score ────────────
        sc.tl.score_genes(adata, available, score_name=score_name, copy=False)
        
        quantile = params.get("quantile", 0.75)
        threshold = float(adata.obs[score_name].quantile(quantile))
        
        adata.obs[is_positive_key] = (
            adata.obs[score_name] > threshold
        ).values
        
        n_positive = int(adata.obs[is_positive_key].sum())
        pct_positive = n_positive / adata.n_obs
        
        report = {
            "cell_type": cell_type,
            "process": process,
            "score_name": score_name,
            "n_markers_requested": len(markers),
            "n_markers_available": len(available),
            "markers_available": available,
            "markers_missing": list(set(markers) - set(available)),
            "pct_available": round(pct_available, 3),
            "quantile_threshold": quantile,
            "score_threshold": round(threshold, 4),
            "n_positive_cells": n_positive,
            "pct_positive_cells": round(pct_positive, 4),
            "status": "success",
        }
        
        print(f"  [MarkerScore] ✅ {cell_type}: {n_positive}/{adata.n_obs} "
              f"({pct_positive:.1%}) cells positive "
              f"(markers: {len(available)}/{len(markers)} found, "
              f"{len(set(markers) - set(available))} missing)")
        
        return {"adata": adata, "annotation_report": report}
