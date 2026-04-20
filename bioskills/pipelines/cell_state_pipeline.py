"""
CellStatePipeline — 通用 Cell State 分析 Pipeline

功能：
- 完全通用：任意细胞类型 + 任意生物过程
- 基因集从 GeneSetDB 动态获取（不再硬编码）
- 所有参数从 state['params'] 读取
- L0 Skill 串接：qc → normalize → hvg → scale → pca → neighbors → leiden → marker_score → gsva → cohens_d
"""

from __future__ import annotations
from typing import Dict, Any, List, Optional
from pathlib import Path

from bioskills.core.base import (
    register, AbstractSkill, Stage, Modality,
    SkillRegistry, State,
)


@register
class CellStatePipeline(AbstractSkill):
    """
    通用 Cell State 分析 Pipeline。
    
    取代 phase4_single_dataset.py 中的硬编码 PARAMS。
    
    输入 State:
      adata: AnnData 对象（必需）
      params.gene_sets: dict[str, list[str]]（可选，动态获取）
      params.cell_type: str（可选，用于 GeneSetDB 查询）
      params.process: str（可选）
      params.hypothesis_text: str（可选）
      params.groups: dict（可选，barcode → tumor/normal）
    
    输出 State:
      effect_size_result: {effect_size, p_value, direction, confidence, ...}
      adata: 处理后的 AnnData
    """
    
    name = "cell_state_pipeline"
    description = (
        "General-purpose Cell State analysis pipeline: "
        "QC→Normalize→HVG→PCA→Cluster→Annotate→GSVA→CohensD. "
        "No hardcoded genes. Resolves markers from GeneSetDB."
    )
    input_contract = ["adata"]
    output_contract = ["adata", "effect_size_result"]
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.GENERIC]
    
    def _run(self, state: State) -> Dict[str, Any]:
        import scanpy as sc
        
        registry = SkillRegistry()
        
        # ── Step 1: 动态获取基因集（替换硬编码）───
        params = state.get("params", {})
        hyp_text = params.get("hypothesis_text", 
                             state.get("hypothesis_text", ""))
        cell_type = params.get("cell_type", 
                               state.get("hypothesis_parser_entity1", "Unknown"))
        process = params.get("process",
                            state.get("hypothesis_parser_relation", ""))
        
        manual_genes = params.get("gene_sets", {})
        
        if manual_genes:
            gene_sets = manual_genes
            print(f"[CellStatePipeline] Using manually provided gene_sets "
                  f"({len(gene_sets)} sets)")
        else:
            from bioskills.knowledge.gene_set_db import resolve_gene_set
            gene_sets = resolve_gene_set(entity=cell_type,
                relation=process,
                hypothesis_text=hyp_text,
            )
            print(f"[CellStatePipeline] Resolved gene_sets from GeneSetDB: "
                  f"{list(gene_sets.keys())}")
        
        # ── Step 2: QC ────────────────────────────
        print(f"[CellStatePipeline] Step 1/10: QC")
        qc = registry.get("qc")
        qc_out = qc.execute(state)
        
        if qc_out["status"] != "success":
            return {
                "effect_size_result": {
                    "status": "failed",
                    "error": f"QC failed: {qc_out.get('error')}"
                }
            }
        
        adata = qc_out["state_updates"]["adata"]
        
        # ── Step 3: Normalize ────────────────────
        print(f"[CellStatePipeline] Step 2/10: Normalize")
        norm_out = registry.get("normalize").execute({"adata": adata, "params": params})
        adata = norm_out["state_updates"].get("adata", adata)
        
        # ── Step 4: HVG ──────────────────────────
        print(f"[CellStatePipeline] Step 3/10: HVG")
        hvg_out = registry.get("hvg").execute({"adata": adata, "params": params})
        adata = hvg_out["state_updates"].get("adata", adata)
        
        # ── Step 5: Scale ─────────────────────────
        print(f"[CellStatePipeline] Step 4/10: Scale")
        scale_out = registry.get("scale").execute({"adata": adata, "params": params})
        adata = scale_out["state_updates"].get("adata", adata)
        
        # ── Step 6: PCA ──────────────────────────
        print(f"[CellStatePipeline] Step 5/10: PCA")
        pca_out = registry.get("pca").execute({"adata": adata, "params": params})
        adata = pca_out["state_updates"].get("adata", adata)
        
        # ── Step 7: Neighbors + Leiden ───────────
        print(f"[CellStatePipeline] Step 6/10: Neighbors + Leiden")
        neighbors_out = registry.get("neighbors").execute({
            "adata": adata,
            "params": params,
        })
        adata = neighbors_out["state_updates"].get("adata", adata)
        
        leiden_out = registry.get("leiden").execute({
            "adata": adata,
            "params": params,
        })
        adata = leiden_out["state_updates"].get("adata", adata)
        n_clusters = leiden_out["state_updates"].get("n_clusters", 0)
        print(f"  Leiden: {n_clusters} clusters")
        
        # ── Step 8: Marker Score（通用细胞注释）──
        print(f"[CellStatePipeline] Step 7/10: Marker Score ({cell_type})")
        marker_out = registry.get("marker_score").execute({
            "adata": adata,
            "params": {
                **params,
                "cell_type": cell_type,
                "process": process,
                "markers": params.get("markers", []),
                "gene_sets": gene_sets,
            },
        })
        
        adata = marker_out["state_updates"].get("adata", adata)
        annot_report = marker_out["state_updates"].get("annotation_report", {})
        
        # 检查细胞数是否足够
        is_positive_key = f"is_{cell_type.lower().replace(' ', '_')}"
        if is_positive_key in adata.obs.columns:
            n_target = int(adata.obs[is_positive_key].sum())
        else:
            n_target = adata.n_obs
        
        if n_target < 100:
            return {
                "effect_size_result": {
                    "status": "skipped",
                    "reason": f"Too few target cells: {n_target} (need ≥100)",
                    "annotation_report": annot_report,
                },
                "adata": adata,
            }
        
        # ── Step 9: GSVA ─────────────────────────
        print(f"[CellStatePipeline] Step 8/10: GSVA")
        gsva_out = registry.get("gsva").execute({
            "adata": adata,
            "params": {
                **params,
                "gene_sets": gene_sets,
            },
        })
        adata = gsva_out["state_updates"].get("adata", adata)
        
        # ── Step 10: Cohen's d ───────────────────
        print(f"[CellStatePipeline] Step 9/10: Cohen's d")
        cohens_out = registry.get("cohens_d").execute({
            "adata": adata,
            "params": params,
        })
        
        effect_result = cohens_out["state_updates"].get("effect_size_result", {})
        
        # ── Step 11: 组装结果 ────────────────────
        final_result = {
            "status": effect_result.get("status", "unknown"),
            **effect_result,
            "cell_type": cell_type,
            "process": process,
            "gene_sets_used": list(gene_sets.keys()),
            "annotation_report": annot_report,
            "n_target_cells": n_target,
            "n_clusters": n_clusters,
            "n_cells_total": adata.n_obs,
        }
        
        print(f"[CellStatePipeline] ✅ Done: "
              f"d={effect_result.get('effect_size', 'N/A')}, "
              f"p={effect_result.get('p_value', 'N/A')}, "
              f"direction={effect_result.get('direction', 'N/A')}")
        
        return {
            "adata": adata,
            "effect_size_result": final_result,
        }
