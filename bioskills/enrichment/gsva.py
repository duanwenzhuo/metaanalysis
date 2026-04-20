"""
L0: GSVA Skill — 通用基因集富集评分（替换硬编码 exhaustion genes）

功能：
1. 接收任意 gene_sets dict（不再硬编码 18 个 exhaustion 基因）
2. 支持 gsva / ssgsea / aucell 三种方法
3. 自动过滤可用基因
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
import numpy as np

@register
class GSVASkill(AbstractSkill):
    """
    通用 GSVA / ssGSEA 基因集富集评分。
    
    Gene sets 来源（优先级）：
    1. state['params']['gene_sets']（Phase 2 / GeneSetDB 产出）
    2. state['literature_miner_gene_sets']（Phase 2 节点产出）
    3. 兜底为空（跳过）
    
    输出:
      - gsva_scores: dict of gene_set → scores dict
      - adata: obs 中增加 gsva_<gene_set_name> 列
    """
    
    name = "gsva"
    description = (
        "GSVA / ssGSEA / AUCell gene set enrichment scoring. "
        "Accepts any gene_sets dict. No hardcoded gene lists."
    )
    input_contract = ["adata"]
    output_contract = ["adata", "gsva_report"]
    stage = Stage.ENRICHMENT
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.SPATIAL, Modality.GENERIC]
    
    tunable_parameters = {
        "gene_sets": {"type": "dict", "default": {}},
        "method": {"type": "str", "default": "ssgsea"},
        "min_size": {"type": "int", "default": 5},
        "max_size": {"type": "int", "default": 500},
        "sample_norm": {"type": "str", "default": "rank"},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        import gseapy as gp
        
        adata = state["adata"]
        params = state.get("params", {})
        
        # 动态获取 gene_sets（优先级：params > state）
        gene_sets = params.get("gene_sets", {})
        
        # 尝试从 state 中的其他字段获取
        if not gene_sets:
            for key in ["literature_miner_gene_sets", 
                        "hypothesis_parser_gene_sets",
                        "gene_sets"]:
                if key in state and state[key]:
                    gene_sets = state[key]
                    break
        
        if not gene_sets:
            print(f"  [{self.name}] ⚠️  No gene_sets provided — skipping")
            return {
                "adata": adata,
                "gsva_report": {"status": "skipped", "reason": "no gene_sets"}
            }
        
        # 过滤可用基因
        expr = adata.to_df().T  # genes × cells
        available_gs = {}
        
        for gs_name, genes in gene_sets.items():
            avail = [g for g in genes if g in expr.index]
            if len(avail) >= params.get("min_size", 5):
                available_gs[gs_name] = avail
            else:
                print(f"  [{self.name}] ⚠️  {gs_name}: only {len(avail)}/{len(genes)} genes found")
        
        if not available_gs:
            return {
                "adata": adata,
                "gsva_report": {"status": "skipped", "reason": "no valid gene sets"}
            }
        
        method = params.get("method", "ssgsea")
        
        try:
            if method == "ssgsea":
                res = gp.ssgsea(
                    data=expr,
                    gene_sets=available_gs,
                    sample_norm_method=params.get("sample_norm", "rank"),
                    min_size=params.get("min_size", 5),
                    max_size=params.get("max_size", 500),
                    no_plot=True,
                    outdir=None,
                )
            elif method == "gsva":
                res = gp.gsva(
                    expr=expr,
                    gene_sets=available_gs,
                    method="gsva",
                    min_size=params.get("min_size", 5),
                    max_size=params.get("max_size", 500),
                    no_plot=True,
                )
            elif method == "aucell":
                res = gp.aucell(
                    expr=expr,
                    gene_sets=available_gs,
                    min_size=params.get("min_size", 5),
                    max_size=params.get("max_size", 500),
                )
            else:
                raise ValueError(f"Unknown method: {method}")
            
            scores = res.resultsOnSamples.T  # cells × gene_sets
            
            # 写回 adata.obs
            score_cols = []
            for gs_name, row in scores.iterrows():
                col_name = f"gsva_{gs_name}"
                adata.obs[col_name] = adata.obs_names.map(row.to_dict()).values
                score_cols.append(col_name)
            
            report = {
                "method": method,
                "gene_sets_requested": list(gene_sets.keys()),
                "gene_sets_used": list(available_gs.keys()),
                "gene_sets_skipped": [
                    k for k in gene_sets if k not in available_gs
                ],
                "score_columns": score_cols,
                "n_cells": adata.n_obs,
                "status": "success",
            }
            
            print(f"  [GSVA] ✅ {method}: {len(available_gs)}/{len(gene_sets)} "
                  f"gene sets scored across {adata.n_obs} cells")
            
            return {"adata": adata, "gsva_report": report}
            
        except Exception as e:
            print(f"  [GSVA] ❌ Failed: {e}")
            return {
                "adata": adata,
                "gsva_report": {"status": "failed", "error": str(e)}
            }
