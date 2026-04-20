"""
L0: GSEASkill — GSEA 富集分析

bioSkills: bio-pathway-analysis-gsea → bioskills L0
支持: GSEA (预排名) / Preranked / Over-Representation Analysis
基因集: MSigDB / Custom / GO / KEGG

bioSkills Method Comparison:
| 方法        | 输入           | 适用场景                     |
|------------|---------------|------------------------------|
| GSEA       | ranked list   | 无阈值偏倚，发现通路         |
| ORA        | 基因列表      | 快速富集，显著基因           |
| ssGSEA     | 表达矩阵      | 单样本评分                   |
| GSVA       | 表达矩阵      | 连续评分（已有gsva skill）  |
"""

import os
from typing import Optional, List, Dict, Any
from bioskills.core.base import AbstractSkill


class GSEASkill(AbstractSkill):
    name = "gsea"
    stage = "enrichment"
    input_contract = ["gene_list", "gene_sets"]
    output_contract = ["gsea_results", "gsea_plot", "enriched_pathways"]

    tunables = {
        "method": "preranked",           # preranked | ora | ssgsea
        "gene_sets": "MSigDB_Hallmark", # MSigDB_Hallmark | GO_BP | KEGG | custom
        "min_size": 15,
        "max_size": 500,
        "permutation_num": 1000,
        "seed": 42,
        "n_jobs": 1,
        "outdir": None,
        "format": "png",                # png | pdf | svg
        "show_plots": False,
        "threshold": 0.05,              # NES FDR threshold
        "gene_sets_custom": None,       # dict: {set_name: [gene_list]}
        "organism": "human",            # human | mouse
    }

    MSIGDB_MAP = {
        "MSigDB_Hallmark": "h.all.v2023.2.Hs.symbols.gmt",
        "GO_BP": "c5.go.bp.v2023.2.Hs.symbols.gmt",
        "GO_MF": "c5.go.mf.v2023.2.Hs.symbols.gmt",
        "GO_CC": "c5.go.cc.v2023.2.Hs.symbols.gmt",
        "KEGG": "c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt",
        "Reactome": "c2.cp.reactome.v2023.2.Hs.symbols.gmt",
        "BioCarta": "c2.cp.biocarta.v2023.2.Hs.symbols.gmt",
    }

    def _get_gene_sets(self, params: dict):
        """Resolve gene sets - use built-in MSigDB or custom."""
        import gseapy as gp

        custom = params.get("gene_sets_custom")
        if custom:
            return custom

        gs_name = params.get("gene_sets", "MSigDB_Hallmark")
        # gseapy supports built-in names directly
        return gs_name

    def _run(self, state: dict) -> dict:
        params = {**self.tunables, **state.get("params", {})}

        method = params.get("method", "preranked")

        try:
            import gseapy as gp
            import pandas as pd
        except ImportError:
            return {"error": "gseapy required: pip install gseapy"}

        try:
            gene_sets = self._get_gene_sets(params)

            if method == "preranked":
                return self._run_preranked(state, params, gene_sets, gp)
            elif method == "ora":
                return self._run_ora(state, params, gene_sets, gp)
            elif method == "ssgsea":
                return self._run_ssgsea(state, params, gene_sets, gp)
            else:
                return {"error": f"Unknown GSEA method: {method}"}

        except Exception as e:
            return {
                "error": f"GSEA failed: {str(e)}",
                "state_updates": {
                    "gsea_report": {"status": "error", "error": str(e), "method": method}
                }
            }

    def _run_preranked(self, state, params, gene_sets, gp):
        import pandas as pd

        # Get ranked gene list
        rnk = state.get("gsea_rnk") or state.get("ranked_genes")
        if rnk is None:
            # Try to build from DE results
            de_results = state.get("de_results") or state.get("deseq2_results") or state.get("edger_results")
            if de_results is not None:
                if isinstance(de_results, list):
                    df = pd.DataFrame(de_results)
                elif hasattr(de_results, "to_dataframe"):
                    df = de_results.to_dataframe()
                else:
                    df = pd.DataFrame(de_results)

                # Try to find ranking column
                rank_col = None
                for col in ["log2FC", "log2FoldChange", "logFC", "stat"]:
                    if col in df.columns:
                        rank_col = col
                        break

                if rank_col:
                    rnk = df[[rank_col]].dropna().sort_values(by=rank_col, ascending=False)
                else:
                    return {"error": "Cannot build ranked list from DE results - no log2FC/stat column"}
            else:
                return {"error": "gsea_rnk or ranked_genes required for preranked GSEA"}

        outdir = params.get("outdir") or "gsea_output"

        gs = gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets,
            min_size=params.get("min_size", 15),
            max_size=params.get("max_size", 500),
            permutation_num=params.get("permutation_num", 1000),
            seed=params.get("seed", 42),
            no_plot=not params.get("show_plots", False),
            outdir=outdir,
            processes=params.get("n_jobs", 1),
            format=params.get("format", "png"),
        )

        # Extract results
        res_df = gs.res2d
        threshold = params.get("threshold", 0.05)

        # Filter significant pathways
        if "FDR q-val" in res_df.columns:
            sig_pathways = res_df[res_df["FDR q-val"].astype(float) < threshold]
        elif "Adjusted P-value" in res_df.columns:
            sig_pathways = res_df[res_df["Adjusted P-value"].astype(float) < threshold]
        else:
            sig_pathways = res_df

        enriched = []
        for _, row in sig_pathways.iterrows():
            enriched.append({
                "pathway": row.get("Term", row.get("Gene_set", "")),
                "nes": float(row.get("NES", 0)),
                "pval": float(row.get("NOM p-val", 1)),
                "fdr": float(row.get("FDR q-val", 1)),
                "size": int(row.get("Tag %", 0)),
            })

        enriched.sort(key=lambda x: x["fdr"])

        return {
            "state_updates": {
                "gsea_results": res_df.to_dict(orient="records"),
                "enriched_pathways": enriched,
                "gsea_report": {
                    "status": "success",
                    "method": "preranked",
                    "n_pathways_tested": len(res_df),
                    "n_enriched": len(enriched),
                    "gene_sets": str(gene_sets),
                    "threshold": threshold,
                    "outdir": outdir,
                }
            }
        }

    def _run_ora(self, state, params, gene_sets, gp):
        # Get gene list (not ranked)
        gene_list = state.get("sig_genes") or state.get("gene_list")
        if gene_list is None:
            return {"error": "sig_genes or gene_list required for ORA"}

        if isinstance(gene_list, str):
            gene_list = [g.strip() for g in gene_list.split(",")]

        outdir = params.get("outdir") or "enrichr_output"

        # Try Enrichr
        try:
            gs = gp.enrichr(
                gene_list=gene_list,
                gene_sets=gene_sets,
                outdir=outdir,
                no_plot=not params.get("show_plots", False),
            )
            res_df = gs.results
        except Exception:
            # Fallback to local ORA
            return {"error": "Enrichr API failed. Use preranked method with ranked gene list."}

        threshold = params.get("threshold", 0.05)
        if "Adjusted P-value" in res_df.columns:
            sig = res_df[res_df["Adjusted P-value"].astype(float) < threshold]
        else:
            sig = res_df

        enriched = []
        for _, row in sig.iterrows():
            enriched.append({
                "pathway": row.get("Term", ""),
                "pval": float(row.get("P-value", 1)),
                "fdr": float(row.get("Adjusted P-value", 1)),
                "genes": row.get("Genes", ""),
                "overlap": row.get("Overlap", ""),
            })

        return {
            "state_updates": {
                "gsea_results": res_df.to_dict(orient="records"),
                "enriched_pathways": enriched,
                "gsea_report": {
                    "status": "success",
                    "method": "ora",
                    "n_enriched": len(enriched),
                    "gene_sets": str(gene_sets),
                }
            }
        }

    def _run_ssgsea(self, state, params, gene_sets, gp):
        import pandas as pd

        adata = state.get("adata")
        if adata is None:
            return {"error": "adata required for ssGSEA"}

        if hasattr(adata, "to_df"):
            expr_df = adata.to_df().T  # genes x samples
        else:
            return {"error": "adata must support to_df() for ssGSEA"}

        outdir = params.get("outdir") or "ssgsea_output"

        gs = gp.ssgsea(
            data=expr_df,
            gene_sets=gene_sets,
            min_size=params.get("min_size", 15),
            max_size=params.get("max_size", 500),
            permutation_num=params.get("permutation_num", 1000),
            seed=params.get("seed", 42),
            no_plot=not params.get("show_plots", False),
            outdir=outdir,
            processes=params.get("n_jobs", 1),
        )

        # ssGSEA scores per sample per pathway
        scores = gs.res2d

        return {
            "state_updates": {
                "gsea_results": scores.to_dict(orient="records"),
                "ssgsea_scores": scores,
                "gsea_report": {
                    "status": "success",
                    "method": "ssgsea",
                    "n_pathways": len(scores["Term"].unique()) if "Term" in scores.columns else 0,
                    "gene_sets": str(gene_sets),
                }
            }
        }
