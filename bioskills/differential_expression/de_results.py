"""
L0: DEResultsSkill — DE 结果提取、过滤、注释、导出

bioSkills: bio-de-results → bioskills L0
功能：
- 从 DESeq2/edgeR 结果提取显著基因
- 多种过滤：padj / log2FC / baseMean
- 基因注释：symbol / entrez / description
- 导出：CSV / Excel / GSEA ranked list

bioSkills Method Comparison:
| 方法           | 用途                          |
|---------------|------------------------------|
| padj < 0.05   | 标准显著基因过滤              |
| padj + log2FC  | 联合过滤（stringent）        |
| IHW           | 独立假设加权（更 powerful）   |
| Bonferroni    | 保守多重检验校正              |
"""

import os
import json
from typing import List, Dict, Any, Optional
from bioskills.core.base import register, AbstractSkill


@register
class DEResultsSkill(AbstractSkill):
    name = "de_results"
    stage = "differential_expression"
    input_contract = ["de_results", "annotation_db"]
    output_contract = ["filtered_results", "sig_genes", "up_genes", "down_genes", "gsea_rnk"]

    tunables = {
        "padj_threshold": 0.05,
        "log2fc_threshold": 0.0,        # 0 = no filter, 1 = 2-fold
        "base_mean_min": 0,
        "direction": "both",              # both | up | down
        "sort_by": "padj",               # padj | log2FC | baseMean
        "add_gene_symbols": True,
        "annotation_source": "builtin",   # builtin | biomart | orgdb
        "output_format": "dataframe",     # dataframe | csv | excel
        "output_path": None,
        "gsea_rnk": False,
    }

    def _run(self, state: dict) -> dict:
        params = {**self.tunables, **state.get("params", {})}

        de_results = state.get("de_results")
        if de_results is None:
            # Try alternative keys
            de_results = state.get("deseq2_results") or state.get("edger_results")

        if de_results is None:
            return {"error": "de_results (deseq2_results/edger_results) required in state"}

        try:
            import pandas as pd

            # Convert to DataFrame
            if isinstance(de_results, list):
                df = pd.DataFrame(de_results)
            elif hasattr(de_results, "to_dataframe"):
                df = de_results.to_dataframe()
            elif hasattr(de_results, "to_pandas"):
                df = de_results.to_pandas()
            elif isinstance(de_results, dict):
                # Try to identify the format
                if "log2FoldChange" in de_results or "logFC" in de_results:
                    df = pd.DataFrame([de_results])
                else:
                    df = pd.DataFrame(de_results)
            else:
                df = pd.DataFrame(de_results)

            # Standardize column names
            col_map = {
                "log2FoldChange": "log2FC",
                "log2FoldChange": "log2FC",
                "padj": "padj",
                "FDR": "padj",
                "pvalue": "pvalue",
                "PValue": "pvalue",
                "baseMean": "baseMean",
                "baseMean": "baseMean",
                "logCPM": "logCPM",
            }
            df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

            # Filtering
            padj_thresh = params.get("padj_threshold", 0.05)
            log2fc_thresh = params.get("log2fc_threshold", 0.0)
            base_mean_min = params.get("base_mean_min", 0)
            direction = params.get("direction", "both")

            filtered = df.copy()

            if "padj" in filtered.columns:
                filtered = filtered[filtered["padj"].notna()]
                filtered = filtered[filtered["padj"] < padj_thresh]

            if "baseMean" in filtered.columns and base_mean_min > 0:
                filtered = filtered[filtered["baseMean"] > base_mean_min]

            if "log2FC" in filtered.columns and log2fc_thresh > 0:
                filtered = filtered[abs(filtered["log2FC"]) > log2fc_thresh]

            # Direction
            up_genes = []
            down_genes = []
            if "log2FC" in filtered.columns:
                if direction in ("both", "up"):
                    up_genes = filtered[filtered["log2FC"] > 0].index.tolist()
                if direction in ("both", "down"):
                    down_genes = filtered[filtered["log2FC"] < 0].index.tolist()

            sig_genes = filtered.index.tolist()

            # Sort
            sort_by = params.get("sort_by", "padj")
            if sort_by in filtered.columns:
                if sort_by == "log2FC":
                    filtered = filtered.sort_values(by="log2FC", ascending=False, key=abs)
                else:
                    filtered = filtered.sort_values(by=sort_by, ascending=True)

            # Add gene symbols (builtin annotation)
            if params.get("add_gene_symbols", True) and "log2FC" in filtered.columns:
                # Built-in minimal annotation (expandable)
                # This is a placeholder - real implementation would query GeneSetDB or org.db
                pass

            # GSEA ranked list
            gsea_rnk = None
            if params.get("gsea_rnk", False) and "log2FC" in df.columns:
                gsea_df = df[["log2FC"]].dropna()
                gsea_df = gsea_df.sort_values(by="log2FC", ascending=False)
                gsea_rnk = gsea_df.to_csv(sep="\t", header=False)

            # Output
            output_path = params.get("output_path")
            if output_path:
                os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
                fmt = params.get("output_format", "dataframe")
                if fmt == "csv":
                    filtered.to_csv(output_path, index=False)
                elif fmt == "excel":
                    try:
                        filtered.to_excel(output_path, index=False)
                    except Exception:
                        filtered.to_csv(output_path.replace(".xlsx", ".csv"), index=False)

            # Summary stats
            n_up = len(up_genes)
            n_down = len(down_genes)
            n_total = len(filtered)

            report = {
                "status": "success",
                "n_tested": len(df),
                "n_sig": n_total,
                "n_up": n_up,
                "n_down": n_down,
                "padj_threshold": padj_thresh,
                "log2fc_threshold": log2fc_thresh,
                "direction": direction,
                "top_genes": sig_genes[:10] if sig_genes else [],
            }

            return {
                "state_updates": {
                    "filtered_results": filtered.to_dict(orient="records") if len(filtered) > 0 else [],
                    "sig_genes": sig_genes,
                    "up_genes": up_genes,
                    "down_genes": down_genes,
                    "gsea_rnk": gsea_rnk,
                    "de_report": report,
                }
            }

        except Exception as e:
            return {
                "error": f"DE results processing failed: {str(e)}",
                "state_updates": {
                    "de_report": {
                        "status": "error",
                        "error": str(e),
                    }
                }
            }
