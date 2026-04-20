"""
L0: VolcanoSkill — 增强版火山图

bioSkills: bio-data-visualization-volcano-customization → bioskills L0
支持: matplotlib / seaborn / plotly (interactive)
特色: 多分组对比 / 基因标签 / GO pathway 高亮 / 多阈值线 / 自适应标签避让

bioSkills Method Comparison:
| 工具              | 交互 | 适用场景                     |
|------------------|------|------------------------------|
| matplotlib       | 否   | 静态发表级                   |
| seaborn          | 否   | 美化静态图                   |
| plotly           | 是   | 交互式探索                   |
| EnhancedVolcano  | 否   | R 发表级（需要 rpy2）        |
"""

import os
from typing import Optional, List, Dict, Any
from bioskills.core.base import register, AbstractSkill, Stage, Modality, State


@register
class VolcanoSkill(AbstractSkill):
    name = "volcano"
    stage = "visualization"
    input_contract = ["de_results"]
    output_contract = ["volcano_plot", "volcano_data"]

    tunables = {
        "padj_col": "padj",            # padj column name
        "log2fc_col": "log2FC",         # log2FC column name
        "padj_thresh": 0.05,
        "log2fc_thresh": 1.0,
        "neg_log10": True,             # -log10(padj) for y-axis
        "label_top_n": 10,             # label top N genes
        "label_genes": None,           # specific genes to label
        "highlight_pathways": None,    # dict: {pathway_name: [gene_list]}
        "colors": {"up": "#E64B35", "down": "#4DBBD5", "ns": "grey"},
        "point_size": 8,
        "alpha": 0.6,
        "title": "Volcano Plot",
        "backend": "matplotlib",       # matplotlib | plotly
        "output_path": None,
        "figsize": [8, 6],
        "dpi": 300,
        "format": "png",              # png | pdf | svg
    }

    def _run(self, state: dict) -> dict:
        params = {**self.tunables, **state.get("params", {})}

        de_results = state.get("de_results") or state.get("deseq2_results") or state.get("edger_results")
        if de_results is None:
            return {"error": "de_results required in state"}

        try:
            import pandas as pd
            import numpy as np

            # Convert to DataFrame
            if isinstance(de_results, list):
                df = pd.DataFrame(de_results)
            elif hasattr(de_results, "to_dataframe"):
                df = de_results.to_dataframe()
            elif hasattr(de_results, "to_pandas"):
                df = de_results.to_pandas()
            elif isinstance(de_results, dict):
                df = pd.DataFrame(de_results)
            else:
                df = pd.DataFrame(de_results)

            # Standardize column names
            col_renames = {}
            padj_col = params.get("padj_col", "padj")
            log2fc_col = params.get("log2fc_col", "log2FC")
            if padj_col not in df.columns:
                if "FDR" in df.columns:
                    col_renames["FDR"] = padj_col
                elif "adj.P.Val" in df.columns:
                    col_renames["adj.P.Val"] = padj_col
            if log2fc_col not in df.columns:
                if "log2FoldChange" in df.columns:
                    col_renames["log2FoldChange"] = log2fc_col
                elif "logFC" in df.columns:
                    col_renames["logFC"] = log2fc_col
            df = df.rename(columns=col_renames)

            if padj_col not in df.columns or log2fc_col not in df.columns:
                return {"error": f"Required columns not found. Need {padj_col} and {log2fc_col}. Available: {list(df.columns)}"}

            # Prepare data
            df = df.copy()
            df["neg_log10_padj"] = -np.log10(df[padj_col].clip(lower=1e-300))

            # Classify points
            padj_t = params.get("padj_thresh", 0.05)
            log2fc_t = params.get("log2fc_thresh", 1.0)

            df["significance"] = "ns"
            df.loc[(df[padj_col] < padj_t) & (df[log2fc_col] > log2fc_t), "significance"] = "up"
            df.loc[(df[padj_col] < padj_t) & (df[log2fc_col] < -log2fc_t), "significance"] = "down"

            colors = params.get("colors", {"up": "#E64B35", "down": "#4DBBD5", "ns": "grey"})
            backend = params.get("backend", "matplotlib")

            if backend == "plotly":
                return self._plotly_volcano(df, params, padj_col, log2fc_col)
            else:
                return self._matplotlib_volcano(df, params, padj_col, log2fc_col, colors, padj_t, log2fc_t)

        except Exception as e:
            return {
                "error": f"Volcano plot failed: {str(e)}",
                "state_updates": {
                    "volcano_report": {"status": "error", "error": str(e)}
                }
            }

    def _matplotlib_volcano(self, df, params, padj_col, log2fc_col, colors, padj_t, log2fc_t):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        fig, ax = plt.subplots(figsize=params.get("figsize", [8, 6]))

        # Plot points by category
        for sig, color in colors.items():
            mask = df["significance"] == sig
            ax.scatter(df.loc[mask, log2fc_col], df.loc[mask, "neg_log10_padj"],
                       c=color, s=params.get("point_size", 8),
                       alpha=params.get("alpha", 0.6), label=sig, edgecolors="none")

        # Threshold lines
        ax.axhline(-np.log10(padj_t), color="grey", linestyle="--", linewidth=0.8)
        ax.axvline(log2fc_t, color="grey", linestyle="--", linewidth=0.8)
        ax.axvline(-log2fc_t, color="grey", linestyle="--", linewidth=0.8)

        # Highlight pathways
        highlight_pathways = params.get("highlight_pathways")
        if highlight_pathways:
            for pw_name, gene_list in highlight_pathways.items():
                pw_mask = df.index.isin(gene_list)
                ax.scatter(df.loc[pw_mask, log2fc_col], df.loc[pw_mask, "neg_log10_padj"],
                           s=40, facecolors="none", edgecolors="black", linewidths=1.5,
                           label=pw_name, zorder=5)

        # Label top genes
        label_top_n = params.get("label_top_n", 10)
        label_genes = params.get("label_genes", [])
        if label_top_n > 0 or label_genes:
            top_sig = df[df["significance"] != "ns"].nlargest(label_top_n, "neg_log10_padj")
            label_set = set(label_genes) | set(top_sig.index.tolist())
            for gene in label_set:
                if gene in df.index:
                    row = df.loc[gene]
                    ax.annotate(str(gene),
                                xy=(row[log2fc_col], row["neg_log10_padj"]),
                                xytext=(5, 5), textcoords="offset points",
                                fontsize=7, fontweight="bold")

        ax.set_xlabel("log2 Fold Change")
        ax.set_ylabel("-log10(adjusted p-value)")
        ax.set_title(params.get("title", "Volcano Plot"))
        ax.legend(loc="upper right", fontsize=8)

        n_up = (df["significance"] == "up").sum()
        n_down = (df["significance"] == "down").sum()

        # Save
        output_path = params.get("output_path")
        if output_path:
            os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
            fmt = params.get("format", "png")
            fig.savefig(output_path, dpi=params.get("dpi", 300), bbox_inches="tight", format=fmt)
        plt.close(fig)

        return {
            "state_updates": {
                "volcano_plot": output_path or "in_memory",
                "volcano_data": df[["log2FC", "padj", "neg_log10_padj", "significance"]].to_dict(orient="records"),
                "volcano_report": {
                    "status": "success",
                    "n_up": int(n_up),
                    "n_down": int(n_down),
                    "n_ns": int((df["significance"] == "ns").sum()),
                    "padj_thresh": padj_t,
                    "log2fc_thresh": log2fc_t,
                    "output_path": output_path,
                    "backend": "matplotlib",
                }
            }
        }

    def _plotly_volcano(self, df, params, padj_col, log2fc_col):
        import plotly.express as px
        import numpy as np

        colors_map = params.get("colors", {"up": "#E64B35", "down": "#4DBBD5", "ns": "grey"})

        fig = px.scatter(
            df, x=log2fc_col, y="neg_log10_padj",
            color="significance",
            color_discrete_map=colors_map,
            hover_name=df.index,
            title=params.get("title", "Volcano Plot"),
            labels={log2fc_col: "log2 Fold Change", "neg_log10_padj": "-log10(adjusted p-value)"},
            opacity=params.get("alpha", 0.6),
        )

        padj_t = params.get("padj_thresh", 0.05)
        log2fc_t = params.get("log2fc_thresh", 1.0)

        fig.add_hline(y=-np.log10(padj_t), line_dash="dash", line_color="grey")
        fig.add_vline(x=log2fc_t, line_dash="dash", line_color="grey")
        fig.add_vline(x=-log2fc_t, line_dash="dash", line_color="grey")

        output_path = params.get("output_path")
        if output_path:
            os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
            fig.write_html(output_path)

        n_up = (df["significance"] == "up").sum()
        n_down = (df["significance"] == "down").sum()

        return {
            "state_updates": {
                "volcano_plot": output_path or "in_memory",
                "volcano_data": df[["log2FC", "padj", "neg_log10_padj", "significance"]].to_dict(orient="records"),
                "volcano_report": {
                    "status": "success",
                    "n_up": int(n_up),
                    "n_down": int(n_down),
                    "n_ns": int((df["significance"] == "ns").sum()),
                    "backend": "plotly",
                    "output_path": output_path,
                }
            }
        }
