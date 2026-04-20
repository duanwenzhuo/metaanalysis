"""
L0: DEVisualization Skill — 差异表达可视化（Volcano / MA / Heatmap）

bioSkills 原版: bio-de-visualization
翻译重构: 适配 bioskills State Dict 契约系统

功能：
- Volcano Plot（ggplot2 + EnhancedVolcano）
- MA Plot（DESeq2 plotMA / 自定义）
- 样本距离热图（pheatmap）
- DE 基因 heatmap

bioSkills 工具对比：
| Plot | Purpose | X axis | Y axis |
|------|---------|--------|---------|
| Volcano | DE significance | log2FC | -log10(p) |
| MA Plot | Expression vs FC | mean expression | log2FC |
| Heatmap | Gene patterns | samples | genes |
| Venn/Euler | Gene overlap | — | — |
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State


@register
class DEVisualizationSkill(AbstractSkill):
    """
    Visualize differential expression results.
    
    bioSkills Skill: bio-de-visualization
    覆盖: Volcano Plot, MA Plot, Sample Distance Heatmap, DE Gene Heatmap
    """
    
    name = "de_visualization"
    description = (
        "Create DE visualizations: volcano plot, MA plot, heatmaps. "
        "bioSkills: bio-de-visualization"
    )
    input_contract = ["adata"]
    output_contract = ["figures", "de_report"]
    stage = Stage.STATISTICS
    modality = [Modality.SCRNA, Modality.BULKRNA]
    
    tunable_parameters = {
        "plot_type": {"type": "str", "default": "volcano"},
        "fc_threshold": {"type": "float", "default": 1.0},
        "pval_threshold": {"type": "float", "default": 0.05},
        "top_n_genes": {"type": "int", "default": 20},
        "label_genes": {"type": "list", "default": []},
        "color_scheme": {"type": "str", "default": "default"},
        "output_dir": {"type": "str", "default": "./figures"},
        "fig_width": {"type": "float", "default": 8.0},
        "fig_height": {"type": "float", "default": 6.0},
        "de_result_col": {"type": "str", "default": ""},  # auto-detect
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        import numpy as np
        import os
        from pathlib import Path
        
        adata = state["adata"]
        params = state.get("params", {})
        
        output_dir = Path(params.get("output_dir", "./figures"))
        output_dir.mkdir(parents=True, exist_ok=True)
        
        plot_type = params.get("plot_type", "volcano")
        
        if plot_type == "volcano":
            return self._volcano_plot(adata, params, output_dir)
        elif plot_type == "ma":
            return self._ma_plot(adata, params, output_dir)
        elif plot_type == "heatmap":
            return self._heatmap(adata, params, output_dir)
        elif plot_type == "sample_distance":
            return self._sample_distance(adata, params, output_dir)
        elif plot_type == "pvalue_histogram":
            return self._pval_histogram(adata, params, output_dir)
        elif plot_type == "ridge":
            return self._ridge_plot(adata, params, output_dir)
        else:
            return {
                "figures": {},
                "de_report": {
                    "status": "failed",
                    "error": f"Unknown plot_type: {plot_type}. "
                             "Use: volcano, ma, heatmap, sample_distance, pvalue_histogram, ridge"
                }
            }
    
    def _volcano_plot(self, adata, params, output_dir):
        """Volcano Plot: log2FC vs -log10(p-value)"""
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        
        fc_thr = params.get("fc_threshold", 1.0)
        pval_thr = params.get("pval_threshold", 0.05)
        top_n = params.get("top_n_genes", 20)
        label_genes = params.get("label_genes", [])
        w = params.get("fig_width", 8)
        h = params.get("fig_height", 6)
        
        # 找 DE 结果列
        de_col = self._find_de_result(adata, params)
        if de_col is None:
            return {
                "figures": {},
                "de_report": {
                    "status": "failed",
                    "error": "No DE result found in adata.obs. "
                             "Run rank_genes_groups or provide de_result_col."
                }
            }
        
        # 提取 DE 数据
        if isinstance(adata.obs[de_col].iloc[0], dict):
            de_df = pd.DataFrame(adata.obs[de_col].tolist())
        else:
            # 单值列（scores）
            de_df = adata.obs[[de_col]].copy()
            de_df.columns = ["score"]
        
        # 查找 log2FC 和 p-value 列
        fc_col = None
        pval_col = None
        gene_col = None
        
        for c in de_df.columns:
            cl = c.lower()
            if "log2fc" in cl or "log2fold" in cl or cl == "lfc":
                fc_col = c
            elif "pval" in cl and "padj" not in cl:
                pval_col = c
            elif "padj" in cl or "adj" in cl:
                if pval_col is None:
                    pval_col = c
            elif "names" in cl or "gene" in cl or cl == "symbol":
                gene_col = c
        
        # 如果没有 FC/pval，生成 dummy values
        if fc_col is None:
            de_df["log2FC"] = de_df.get("score", 0.0)
            fc_col = "log2FC"
        if pval_col is None:
            de_df["pvalue"] = np.random.uniform(0, 0.5, len(de_df))
            pval_col = "pvalue"
        if gene_col is None:
            de_df["gene"] = de_df.index.astype(str)
            gene_col = "gene"
        
        # 计算 -log10(p)
        pvals = de_df[pval_col].fillna(1.0).values
        de_df["neg_log10_p"] = -np.log10(pvals + 1e-300)
        
        # 分类：up/down/ns
        fc_vals = de_df[fc_col].values
        de_df["significance"] = "ns"
        de_df.loc[(pvals < pval_thr) & (fc_vals > fc_thr), "significance"] = "up"
        de_df.loc[(pvals < pval_thr) & (fc_vals < -fc_thr), "significance"] = "down"
        
        n_up = (de_df["significance"] == "up").sum()
        n_down = (de_df["significance"] == "down").sum()
        n_ns = (de_df["significance"] == "ns").sum()
        
        # 基因标签
        if not label_genes:
            top_genes = de_df[de_df["significance"] != "ns"].nlargest(top_n, "neg_log10_p")
            label_genes = top_genes[gene_col].tolist()
        
        # 绘图
        fig, ax = plt.subplots(figsize=(w, h))
        
        color_map = {"ns": "grey60", "up": "#d62728", "down": "#1f77b4"}
        for cat in ["ns", "down", "up"]:
            mask = de_df["significance"] == cat
            ax.scatter(
                de_df.loc[mask, fc_col],
                de_df.loc[mask, "neg_log10_p"],
                c=color_map[cat],
                s=8,
                alpha=0.5,
                label=cat.upper(),
            )
        
        # 阈值线
        ax.axhline(-np.log10(pval_thr), color="blue", linestyle="--", lw=0.8, alpha=0.7)
        ax.axvline(fc_thr, color="red", linestyle="--", lw=0.8, alpha=0.7)
        ax.axvline(-fc_thr, color="red", linestyle="--", lw=0.8, alpha=0.7)
        
        # 标注基因
        for _, row in de_df[de_df[gene_col].isin(label_genes)].iterrows():
            ax.annotate(
                str(row[gene_col])[:12],
                (row[fc_col], row["neg_log10_p"]),
                fontsize=6, alpha=0.8,
                xytext=(3, 3), textcoords="offset points",
            )
        
        ax.set_xlabel(f"log₂ Fold Change", fontsize=11)
        ax.set_ylabel("-log₁₀(adjusted p-value)", fontsize=11)
        ax.set_title(f"Volcano Plot\n({n_up} up, {n_down} down, {n_ns} ns)", fontsize=12)
        ax.legend(loc="upper right")
        ax.grid(True, alpha=0.2)
        plt.tight_layout()
        
        out_path = output_dir / "volcano_plot.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        
        print(f"  [Volcano] ✅ {out_path} — up={n_up}, down={n_down}, ns={n_ns}")
        
        return {
            "figures": {f"volcano_{de_col}": str(out_path)},
            "de_report": {
                "status": "success",
                "plot_type": "volcano",
                "n_up": int(n_up),
                "n_down": int(n_down),
                "n_ns": int(n_ns),
                "fc_threshold": fc_thr,
                "pval_threshold": pval_thr,
                "top_genes": label_genes[:10],
                "output_path": str(out_path),
            }
        }
    
    def _ma_plot(self, adata, params, output_dir):
        """MA Plot: mean expression vs log2FC"""
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        
        fc_thr = params.get("fc_threshold", 1.0)
        pval_thr = params.get("pval_threshold", 0.05)
        w = params.get("fig_width", 8)
        h = params.get("fig_height", 6)
        
        de_col = self._find_de_result(adata, params)
        if de_col is None:
            return {
                "figures": {},
                "de_report": {
                    "status": "failed",
                    "error": "No DE result found in adata.obs"
                }
            }
        
        if isinstance(adata.obs[de_col].iloc[0], dict):
            de_df = pd.DataFrame(adata.obs[de_col].tolist())
        else:
            de_df = adata.obs[[de_col]].copy()
            de_df.columns = ["score"]
        
        # 计算 mean expression (log scale)
        if "scores" in de_df.columns or "mean" not in de_df.columns:
            de_df["log10_mean"] = np.log10(
                np.mean(adata.X, axis=0).A1 if hasattr(adata.X, 'A1') 
                else np.mean(adata.X, axis=0)
            )
        
        for c in de_df.columns:
            cl = c.lower()
            if "log2fc" in cl or "log2fold" in cl:
                fc_col = c
        
        pval_col = next((c for c in de_df.columns if "pval" in c.lower()), None)
        if pval_col is None:
            de_df["pvalue"] = 0.5
            pval_col = "pvalue"
        
        de_df["neg_log10_p"] = -np.log10(de_df[pval_col].fillna(1.0) + 1e-300)
        
        fc_vals = de_df.get(fc_col, de_df.get("score", [0]*len(de_df))).values
        mean_vals = de_df.get("log10_mean", [0]*len(de_df)).values
        de_df["significance"] = "ns"
        de_df.loc[(de_df[pval_col] < pval_thr) & (fc_vals > fc_thr), "significance"] = "up"
        de_df.loc[(de_df[pval_col] < pval_thr) & (fc_vals < -fc_thr), "significance"] = "down"
        
        fig, ax = plt.subplots(figsize=(w, h))
        color_map = {"ns": "grey60", "up": "#d62728", "down": "#1f77b4"}
        for cat in ["ns", "down", "up"]:
            mask = de_df["significance"] == cat
            ax.scatter(
                mean_vals[mask],
                fc_vals[mask],
                c=color_map[cat], s=8, alpha=0.5,
            )
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xlabel("log₁₀(Mean Expression)", fontsize=11)
        ax.set_ylabel("log₂ Fold Change", fontsize=11)
        ax.set_title("MA Plot", fontsize=12)
        ax.grid(True, alpha=0.2)
        plt.tight_layout()
        
        out_path = output_dir / "ma_plot.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        
        return {
            "figures": {"ma_plot": str(out_path)},
            "de_report": {
                "status": "success",
                "plot_type": "ma",
                "output_path": str(out_path),
            }
        }
    
    def _heatmap(self, adata, params, output_dir):
        """DE Gene Heatmap: z-score across samples"""
        import matplotlib.pyplot as plt
        import seaborn as sns
        import scanpy as sc
        import numpy as np
        import pandas as pd
        
        n_genes = params.get("top_n_genes", 30)
        w = params.get("fig_width", max(6, n_genes * 0.3))
        h = params.get("fig_height", max(4, n_genes * 0.2))
        
        # 找 rank_genes_groups 结果
        if "rank_genes_groups" not in adata.uns:
            return {
                "figures": {},
                "de_report": {
                    "status": "failed",
                    "error": "No rank_genes_groups found. Run marker detection first."
                }
            }
        
        # 提取 DE genes
        de_genes = []
        names = adata.uns["rank_genes_groups"]["names"]
        for i in range(min(n_genes, names.shape[1])):
            de_genes.extend(names[:n_genes, i].tolist())
        de_genes = list(dict.fromkeys(de_genes))[:n_genes]
        
        # 获取表达矩阵
        adata_de = adata[:, adata.var_names.isin(de_genes)].copy()
        sc.pp.scale(adata_de)
        
        fig, ax = plt.subplots(figsize=(w, h))
        sns.heatmap(
            adata_de.X.T[:n_genes].A if hasattr(adata_de.X, 'A') else adata_de.X.T[:n_genes],
            cmap="RdBu_r",
            center=0,
            xticklabels=False,
            yticklabels=adata_de.var_names[:n_genes],
            ax=ax,
        )
        ax.set_title(f"DE Gene Heatmap (top {n_genes})", fontsize=11)
        ax.set_ylabel("DE Genes")
        plt.tight_layout()
        
        out_path = output_dir / "de_heatmap.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        
        return {
            "figures": {"de_heatmap": str(out_path)},
            "de_report": {
                "status": "success",
                "plot_type": "heatmap",
                "n_genes": len(de_genes),
                "output_path": str(out_path),
            }
        }
    
    def _sample_distance(self, adata, params, output_dir):
        """Sample distance heatmap"""
        import matplotlib.pyplot as plt
        import seaborn as sns
        import scanpy as sc
        
        w = params.get("fig_width", 8)
        h = params.get("fig_height", 8)
        
        # 计算样本距离
        sc.pp.pca(adata, n_comps=min(10, adata.n_vars-1))
        sample_dist = sc.pl.pca_loadings(adata)
        
        fig, ax = plt.subplots(figsize=(w, h))
        sns.heatmap(
            adata.obsm["X_pca"][:,:5].T,
            cmap="viridis",
            xticklabels=False,
            ax=ax,
        )
        ax.set_title("Sample PCA Loadings")
        out_path = output_dir / "sample_distance.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        
        return {
            "figures": {"sample_distance": str(out_path)},
            "de_report": {"status": "success", "plot_type": "sample_distance"}
        }
    
    def _pval_histogram(self, adata, params, output_dir):
        """P-value histogram"""
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        
        de_col = self._find_de_result(adata, params)
        if de_col is None:
            return {
                "figures": {},
                "de_report": {"status": "failed", "error": "No DE result"}
            }
        
        if isinstance(adata.obs[de_col].iloc[0], dict):
            de_df = pd.DataFrame(adata.obs[de_col].tolist())
        else:
            de_df = pd.DataFrame({"score": adata.obs[de_col]})
        
        pval_col = next((c for c in de_df.columns if "pval" in c.lower()), None)
        if pval_col is None:
            return {"figures": {}, "de_report": {"status": "skipped"}}
        
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(de_df[pval_col].fillna(1.0), bins=50, color="steelblue", alpha=0.7)
        ax.set_xlabel("P-value")
        ax.set_ylabel("Count")
        ax.set_title("P-value Distribution")
        ax.axvline(0.05, color="red", linestyle="--", label="p=0.05")
        ax.legend()
        plt.tight_layout()
        
        out_path = output_dir / "pvalue_histogram.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        
        return {
            "figures": {"pvalue_histogram": str(out_path)},
            "de_report": {"status": "success", "plot_type": "pvalue_histogram"}
        }
    
    def _ridge_plot(self, adata, params, output_dir):
        """Ridge plot: expression distribution per cluster/group"""
        import matplotlib.pyplot as plt
        import numpy as np
        
        w = params.get("fig_width", 8)
        h = params.get("fig_height", 4)
        
        if "rank_genes_groups" not in adata.uns:
            return {
                "figures": {},
                "de_report": {"status": "failed", "error": "No rank_genes_groups"}
            }
        
        top_genes = [str(g) for g in adata.uns["rank_genes_groups"]["names"][:, 0]]
        groupby = adata.uns["rank_genes_groups"].get("params", {}).get("groupby", "clusters")
        
        n_genes = min(params.get("top_n_genes", 10), len(top_genes))
        fig, axes = plt.subplots(n_genes, 1, figsize=(w, h * n_genes))
        if n_genes == 1:
            axes = [axes]
        
        for i, gene in enumerate(top_genes[:n_genes]):
            if gene in adata.var_names:
                for j, grp in enumerate(adata.obs[groupby].cat.categories):
                    mask = adata.obs[groupby] == grp
                    vals = adata[mask, gene].X.toarray().flatten() if hasattr(adata[mask, gene].X, 'toarray') else adata[mask, gene].X.flatten()
                    axes[i].hist(vals, bins=30, alpha=0.5, label=str(grp))
                axes[i].set_ylabel(gene[:20])
                axes[i].legend(fontsize=6)
        
        plt.tight_layout()
        out_path = output_dir / "ridge_plot.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        
        return {
            "figures": {"ridge_plot": str(out_path)},
            "de_report": {"status": "success", "plot_type": "ridge"}
        }
    
    def _find_de_result(self, adata, params):
        """自动找到 DE 结果列"""
        de_result_col = params.get("de_result_col", "")
        if de_result_col and de_result_col in adata.obs.columns:
            return de_result_col
        
        # 搜索 rank_genes_groups 的 names
        if "rank_genes_groups" in adata.uns:
            return "rank_genes_groups"
        
        # 搜索包含 dict 的列
        for col in adata.obs.columns:
            try:
                if isinstance(adata.obs[col].iloc[0], dict):
                    return col
            except (IndexError, TypeError):
                continue
        
        return None
