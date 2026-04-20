"""L2 Synthesized Skill: custom_gsva_scores
GSVA / Gene Set Enrichment Analysis Skill
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
import numpy as np

@register
class GsvaScoresSkill(AbstractSkill):
    name = "custom_gsva_scores"
    description = "Gene Set Variation Analysis (GSVA) for gene set scoring"
    input_contract = ["adata", "params"]
    output_contract = ["gsva_scores", "adata"]
    stage = Stage.ENRICHMENT
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.GENERIC]

    def _run(self, state: State) -> dict:
        import scanpy as sc
        import gseapy as gp
        adata = state["adata"]
        params = state.get("params", {})
        gene_sets = params.get("gene_sets", {"CUSTOM": []})
        expr = adata.to_df().T
        available_gs = {}
        for gs_name, genes in gene_sets.items():
            avail = [g for g in genes if g in expr.index]
            if len(avail) >= 5:
                available_gs[gs_name] = avail
        if not available_gs:
            return {"gsva_scores": {}, "adata": adata}
        res = gp.ssgsea(
            data=expr, gene_sets=available_gs,
            sample_norm_method=params.get("sample_norm", "rank"),
            min_size=5, max_size=500, no_plot=True, outdir=None,
        )
        scores = res.resultsOnSamples.T
        for gs_name, row in scores.iterrows():
            adata.obs[f"gsva_{gs_name}"] = adata.obs_names.map(row.to_dict()).values
        return {"gsva_scores": scores.to_dict(), "adata": adata}
