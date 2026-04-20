"""L0: MannWhitney Skill — Mann-Whitney U 检验"""
from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
import numpy as np

@register
class MannWhitneySkill(AbstractSkill):
    name = "mann_whitney"
    description = "Mann-Whitney U non-parametric test between two groups"
    input_contract = ["adata"]
    output_contract = ["mann_whitney_result"]
    stage = Stage.STATISTICS
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.GENERIC]
    
    def _run(self, state: State) -> dict:
        from scipy.stats import mannwhitneyu
        adata = state["adata"]
        params = state.get("params", {})
        score_col = params.get("score_column", 
                               [c for c in adata.obs.columns if c.startswith("gsva_")][0]
                               if any(c.startswith("gsva_") for c in adata.obs.columns) else None)
        if not score_col:
            return {"mann_whitney_result": {"status": "skipped", "reason": "no score col"}}
        
        groups = self._infer_groups(adata, params)
        if not groups:
            return {"mann_whitney_result": {"status": "failed", "error": "no groups"}}
        
        a = adata.obs.loc[groups.get("tumor", []), score_col].dropna().values
        b = adata.obs.loc[groups.get("normal", []), score_col].dropna().values
        
        try:
            stat, p = mannwhitneyu(a, b, alternative='two-sided')
            return {
                "mann_whitney_result": {
                    "status": "success",
                    "statistic": float(stat),
                    "p_value": float(p),
                    "n_tumor": len(a), "n_normal": len(b),
                }
            }
        except Exception as e:
            return {"mann_whitney_result": {"status": "failed", "error": str(e)}}
    
    def _infer_groups(self, adata, params):
        if "groups" in params: return params["groups"]
        for col in ["group", "condition"]:
            if col in adata.obs.columns:
                return {str(c).lower(): adata.obs_names[adata.obs[col]==c].tolist()
                         for c in adata.obs[col].cat.categories}
        return None
