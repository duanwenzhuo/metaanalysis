"""
L0: Normalize Skill — 归一化

输入: adata
输出: adata (normalized + log1p)
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class NormalizeSkill(AbstractSkill):
    name = "normalize"
    description = "Total count normalization + log1p transformation"
    input_contract = ["adata"]
    output_contract = ["adata"]
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.GENERIC]
    
    tunable_parameters = {
        "normalize_target_sum": {"type": "float", "default": 1e4},
        "exclude_highly_expressed": {"type": "bool", "default": False},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        target_sum = params.get("normalize_target_sum", 1e4)
        exclude_high = params.get("exclude_highly_expressed", False)
        
        if exclude_high:
            sc.pp.normalize_total(adata, target_sum=target_sum,
                                  exclude_highly_expressed=True)
        else:
            sc.pp.normalize_total(adata, target_sum=target_sum)
        
        sc.pp.log1p(adata)
        print(f"  Normalize: target_sum={target_sum}, log1p applied")
        
        return {"adata": adata}
