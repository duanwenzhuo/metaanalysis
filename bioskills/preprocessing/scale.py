"""
L0: Scale Skill — Z-score 标准化
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class ScaleSkill(AbstractSkill):
    name = "scale"
    description = "Z-score scaling per gene (zero mean, unit variance)"
    input_contract = ["adata"]
    output_contract = ["adata"]
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.GENERIC]
    
    tunable_parameters = {
        "max_value": {"type": "float", "default": 10.0},
        "zero_center": {"type": "bool", "default": True},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        sc.pp.scale(
            adata,
            max_value=params.get("max_value", 10.0),
            zero_center=params.get("zero_center", True),
        )
        
        print(f"  Scale: max_value={params.get('max_value', 10.0)}")
        return {"adata": adata}
