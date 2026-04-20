"""
L0: HVG Skill — 高度可变基因筛选
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class HVGSkill(AbstractSkill):
    name = "hvg"
    description = "Identify highly variable genes (HVG) using Seurat flavor"
    input_contract = ["adata"]
    output_contract = ["adata"]
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.GENERIC]
    
    tunable_parameters = {
        "n_top_genes": {"type": "int", "default": 2000, "min": 500, "max": 5000},
        "flavor": {"type": "str", "default": "seurat"},
        "min_mean": {"type": "float", "default": 0.0125},
        "max_mean": {"type": "float", "default": 3.0},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        n_hvg = params.get("n_top_genes", 2000)
        flavor = params.get("flavor", "seurat")
        
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_hvg,
            flavor=flavor,
            min_mean=params.get("min_mean", 0.0125),
            max_mean=params.get("max_mean", 3.0),
            subset=True,
        )
        
        n_hvg_actual = adata.n_vars
        print(f"  HVG: selected {n_hvg_actual} genes (requested {n_hvg})")
        
        return {"adata": adata}
