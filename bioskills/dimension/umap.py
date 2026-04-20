"""L0: UMAP Skill"""
from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class UMAPSkill(AbstractSkill):
    name = "umap"
    description = "UMAP 2D embedding visualization"
    input_contract = ["adata"]
    output_contract = ["adata"]
    stage = Stage.DIMENSION
    modality = [Modality.SCRNA, Modality.SPATIAL, Modality.GENERIC]
    
    tunable_parameters = {
        "min_dist": {"type": "float", "default": 0.5, "min": 0.0, "max": 1.0},
        "n_components": {"type": "int", "default": 2},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        if "neighbors" not in adata.uns:
            sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(50, adata.n_vars-1))
        
        sc.tl.umap(
            adata,
            min_dist=params.get("min_dist", 0.5),
            n_components=params.get("n_components", 2),
        )
        
        print(f"  UMAP: min_dist={params.get('min_dist', 0.5)}")
        return {"adata": adata}
