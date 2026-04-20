"""L0: Neighbors Skill"""
from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class NeighborsSkill(AbstractSkill):
    name = "neighbors"
    description = "Build neighbor graph (required for clustering and UMAP)"
    input_contract = ["adata"]
    output_contract = ["adata"]
    stage = Stage.DIMENSION
    modality = [Modality.SCRNA, Modality.SPATIAL, Modality.GENERIC]
    
    tunable_parameters = {
        "n_neighbors": {"type": "int", "default": 15, "min": 5, "max": 100},
        "n_pcs": {"type": "int", "default": 50},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        n_neighbors = params.get("n_neighbors", 15)
        n_pcs = params.get("n_pcs", 50)
        
        if "X_pca" not in adata.obsm:
            print(f"  [{self.name}] ⚠️  X_pca not found, running PCA first")
            sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1))
        
        actual_pcs = min(n_pcs, adata.obsm["X_pca"].shape[1])
        
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=actual_pcs)
        print(f"  Neighbors: k={n_neighbors}, n_pcs={actual_pcs}")
        return {"adata": adata}
