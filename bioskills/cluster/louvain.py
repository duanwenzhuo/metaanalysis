"""L0: Louvain Clustering Skill"""
from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class LouvainSkill(AbstractSkill):
    name = "louvain"
    description = "Louvain community detection clustering"
    input_contract = ["adata"]
    output_contract = ["adata", "clusters", "n_clusters"]
    stage = Stage.CLUSTER
    modality = [Modality.SCRNA, Modality.GENERIC]
    
    tunable_parameters = {
        "resolution": {"type": "float", "default": 0.5, "min": 0.1, "max": 2.0},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {})
        resolution = params.get("resolution", 0.5)
        
        if "neighbors" not in adata.uns:
            sc.pp.neighbors(adata, n_neighbors=15)
        
        sc.tl.louvain(adata, resolution=resolution, key_added="clusters")
        
        n_clusters = int(adata.obs["clusters"].cat.categories.size)
        clusters = adata.obs["clusters"].to_dict()
        
        print(f"  Louvain: {n_clusters} clusters (resolution={resolution})")
        return {"adata": adata, "clusters": clusters, "n_clusters": n_clusters}
