"""L0: Leiden Clustering Skill"""
from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class LeidenSkill(AbstractSkill):
    name = "leiden"
    description = "Leiden community detection clustering (recommended over Louvain)"
    input_contract = ["adata"]
    output_contract = ["adata", "clusters", "n_clusters"]
    stage = Stage.CLUSTER
    modality = [Modality.SCRNA, Modality.SPATIAL, Modality.GENERIC]
    
    tunable_parameters = {
        "resolution": {"type": "float", "default": 0.5, "min": 0.1, "max": 2.0},
        "key": {"type": "str", "default": "clusters"},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {})
        resolution = params.get("resolution", 0.5)
        key = params.get("key", "clusters")
        
        if "neighbors" not in adata.uns:
            sc.pp.neighbors(adata, n_neighbors=15)
        
        sc.tl.leiden(adata, resolution=resolution, key_added=key)
        
        n_clusters = int(adata.obs[key].cat.categories.size)
        clusters = adata.obs[key].to_dict()
        
        print(f"  Leiden: {n_clusters} clusters (resolution={resolution})")
        return {"adata": adata, "clusters": clusters, "n_clusters": n_clusters}
