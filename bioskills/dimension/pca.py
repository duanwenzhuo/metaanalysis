"""L0: PCA Skill"""
from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class PCASkill(AbstractSkill):
    name = "pca"
    description = "Principal Component Analysis (PCA) for dimensionality reduction"
    input_contract = ["adata"]
    output_contract = ["adata", "pca_variance"]
    stage = Stage.DIMENSION
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.GENERIC]
    
    tunable_parameters = {
        "n_comps": {"type": "int", "default": 50, "min": 10, "max": 100},
        "svd_solver": {"type": "str", "default": "arpack"},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        import numpy as np
        adata = state["adata"].copy()
        params = state.get("params", {})
        n_comps = params.get("n_comps", 50)
        
        sc.tl.pca(adata, n_comps=n_comps, svd_solver=params.get("svd_solver", "arpack"))
        
        var = adata.uns["pca"]["variance"]
        var_ratio = var / var.sum() if var.sum() > 0 else var
        cumsum = np.cumsum(var_ratio)
        
        pca_variance = {
            "variance": var.tolist(),
            "variance_ratio": var_ratio.tolist(),
            "cumulative_variance_ratio": cumsum.tolist(),
            "n_comps": n_comps,
            "elbow_component": int(np.argmax(cumsum >= 0.9)) + 1 if len(cumsum) > 0 else n_comps,
        }
        
        print(f"  PCA: {n_comps} components, elbow@{{pca_variance['elbow_component']}} "
              f"(90% var)")
        return {"adata": adata, "pca_variance": pca_variance}
