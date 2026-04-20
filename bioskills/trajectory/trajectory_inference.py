"""
L0: TrajectoryInference Skill — 轨迹推断（Monocle3 / scVelo / PAGA）

bioSkills 原版: bio-single-cell-trajectory-inference
翻译重构: 适配 bioskills State Dict 契约系统

功能：
- Monocle3（R）：principal graph 轨迹推断
- scVelo（Python）：RNA velocity 动力学推断
- PAGA（Python）：partition-based graph abstraction

输出:
  - adata: obs 增加 pseudotime 列
  - trajectory_report: {trajectory_type, branch_points, endpoints, pseudotime_range}
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
from pathlib import Path


@register
class TrajectoryInferenceSkill(AbstractSkill):
    """
    Infer developmental trajectories and pseudotime from scRNA-seq.
    
    工具选择策略（bioSkills）：
    | Tool      | Speed | Scalability | Best For |
    |-----------|-------|-------------|----------|
    | Monocle3  | Slow  | Good        | Trajectory + pseudotime |
    | scVelo    | Slow  | Moderate    | RNA velocity (dynamics) |
    | PAGA      | Fast  | Excellent   | Large-scale trajectories |
    | Slingshot | Fast  | Good        | Simple trajectories |
    
    bioSkills 触发词：
    "Find the developmental trajectory in my data"
    "Infer pseudotime ordering"
    """
    
    name = "trajectory_inference"
    description = (
        "Infer developmental trajectories and pseudotime using Monocle3, scVelo, or PAGA. "
        "bioSkills: bio-single-cell-trajectory-inference"
    )
    input_contract = ["adata"]
    output_contract = ["adata", "trajectory_report"]
    stage = Stage.CLUSTER
    modality = [Modality.SCRNA]
    
    tunable_parameters = {
        "method": {"type": "str", "default": "monocle3"},
        "n_dims": {"type": "int", "default": 50},  # num_dim for PCA
        "resolution": {"type": "float", "default": 1e-5},  # for clustering
        "root_cell_col": {"type": "str", "default": ""},  # empty = auto
        "velocity_mode": {"type": "str", "default": "stochastic"},
        "n_neighbors_paga": {"type": "int", "default": 15},
        "clustering_resolution": {"type": "float", "default": 0.5},
    }
    
    def _run(self, state: State) -> dict:
        params = state.get("params", {})
        method = params.get("method", "monocle3")
        
        if method in ["monocle3", "slingshot"]:
            return self._run_monocle3(state)
        elif method in ["scvelo", "velocity", "paga"]:
            return self._run_scvelo(state)
        else:
            return {
                "adata": state["adata"],
                "trajectory_report": {
                    "status": "failed",
                    "error": f"Unknown method: {method}. "
                             "Use 'monocle3', 'slingshot', or 'scvelo'."
                }
            }
    
    def _run_monocle3(self, state: State) -> dict:
        """Monocle3: Trajectory inference via principal graph"""
        import scanpy as sc
        import numpy as np
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        n_dims = params.get("n_dims", 50)
        
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            pandas2ri.activate()
        except ImportError:
            return {
                "adata": adata,
                "trajectory_report": {
                    "status": "failed",
                    "error": "rpy2 not installed. Use scvelo method instead."
                }
            }
        
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                
                # 保存临时文件
                import tempfile, os
                tmpdir = tempfile.mkdtemp()
                h5ad_path = os.path.join(tmpdir, "adata_temp.h5ad")
                adata.write_h5ad(h5ad_path)
                
                root_col = params.get("root_cell_col", "")
                clustering_res = params.get("clustering_resolution", 0.5)
                n_dims = params.get("n_dims", 50)
                
                r_script = f"""
                library(monocle3)
                library(Seurat)
                library(SeuratDisk)
                
                # Read adata
                cds <- as.cell_data_set(Seurat::LoadH5Seurat('{h5ad_path}'))
                
                # Preprocess CDS
                cds <- preprocess_cds(cds, num_dim = {n_dims})
                cds <- reduce_dimension(cds, reduction_method = 'UMAP')
                
                # Cluster
                cds <- cluster_cells(cds, resolution = {clustering_res})
                
                # Learn trajectory graph
                cds <- learn_graph(cds)
                
                # Order cells (use cluster 1 as root if not specified)
                {f"root_cells <- colnames(cds)@colData@listData[['{root_col}']][1]" if root_col else "root_cells <- NULL"}
                cds <- order_cells(cds, root_cells = root_cells)
                
                # Get pseudotime
                pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime
                names(pseudotime) <- colnames(cds)
                
                # Get branch points and endpoints
                br_pts <- cds@principal_graph_aux$UMAP$branch_points
                end_pts <- cds@principal_graph_aux$UMAP$end_points
                
                # Write results
                saveRDS(pseudotime, '{tmpdir}/pseudotime.rds')
                saveRDS(br_pts, '{tmpdir}/branch_points.rds')
                saveRDS(end_pts, '{tmpdir}/endpoints.rds')
                """
                
                ro.r(r_script)
                
                # 读取结果
                pseudotime_r = ro.r(f"readRDS('{tmpdir}/pseudotime.rds')")
                pseudotime = np.array(pseudotime_r)
                
                if len(pseudotime) == adata.n_obs:
                    adata.obs["pseudotime_monocle3"] = pseudotime
                
                try:
                    br_pts_r = ro.r(f"readRDS('{tmpdir}/branch_points.rds')")
                    branch_points = list(br_pts_r)
                except:
                    branch_points = []
                
                try:
                    end_pts_r = ro.r(f"readRDS('{tmpdir}/endpoints.rds')")
                    endpoints = list(end_pts_r)
                except:
                    endpoints = []
                
                n_ordered = int(np.sum(~np.isnan(pseudotime)))
                print(f"  [Monocle3] ✅ Pseudotime computed for {n_ordered}/{adata.n_obs} cells")
                print(f"  [Monocle3] Branch points: {branch_points}")
                print(f"  [Monocle3] Endpoints: {endpoints}")
                
                return {
                    "adata": adata,
                    "trajectory_report": {
                        "method": "Monocle3",
                        "n_cells_ordered": n_ordered,
                        "branch_points": branch_points,
                        "endpoints": endpoints,
                        "pseudotime_range": (
                            float(np.nanmin(pseudotime)),
                            float(np.nanmax(pseudotime))
                        ) if len(pseudotime) > 0 else (None, None),
                        "status": "success",
                    }
                }
                
        except Exception as e:
            return {
                "adata": state["adata"],
                "trajectory_report": {
                    "status": "failed",
                    "error": f"Monocle3 R failed: {e}",
                    "suggestion": "Install Monocle3: R -e 'BiocManager::install(c('monocle3'))'"
                }
            }
    
    def _run_scvelo(self, state: State) -> dict:
        """scVelo: RNA velocity dynamics inference"""
        import scanpy as sc
        import numpy as np
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        try:
            import scvelo as scv
        except ImportError:
            print("  [scVelo] ⚠️  Installing scvelo...")
            import subprocess
            subprocess.run(["pip3", "install", "scvelo", "-q"], check=True)
            import scvelo as scv
        
        try:
            import scvelo as scv
            
            # 1. 检查是否有 velocity layers (spliced/unspliced)
            if "Ms" not in adata.layers and "spliced" not in adata.layers:
                return {
                    "adata": adata,
                    "trajectory_report": {
                        "status": "failed",
                        "error": "scVelo requires spliced/unspliced counts. "
                                 "Ensure your data has adata.layers['spliced'] and adata.layers['unspliced']."
                    }
                }
            
            # 2. Velocity estimation
            scv.settings.set_figure_params("scvelo")
            
            scv.pp.filter_genes(adata, min_shared_counts=20)
            scv.pp.normalize_per_cell(adata)
            scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
            scv.pp.highly_variable_genes(adata)
            scv.pp.moments(adata, n_neighbors=30, n_pcs=30)
            
            velocity_mode = params.get("velocity_mode", "stochastic")
            
            if velocity_mode == "stochastic":
                scv.tl.velocity(adata, mode="stochastic")
            elif velocity_mode == "deterministic":
                scv.tl.velocity(adata, mode="deterministic")
            elif velocity_mode == "steady_state":
                scv.tl.velocity(adata, mode="steady_state")
            
            scv.tl.velocity_graph(adata, n_jobs=4)
            
            # 3. Pseudotime via velocity
            scv.tl.velocity_pseudotime(adata)
            
            # 4. PAGA trajectory (if UMAP already exists or compute it)
            if "X_umap" not in adata.obsm:
                scv.tl.umap(adata)
            
            try:
                scv.tl.paga(adata, groups="leiden" if "leiden" in adata.obs.columns else None)
            except Exception:
                pass
            
            # 写回 adata.obs
            if "velocity_pseudotime" in adata.obs.columns:
                adata.obs["pseudotime_scvelo"] = adata.obs["velocity_pseudotime"]
            
            print(f"  [scVelo] ✅ Velocity computed ({velocity_mode} mode)")
            print(f"  [scVelo] Velocity graph: {adata.n_obs} cells × {adata.n_vars} genes")
            
            return {
                "adata": adata,
                "trajectory_report": {
                    "method": f"scVelo ({velocity_mode})",
                    "n_cells": adata.n_obs,
                    "n_velocity_genes": adata.n_vars,
                    "has_paga": "paga_connectivities" in adata.uns,
                    "pseudotime_key": "velocity_pseudotime",
                    "status": "success",
                }
            }
            
        except Exception as e:
            return {
                "adata": state["adata"],
                "trajectory_report": {
                    "status": "failed",
                    "error": f"scVelo failed: {e}"
                }
            }
