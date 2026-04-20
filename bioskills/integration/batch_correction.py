"""
L0: BatchCorrection Skill — 批次校正（Harmony / scVI / ComBat-Seq）

bioSkills 原版: bio-differential-expression-batch-correction + bio-single-cell-batch-integration
翻译重构: 适配 bioskills State Dict 契约系统

功能：
- Harmony（Python/R）：快速批次校正
- scVI（Python）：深度生成模型批次校正
- ComBat-Seq（R）：负二项回归批次校正（保留 counts）
- Seurat CCA/RPCA（R）：锚点方法

输出:
  - adata: obsm['X_harmony'] / obsm['X_scVI'] / layers['corrected']
  - batch_report: {method, n_batches, batch_metric}
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State


@register
class BatchCorrectionSkill(AbstractSkill):
    """
    Remove batch effects from merged scRNA-seq datasets.
    
    工具选择策略（bioSkills）：
    | Tool        | Speed  | Scalability | Best For                |
    |-------------|--------|-------------|-------------------------|
    | Harmony     | Fast   | Good        | Quick integration        |
    | scVI        | Slow   | Excellent   | Large datasets (>500k)  |
    | ComBat-Seq  | Fast   | Good        | Count data (RNA-seq)    |
    | Seurat CCA  | Slow   | Good        | Conserved biology       |
    | fastMNN     | Fast   | Good        | MNN-based correction    |
    
    bioSkills 触发词："Remove batch effects from my data"
    """
    
    name = "batch_correction"
    description = (
        "Remove batch effects using Harmony, scVI, ComBat-Seq, or Seurat CCA. "
        "bioSkills: bio-batch-integration + bio-batch-correction"
    )
    input_contract = ["adata"]
    output_contract = ["adata", "batch_report"]
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.BULKRNA]
    
    tunable_parameters = {
        "method": {"type": "str", "default": "harmony"},
        "batch_key": {"type": "str", "default": "batch"},
        "n_harmony_dims": {"type": "int", "default": 30},
        "max_iter_harmony": {"type": "int", "default": 20},
        "scvi_n_latent": {"type": "int", "default": 30},
        "scvi_n_layers": {"type": "int", "default": 2},
        "scvi_max_epochs": {"type": "int", "default": 100},
        "continuous_covariates": {"type": "list", "default": []},
    }
    
    def _run(self, state: State) -> dict:
        params = state.get("params", {})
        method = params.get("method", "harmony")
        
        if method == "harmony":
            return self._run_harmony(state)
        elif method == "scvi":
            return self._run_scvi(state)
        elif method == "combat_seq":
            return self._run_combat_seq(state)
        elif method == "seurat_cca":
            return self._run_seurat_cca(state)
        else:
            return {
                "adata": state["adata"],
                "batch_report": {
                    "status": "failed",
                    "error": f"Unknown method: {method}. "
                             "Use: harmony, scvi, combat_seq, seurat_cca"
                }
            }
    
    def _run_harmony(self, state: State) -> dict:
        """Harmony: 快速迭代批次校正（推荐首选）"""
        import scanpy as sc
        import numpy as np
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        batch_key = params.get("batch_key", "batch")
        n_dims = params.get("n_harmony_dims", 30)
        
        if batch_key not in adata.obs.columns:
            return {
                "adata": adata,
                "batch_report": {
                    "status": "failed",
                    "error": f"Batch key '{batch_key}' not found in adata.obs. "
                             "Available columns: {list(adata.obs.columns)}"
                }
            }
        
        n_batches = adata.obs[batch_key].nunique()
        print(f"  [Harmony] Correcting {n_batches} batches, {adata.n_obs} cells")
        
        # ── Step 1: 预处理 ──────────────────────
        adata_hvg = adata.copy()
        sc.pp.normalize_total(adata_hvg, target_sum=1e4)
        sc.pp.log1p(adata_hvg)
        sc.pp.highly_variable_genes(
            adata_hvg,
            batch_key=batch_key,
            n_top_genes=2000,
        )
        
        # 检查 HVG 数量
        n_hvg = int(adata_hvg.var.highly_variable.sum())
        if n_hvg < 50:
            print(f"  [Harmony] ⚠️  Only {n_hvg} HVGs — batch correction may be weak")
        
        adata_hvg = adata_hvg[:, adata_hvg.var.highly_variable].copy()
        sc.pp.scale(adata_hvg)
        sc.tl.pca(adata_hvg, n_comps=min(n_dims, adata_hvg.n_vars-1))
        
        # ── Step 2: Harmony 校正 ─────────────────
        try:
            import scanpy.external as sce
        except ImportError:
            print("  [Harmony] ⚠️  Installing harmony Python package...")
            import subprocess
            subprocess.run(["pip3", "install", "harmony", "-q"], check=True)
            import scanpy.external as sce
        
        try:
            sce.pp.harmony_integrate(
                adata_hvg,
                key=batch_key,
                max_iter_harmony=params.get("max_iter_harmony", 20),
            )
            
            # 复制结果到原始 adata
            if "X_pca_harmony" in adata_hvg.obsm:
                adata.obsm["X_pca_harmony"] = adata_hvg.obsm["X_pca_harmony"]
                adata.obsm["X_pca"] = adata_hvg.obsm["X_pca_harmony"]  # 替换默认 PCA
            
            print(f"  [Harmony] ✅ Batch correction complete")
            print(f"  [Harmony] Corrected embedding in adata.obsm['X_pca_harmony']")
            
            # ── Step 3: 评估 ────────────────────
            batch_metrics = self._evaluate_integration(adata, batch_key, "X_pca_harmony")
            
            return {
                "adata": adata,
                "batch_report": {
                    "method": "Harmony",
                    "batch_key": batch_key,
                    "n_batches": n_batches,
                    "n_hvg_used": n_hvg,
                    "embedding_key": "X_pca_harmony",
                    **batch_metrics,
                    "status": "success",
                }
            }
            
        except Exception as e:
            return {
                "adata": adata,
                "batch_report": {
                    "status": "failed",
                    "error": f"Harmony failed: {e}",
                    "suggestion": "Ensure harmony R package is installed: "
                                "R -e 'BiocManager::install(\"harmony\")'"
                }
            }
    
    def _run_scvi(self, state: State) -> dict:
        """scVI: 深度生成模型批次校正"""
        import scanpy as sc
        import numpy as np
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        batch_key = params.get("batch_key", "batch")
        n_latent = params.get("scvi_n_latent", 30)
        n_layers = params.get("scvi_n_layers", 2)
        
        try:
            import scvi
        except ImportError:
            print("  [scVI] ⚠️  Installing scvi-tools...")
            import subprocess
            subprocess.run(["pip3", "install", "scvi-tools", "-q"], check=True)
            import scvi
        
        try:
            # Setup scVI
            scvi.model.SCVI.setup_anndata(
                adata,
                batch_key=batch_key,
                continuous_covariate_keys=params.get("continuous_covariates", []),
            )
            
            # Train model
            model = scvi.model.SCVI(
                adata,
                n_latent=n_latent,
                n_layers=n_layers,
            )
            model.train(
                max_epochs=params.get("scvi_max_epochs", 100),
                early_stopping=True,
            )
            
            # 获取 latent representation
            adata.obsm["X_scVI"] = model.get_latent_representation()
            adata.obsm["X_pca"] = adata.obsm["X_scVI"]  # 替换默认
            
            print(f"  [scVI] ✅ Trained (latent dim={n_latent})")
            print(f"  [scVI] Embedding in adata.obsm['X_scVI']")
            
            batch_metrics = self._evaluate_integration(adata, batch_key, "X_scVI")
            
            return {
                "adata": adata,
                "batch_report": {
                    "method": "scVI",
                    "batch_key": batch_key,
                    "n_latent": n_latent,
                    "n_layers": n_layers,
                    "embedding_key": "X_scVI",
                    **batch_metrics,
                    "status": "success",
                }
            }
            
        except Exception as e:
            return {
                "adata": adata,
                "batch_report": {
                    "status": "failed",
                    "error": f"scVI failed: {e}"
                }
            }
    
    def _run_combat_seq(self, state: State) -> dict:
        """ComBat-Seq: 负二项回归批次校正（保留 counts）"""
        import scanpy as sc
        import numpy as np
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        batch_key = params.get("batch_key", "batch")
        
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            pandas2ri.activate()
        except ImportError:
            return {
                "adata": adata,
                "batch_report": {
                    "status": "failed",
                    "error": "rpy2 not installed"
                }
            }
        
        try:
            import tempfile, os
            tmpdir = tempfile.mkdtemp()
            h5ad_path = os.path.join(tmpdir, "adata_temp.h5ad")
            adata.write_h5ad(h5ad_path)
            
            condition_key = params.get("condition_key", "")
            full_mod = "TRUE" if condition_key else "FALSE"
            cond_arg = f", group={condition_key}" if condition_key else ""
            
            r_script = f"""
            library(sva)
            library(SeuratDisk)
            
            adata <- ReadH5AD('{h5ad_path}')
            counts <- as.matrix(adata@assays@data$counts)
            batch <- adata@meta.data[['{batch_key}']]
            {f"condition <- adata@meta.data[['{condition_key}']]" if condition_key else ""}
            
            # ComBat-Seq (requires counts, adjusts for biological groups)
            corrected <- ComBat_seq(
                counts = counts,
                batch = batch,
                group = {f'condition' if condition_key else 'NULL'},
                full_mod = {full_mod}
            )
            
            # Save corrected counts
            corrected_matrix <- as(combat_corrected, 'dgCMatrix')
            saveRDS(corrected, '{tmpdir}/corrected_counts.rds')
            """
            
            if condition_key:
                r_script = r_script.replace("combat_corrected", "corrected")
            else:
                r_script = r_script.replace("combat_corrected", "corrected")
            
            ro.r(r_script)
            
            from rpy2.robjects.vectors import IntVector, FloatVector
            corrected_r = ro.r(f"readRDS('{tmpdir}/corrected_counts.rds')")
            
            # 写回 adata
            from scipy import sparse
            if sparse.issparse(adata.X):
                adata.X = sparse.csr_matrix(corrected_r)
            else:
                adata.X = np.array(corrected_r)
            
            print(f"  [ComBat-Seq] ✅ Batch correction complete")
            
            return {
                "adata": adata,
                "batch_report": {
                    "method": "ComBat-Seq",
                    "batch_key": batch_key,
                    "preserves_counts": True,
                    "status": "success",
                }
            }
            
        except Exception as e:
            return {
                "adata": adata,
                "batch_report": {
                    "status": "failed",
                    "error": f"ComBat-Seq failed: {e}"
                }
            }
    
    def _run_seurat_cca(self, state: State) -> dict:
        """Seurat CCA: 锚点方法批次校正"""
        return {
            "adata": state["adata"],
            "batch_report": {
                "status": "skipped",
                "method": "Seurat CCA",
                "error": "Seurat CCA requires R. Use harmony or scVI (Python) instead.",
            }
        }
    
    def _evaluate_integration(self, adata, batch_key, embedding_key):
        """评估批次校正效果：计算 silhouette scores"""
        from sklearn.metrics import silhouette_score
        import numpy as np
        
        try:
            if embedding_key not in adata.obsm:
                return {"integration_evaluation": "skipped (no embedding)"}
            
            emb = adata.obsm[embedding_key]
            batch_labels = adata.obs[batch_key].cat.codes.values
            
            # Batch silhouette (低 = 批次混合好)
            batch_sil = silhouette_score(emb, batch_labels)
            
            # 如果有细胞类型，也计算
            cell_type_col = next(
                (c for c in ["cell_type", "celltype", "cell_type_annotations"] 
                 if c in adata.obs.columns), None
            )
            
            eval_result = {
                "batch_silhouette_score": round(batch_sil, 4),
                "batch_silhouette_interpretation": (
                    "GOOD (well mixed)" if abs(batch_sil) < 0.1 else
                    "FAIR" if abs(batch_sil) < 0.3 else
                    "POOR (batch effect remains)"
                ),
            }
            
            if cell_type_col:
                type_labels = adata.obs[cell_type_col].cat.codes.values
                type_sil = silhouette_score(emb, type_labels)
                eval_result["celltype_silhouette_score"] = round(type_sil, 4)
                print(f"  [Harmony] Batch silhouette: {batch_sil:.3f} "
                      f"(lower=better mixing)")
                print(f"  [Harmony] Cell type silhouette: {type_sil:.3f} "
                      f"(higher=better separation)")
            else:
                print(f"  [Harmony] Batch silhouette: {batch_sil:.3f} "
                      f"(lower=better mixing)")
            
            return eval_result
            
        except Exception as e:
            return {"integration_evaluation": f"error: {e}"}
