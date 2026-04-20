"""
L0: CellAnnotationRef Skill — 参考数据库自动细胞注释（CellTypist / ScPred / SingleR）

bioSkills 原版: bio-single-cell-cell-annotation (部分)
翻译重构: 适配 bioskills State Dict 契约系统

功能：
- CellTypist：参考数据库自动注释 + majority voting
- ScPred：基于表达签名的分类
- SingleR：参考驱动精细注释

bioSkills 对比：
| Tool    | Speed | Accuracy | Requires Reference |
|---------|-------|----------|------------------|
| CellTypist | Fast | Good   | Model file      |
| ScPred   | Fast | Fair    | Trained model   |
| SingleR  | Slow | Excellent| Reference data  |
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State


@register
class CellAnnotationRefSkill(AbstractSkill):
    """
    Automated cell type annotation using reference-based methods.
    
    bioSkills Skill: bio-single-cell-cell-annotation
    覆盖: CellTypist, ScPred, SingleR
    """
    
    name = "cell_annotation_ref"
    description = (
        "Automated cell type annotation using CellTypist, ScPred, or SingleR. "
        "bioSkills: bio-single-cell-cell-annotation"
    )
    input_contract = ["adata"]
    output_contract = ["adata", "annotation_report"]
    stage = Stage.ANNOTATION
    modality = [Modality.SCRNA]
    
    tunable_parameters = {
        "method": {"type": "str", "default": "celltypist"},
        "model": {"type": "str", "default": "Human_G1_G2_M_premium"},
        "majority_voting": {"type": "bool", "default": True},
        "min_confidence": {"type": "float", "default": 0.0},
        "ref_label_col": {"type": "str", "default": ""},
    }
    
    def _run(self, state: State) -> dict:
        params = state.get("params", {})
        method = params.get("method", "celltypist")
        
        if method == "celltypist":
            return self._run_celltypist(state)
        elif method == "singleR":
            return self._run_singleR(state)
        elif method == "scpred":
            return self._run_scpred(state)
        else:
            return {
                "adata": state["adata"],
                "annotation_report": {
                    "status": "failed",
                    "error": f"Unknown method: {method}. "
                             "Use: celltypist, singleR, scpred"
                }
            }
    
    def _run_celltypist(self, state: State) -> dict:
        """CellTypist: 参考模型自动注释"""
        import scanpy as sc
        import numpy as np
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        model = params.get("model", "Human_G1_G2_M_premium")
        majority_voting = params.get("majority_voting", True)
        
        try:
            import celltypist
        except ImportError:
            print("  [CellTypist] ⚠️  Installing celltypist...")
            import subprocess
            subprocess.run(["pip3", "install", "celltypist", "-q"], check=True)
            import celltypist
        
        try:
            # 下载模型（如果需要）
            try:
                model_path = celltypist.models.get_model(model)
            except Exception:
                print(f"  [CellTypist] Downloading model '{model}'...")
                model_path = celltypist.models.Model(model)
                model_path.download()
            
            # 预测
            if majority_voting:
                predictions = celltypist.annotate(
                    adata,
                    model=model_path,
                    majority_voting=True,
                    min_confidence=params.get("min_confidence", 0.0),
                )
                # 提取结果
                adata.obs["celltypist_prediction"] = predictions.predicted_labels[
                    "predicted_labels"
                ].values
                adata.obs["celltypist_majority_voting"] = predictions.predicted_labels[
                    "majority_voting"
                ].values
                
                pred_col = "celltypist_majority_voting"
            else:
                predictions = celltypist.annotate(
                    adata, model=model_path, majority_voting=False,
                )
                adata.obs["celltypist_prediction"] = predictions.predicted_labels[
                    "predicted_labels"
                ].values
                pred_col = "celltypist_prediction"
            
            # 统计注释结果
            cell_type_counts = adata.obs[pred_col].value_counts()
            
            print(f"  [CellTypist] ✅ Annotated {adata.n_obs} cells ({len(cell_type_counts)} types)")
            for ct, n in cell_type_counts.head(10).items():
                print(f"    {ct}: {n} cells ({n/adata.n_obs:.1%})")
            
            return {
                "adata": adata,
                "annotation_report": {
                    "method": "CellTypist",
                    "model": model,
                    "n_cell_types": len(cell_type_counts),
                    "cell_type_distribution": cell_type_counts.head(20).to_dict(),
                    "prediction_col": pred_col,
                    "majority_voting": majority_voting,
                    "status": "success",
                }
            }
            
        except Exception as e:
            return {
                "adata": state["adata"],
                "annotation_report": {
                    "status": "failed",
                    "error": f"CellTypist failed: {e}",
                    "suggestion": "Try a different model: "
                                 "celltypist.models.get_model('Human_All_Myeloid')"
                }
            }
    
    def _run_singleR(self, state: State) -> dict:
        """SingleR: Bioconductor 参考驱动注释"""
        import scanpy as sc
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            pandas2ri.activate()
        except ImportError:
            return {
                "adata": adata,
                "annotation_report": {
                    "status": "failed",
                    "error": "rpy2 not installed. Use CellTypist (Python) instead."
                }
            }
        
        try:
            import tempfile, os
            tmpdir = tempfile.mkdtemp()
            h5ad_path = os.path.join(tmpdir, "adata_temp.h5ad")
            adata.write_h5ad(h5ad_path)
            
            r_script = f"""
            library(SingleCellExperiment)
            library(SeuratDisk)
            library(celldex)
            
            adata <- ReadH5AD('{h5ad_path}')
            
            # Use reference datasets from celldex
            ref <- HumanPrimaryCellAtlasData()
            
            # Run SingleR
            pred <- SingleR(
                test = adata,
                ref = ref,
                labels = ref$label.main
            )
            
            # Extract labels
            labels <- as.character(pred$labels)
            scores <- pred@metadata$scores
            
            saveRDS(labels, '{tmpdir}/labels.rds')
            saveRDS(scores, '{tmpdir}/scores.rds')
            """
            
            ro.r(r_script)
            
            labels_r = ro.r(f"readRDS('{tmpdir}/labels.rds')")
            labels = list(labels_r)
            
            adata.obs["singleR_prediction"] = labels[:adata.n_obs]
            
            counts = adata.obs["singleR_prediction"].value_counts()
            print(f"  [SingleR] ✅ {len(counts)} cell types annotated")
            
            return {
                "adata": adata,
                "annotation_report": {
                    "method": "SingleR",
                    "reference": "HumanPrimaryCellAtlasData",
                    "n_cell_types": len(counts),
                    "cell_type_distribution": counts.head(20).to_dict(),
                    "prediction_col": "singleR_prediction",
                    "status": "success",
                }
            }
            
        except Exception as e:
            return {
                "adata": state["adata"],
                "annotation_report": {
                    "status": "failed",
                    "error": f"SingleR failed: {e}"
                }
            }
    
    def _run_scpred(self, state: State) -> dict:
        """ScPred: 基于表达签名分类"""
        import scanpy as sc
        import numpy as np
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        try:
            import scpred as sp
        except ImportError:
            print("  [ScPred] ⚠️  Installing scpred...")
            import subprocess
            subprocess.run(["pip3", "install", "scpred", "-q"], check=True)
            import scpred as sp
        
        ref_label_col = params.get("ref_label_col", "")
        if not ref_label_col:
            ref_label_col = next(
                (c for c in ["cell_type", "celltype", "cell_type_annotations"]
                 if c in adata.obs.columns), None
            )
        
        if not ref_label_col or ref_label_col not in adata.obs.columns:
            return {
                "adata": adata,
                "annotation_report": {
                    "status": "failed",
                    "error": f"Reference label column not found. "
                             "Provide ref_label_col or label a subset of cells first."
                }
            }
        
        # 训练 ScPred 模型
        sp.tl.train(adata, labels=ref_label_col, features="hvg")
        adata.obs["scpred_prediction"] = sp.tl.predict(adata)
        
        counts = adata.obs["scpred_prediction"].value_counts()
        print(f"  [ScPred] ✅ {len(counts)} cell types predicted")
        
        return {
            "adata": adata,
            "annotation_report": {
                "method": "ScPred",
                "ref_label_col": ref_label_col,
                "n_cell_types": len(counts),
                "cell_type_distribution": counts.head(20).to_dict(),
                "prediction_col": "scpred_prediction",
                "status": "success",
            }
        }
