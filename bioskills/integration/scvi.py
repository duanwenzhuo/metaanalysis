"""
scVI / scANVI Integration Skill — 深度学习批次校正与整合

特点:
- scVI: Variational Inference for single-cell
- scANVI: 对标注数据进行半监督整合（transfer learning）
- 无需假设批次效应线性，能处理复杂效应
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
from typing import Dict, Any, Optional
import numpy as np


@register
class ScVISkill(AbstractSkill):
    """
    L1: scVI Integration Skill
    
    使用 scVI 进行批次校正和整合。
    适用于大规模数据，能学习非线性批次效应。
    
    输入:
      - adata: AnnData — expression matrix
      - batch_key: str — adata.obs 中的批次字段名
      - n_latent: int — 隐空间维度（默认 10）
      - n_layers: int — encoder/decoder 层数（默认 1）
      - max_epochs: int — 最大训练轮数（默认 100）
    
    输出:
      - adata: AnnData — 包含 X_scVI (latent embedding)
      - scvi_model: scvi.model.SCVI — 训练好的模型
      - latent: np.ndarray — 隐表示矩阵
    """
    name = "scvi"
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.SPATIAL]
    
    input_contract = {
        "adata": "anndata.AnnData",
        "batch_key": str,
    }
    output_contract = {
        "adata": "anndata.AnnData",
        "latent": np.ndarray,
    }
    
    required_state_keys = {"adata", "batch_key"}
    
    def _run(self, state: State) -> Dict[str, Any]:
        import anndata
        import numpy as np
        
        adata = state["adata"]
        batch_key = state.get("batch_key", "batch")
        n_latent = state.get("n_latent", 10)
        n_layers = state.get("n_layers", 1)
        max_epochs = state.get("max_epochs", 100)
        
        # 检查 batch_key 存在
        if batch_key not in adata.obs.columns:
            raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")
        
        # 尝试导入 scvi
        try:
            import scvi
        except ImportError:
            # 如果没有安装 scvi，返回降级结果
            print("⚠️  scvi not installed, skipping scVI integration")
            return {
                "adata": adata,
                "latent": np.zeros((adata.n_obs, n_latent)),
                "error": "scvi_not_installed",
            }
        
        # 设置 scVI
        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
        
        model = scvi.model.SCVI(
            adata,
            n_latent=n_latent,
            n_layers=n_layers,
        )
        
        # 训练
        model.train(max_epochs=max_epochs)
        
        # 获取隐表示
        latent = model.get_latent_representation()
        adata.obsm["X_scVI"] = latent
        
        return {
            "adata": adata,
            "latent": latent,
            "scvi_model": model,
        }


@register
class ScanVISkill(AbstractSkill):
    """
    L1: scANVI Integration Skill
    
    使用 scANVI 进行半监督批次整合（考虑细胞类型标注）。
    适用于有部分标注的数据，能实现跨数据集 transfer learning。
    
    输入:
      - adata: AnnData
      - batch_key: str
      - labels_key: str — 细胞类型标注字段（未知用 "Unknown" 填充）
      - n_latent: int
    
    输出:
      - adata: AnnData — 包含 X_scANVI + 预测标签
      - latent: np.ndarray
      - predicted_labels: np.ndarray
    """
    name = "scanvi"
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.SPATIAL]
    
    input_contract = {
        "adata": "anndata.AnnData",
        "batch_key": str,
        "labels_key": str,
    }
    output_contract = {
        "adata": "anndata.AnnData",
        "latent": np.ndarray,
        "predicted_labels": np.ndarray,
    }
    
    required_state_keys = {"adata", "batch_key", "labels_key"}
    
    def _run(self, state: State) -> Dict[str, Any]:
        import anndata
        import numpy as np
        
        adata = state["adata"]
        batch_key = state.get("batch_key", "batch")
        labels_key = state.get("labels_key", "cell_type")
        n_latent = state.get("n_latent", 10)
        max_epochs = state.get("max_epochs", 100)
        
        if batch_key not in adata.obs.columns:
            raise ValueError(f"batch_key '{batch_key}' not found")
        if labels_key not in adata.obs.columns:
            raise ValueError(f"labels_key '{labels_key}' not found")
        
        try:
            import scvi
        except ImportError:
            print("⚠️  scvi not installed, skipping scANVI")
            return {
                "adata": adata,
                "latent": np.zeros((adata.n_obs, n_latent)),
                "predicted_labels": np.array(["Unknown"] * adata.n_obs),
                "error": "scvi_not_installed",
            }
        
        # 先训练 scVI
        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
        scvi_model = scvi.model.SCVI(adata, n_latent=n_latent)
        scvi_model.train(max_epochs=max_epochs)
        
        # 再训练 scANVI
        scvi.model.SCANVI.setup_anndata(
            adata,
            batch_key=batch_key,
            labels_key=labels_key,
            unlabeled_category="Unknown",
        )
        scanvi_model = scvi.model.SCANVI(adata, n_latent=n_latent)
        scanvi_model.train(max_epochs=max_epochs)
        
        # 获取结果
        latent = scanvi_model.get_latent_representation()
        predicted = scanvi_model.predict()
        
        adata.obsm["X_scANVI"] = latent
        adata.obs["predicted_cell_type"] = predicted
        
        return {
            "adata": adata,
            "latent": latent,
            "predicted_labels": predicted,
            "scanvi_model": scanvi_model,
        }
