"""
SCENIC Gene Regulatory Network Skill — 基因调控网络推断

SCENIC (Single-Cell rEglatory Network Inference and Clustering) 用于:
- 推断转录因子 (TF) 的调控网络
- 计算 regulon 活性评分 (AUCell)
- 基于调控网络的细胞聚类
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
from typing import Dict, Any, List, Optional
import numpy as np


@register
class SCENICSkill(AbstractSkill):
    """
    L1: SCENIC Gene Regulatory Network Skill
    
    推断基因调控网络 (GRN) 并计算 regulon 活性。
    
    输入:
      - adata: AnnData
      - organism: str — "human" | "mouse" | "fly"
      - tf_list: List[str] — 自定义 TF 列表（可选）
    
    输出:
      - adata: AnnData — 包含 regulon 活性评分
      - regulons: dict — {regulon_name: target_genes}
      - auc_matrix: DataFrame — cells × regulons AUCell 评分
    
    注意:
      - 需要安装 pySCENIC: pip install pyscenic
      - 需要 motif 数据库和 TF-target 数据库
    """
    name = "scenic"
    stage = Stage.ENRICHMENT
    modality = [Modality.SCRNA]
    
    input_contract = {
        "adata": "anndata.AnnData",
        "organism": str,
    }
    output_contract = {
        "adata": "anndata.AnnData",
        "regulons": dict,
        "auc_matrix": "pandas.DataFrame",
    }
    
    required_state_keys = {"adata"}
    
    def _run(self, state: State) -> Dict[str, Any]:
        import anndata
        import pandas as pd
        import numpy as np
        
        adata = state["adata"]
        organism = state.get("organism", "human")
        tf_list = state.get("tf_list", None)
        n_cpu = state.get("n_cpu", 4)
        
        # 检查 pySCENIC 是否安装
        try:
            from pyscenic.utils import modules_from_adjacencies
            from pyscenic.prune import prune2df
            from pyscenic.aucell import aucell
            from pyscenic.scores import score_cells
        except ImportError:
            print("⚠️  pyscenic not installed, returning empty GRN")
            return {
                "adata": adata,
                "regulons": {},
                "auc_matrix": pd.DataFrame(),
                "error": "pyscenic_not_installed",
            }
        
        # SCENIC 三步流程
        # Step 1: GRN inference (GENIE3 / GRNBoost2)
        try:
            from arboreto.utils import load_tf_names
            from arboreto.algo import grnboost2
            
            # 获取 TF 列表
            if tf_list is None:
                tf_names = load_tf_names(f"https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/tf_names_{organism}.txt")
            else:
                tf_names = tf_list
            
            # 运行 GRNBoost2
            # 注意：这是简化版本，完整 SCENIC 需要更多数据库
            adjacencies = grnboost2(
                expression_matrix=adata.to_df().T,# genes x cells
                tf_names=tf_names,
                verbose=True,
                client_or_address="local",  # 本地计算
                n_estimators=100,
            )
            
            # 构建 modules
            modules = modules_from_adjacencies(
                adjacencies,
                adata.to_df().T,
            )
        except Exception as e:
            print(f"⚠️  GRN inference failed: {e}")
            modules = []
        
        # Step 2: Motif pruning (cisTarget)
        # 简化：跳过完整 motif 分析，直接用 TF-target 边
        
        regulons = {}
        for module in modules:
            tf_name = module.transcription_factor
            target_genes = list(module.genes)
            regulons[f"{tf_name}(+)"] = target_genes
        
        # Step 3: AUCell scoring
        if regulons:
            try:
                # 简化版 AUCell
                from pyscenic.aucell import aucell as aucell_func
                auc_mtx = aucell_func(
                    adata.to_df().T,# genes × cells
                    regulons,
                ).T# regulons × cells
            except Exception as e:
                print(f"⚠️  AUCell failed: {e}")
                # 手动计算简化版 AUCell
                auc_mtx = pd.DataFrame(
                    index=regulons.keys(),
                    columns=adata.obs_names,
                )
                for regulon, genes in regulons.items():
                    # 计算 regulon 中基因的平均表达
                    available_genes = [g for g in genes if g in adata.var_names]
                    if available_genes:
                        auc_mtx.loc[regulon] = adata[:, available_genes].X.mean(axis=1)
            else:
                auc_mtx = auc_mtx.T
            
            # 将 regulon 活性添加到 adata.obs
            for regulon in regulons:
                if regulon in auc_mtx.columns:
                    adata.obs[f"regulon_{regulon}"] = auc_mtx[regulon]
                elif regulon in auc_mtx.index:
                    adata.obs[f"regulon_{regulon}"] = auc_mtx.loc[regulon]
        else:
            auc_mtx = pd.DataFrame()
        
        return {
            "adata": adata,
            "regulons": regulons,
            "auc_matrix": auc_mtx,
        }


@register
class RegulonScoreSkill(AbstractSkill):
    """
    L1: Regulon Scoring Skill
    
    为给定的 regulon 计算 AUCell 评分。
    简化版，不需要完整 SCENIC 流程。
    
    输入:
      - adata: AnnData
      - regulon_genes: List[str] — 目标基因列表
    
    输出:
      - adata: AnnData — 包含 regulon 评分
      - score: np.ndarray
    """
    name = "regulon_score"
    stage = Stage.ENRICHMENT
    modality = [Modality.SCRNA]
    
    input_contract = {
        "adata": "anndata.AnnData",
        "regulon_genes": list,
    }
    output_contract = {
        "adata": "anndata.AnnData",
        "score": np.ndarray,
    }
    
    required_state_keys = {"adata", "regulon_genes"}
    
    def _run(self, state: State) -> Dict[str, Any]:
        import scanpy as sc
        import numpy as np
        
        adata = state["adata"]
        regulon_genes = state["regulon_genes"]
        regulon_name = state.get("regulon_name", "custom_regulon")
        
        # 过滤可用基因
        available = [g for g in regulon_genes if g in adata.var_names]
        
        if not available:
            print(f"⚠️  No genes from regulon found in adata")
            return {
                "adata": adata,
                "score": np.zeros(adata.n_obs),
            }
        
        # 使用 scanpy score_genes
        sc.tl.score_genes(
            adata,
            gene_list=available,
            score_name=f"regulon_{regulon_name}",
            copy=False,
        )
        
        score = adata.obs[f"regulon_{regulon_name}"].values
        
        return {
            "adata": adata,
            "score": score,
        }
