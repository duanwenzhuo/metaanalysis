"""
PAGA Trajectory Skill — 轨迹结构可视化与推断

PAGA (Partition-based Graph Abstraction) 用于:
- 推断细胞类型的拓扑连接关系
- 可视化发育轨迹
- 检测分支点
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
from typing import Dict, Any, Optional
import numpy as np


@register
class PAGASkill(AbstractSkill):
    """
    L1: PAGA Trajectory Skill
    
    使用 PAGA 分析细胞类型的拓扑结构。
    需要先完成聚类（leiden/louvain）和 neighbors 计算。
    
    输入:
      - adata: AnnData
      - groups: str — 分组字段（默认 "clusters"）
    
    输出:
      - adata: AnnData — 包含 paga 连接信息
      - paga_graph: dict — PAGA 图结构
      - transitions: DataFrame — 细胞类型转移概率
    
    依赖:
      - neighbors
      - clustering（leiden/louvain）
    """
    name = "paga"
    stage = Stage.ANNOTATION
    modality = [Modality.SCRNA, Modality.SPATIAL]
    
    input_contract = {
        "adata": "anndata.AnnData",
        "groups": str,
    }
    output_contract = {
        "adata": "anndata.AnnData",
        "paga_graph": dict,
        "transitions": "pandas.DataFrame",
    }
    
    required_state_keys = {"adata"}
    
    def _run(self, state: State) -> Dict[str, Any]:
        import scanpy as sc
        import pandas as pd
        import numpy as np
        
        adata = state["adata"]
        groups = state.get("groups", "clusters")
        
        # 检查是否有 neighbors
        if "neighbors" not in adata.uns:
            print("⚠️  neighbors not found, computing...")
            sc.pp.neighbors(adata)
        
        # 检查分组字段
        if groups not in adata.obs.columns:
            raise ValueError(f"groups '{groups}' not found in adata.obs")
        
        # 运行 PAGA
        sc.tl.paga(adata, groups=groups)
        
        # 提取图结构
        paga_graph = {
            "connectivities": adata.uns["paga"]["connectivities"].toarray().tolist(),
            "connectivities_tree": adata.uns["paga"]["connectivities_tree"].toarray().tolist() if "connectivities_tree" in adata.uns["paga"] else None,
            "n_nodes": len(adata.obs[groups].unique()),
            "node_names": list(adata.obs[groups].unique()),
        }
        
        # 计算转移概率（基于 diffusion pseudotime）
        try:
            sc.tl.diffmap(adata)
            sc.tl.dpt(adata)
            
            transitions = pd.DataFrame(
                paga_graph["connectivities"],
                index=paga_graph["node_names"],
                columns=paga_graph["node_names"],
            )
        except Exception as e:
            transitions = pd.DataFrame()
            print(f"⚠️  Diffusion pseudotime failed: {e}")
        
        # 生成 PAGA 坐标
        sc.tl.draw_graph(adata, init_pos="paga")
        
        return {
            "adata": adata,
            "paga_graph": paga_graph,
            "transitions": transitions,
        }


@register
class PAGAForkSkill(AbstractSkill):
    """
    L1: PAGA Fork Detection Skill
    
    检测发育轨迹中的分支点。
    基于 PAGA 图的度中心性识别分支。
    
    输入:
      - adata: AnnData（需已运行 PAGA）
    
    输出:
      - forks: List[str] — 分支节点ID
      - endpoints: List[str] — 终端节点ID
      - root: str — 推断的根节点
    """
    name = "paga_fork"
    stage = Stage.ANNOTATION
    modality = [Modality.SCRNA]
    
    input_contract = {"adata": "anndata.AnnData"}
    output_contract = {
        "forks": list,
        "endpoints": list,
        "root": str,
    }
    
    required_state_keys = {"adata"}
    
    def _run(self, state: State) -> Dict[str, Any]:
        import scanpy as sc
        import numpy as np
        import networkx as nx
        
        adata = state["adata"]
        
        if "paga" not in adata.uns:
            raise ValueError("PAGA not computed. Run PAGA skill first.")
        
        # 从连接矩阵构建图
        conn = adata.uns["paga"]["connectivities"]
        G = nx.from_scipy_sparse_array(conn)
        
        # 计算度中心性
        degrees = dict(G.degree())
        
        # 找分支点（度 > 2）
        forks = [n for n, d in degrees.items() if d > 2]
        endpoints = [n for n, d in degrees.items() if d == 1]
        
        # 推断根节点（最高度 + 最高 pagerank）
        try:
            pagerank = nx.pagerank(G)
            root = max(pagerank, key=pagerank.get)
        except:
            root = max(degrees, key=degrees.get)
        
        # 转换为节点名（如果有）
        node_names = list(adata.obs["clusters"].unique()) if "clusters" in adata.obs.columns else None
        
        return {
            "forks": [node_names[f] if node_names else str(f) for f in forks] if forks else [],
            "endpoints": [node_names[e] if node_names else str(e) for e in endpoints] if endpoints else [],
            "root": node_names[root] if node_names else str(root),
        }
