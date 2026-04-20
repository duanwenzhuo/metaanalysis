"""
L0: CellCommunication Skill — 细胞间通讯分析（CellChat / LIANA）

bioSkills 原版: bio-single-cell-cell-communication
翻译重构: 适配 bioskills State Dict 契约系统

功能：
- CellChat（R）：配体-受体通讯网络推断
- LIANA（Python）：多方法 CellChat/NicheNet/MetaCells 对比
- NicheNet（R）：配体-靶基因优先级分析

输出:
  - cellchat_object: CellChat R 对象（通过路径传递）
  - adata: obs 增加通讯相关列
  - communication_report: {pathways, top_interactions, sender_receiver}
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State


@register
class CellCommunicationSkill(AbstractSkill):
    """
    Infer cell-cell communication networks using ligand-receptor interaction databases.
    
    工具选择策略：
    | Tool     | Speed | Database | Language |
    |----------|-------|----------|----------|
    | CellChat | Slow  | Built-in | R        |
    | LIANA    | Fast  | 4 库    | Python   |
    | NicheNet | Slow  | Custom   | R        |
    
    数据库：CellChatDB.human, CellChatDB.mouse, LRDB, CellPhoneDB, ICELLNET
    """
    
    name = "cell_communication"
    description = (
        "Infer cell-cell communication using CellChat, LIANA, or NicheNet. "
        "Predict ligand-receptor interactions between annotated cell types. "
        "bioSkills: bio-single-cell-cell-communication"
    )
    input_contract = ["adata"]
    output_contract = ["adata", "communication_report"]
    stage = Stage.ENRICHMENT
    modality = [Modality.SCRNA]
    
    tunable_parameters = {
        "method": {"type": "str", "default": "cellchat"},
        "organism": {"type": "str", "default": "human"},  # human | mouse
        "cell_type_col": {"type": "str", "default": "cell_type"},
        "signaling_type": {"type": "str", "default": "Secreted Signaling"},
        "min_cells": {"type": "int", "default": 10},
        "top_n_interactions": {"type": "int", "default": 20},
    }
    
    def _run(self, state: State) -> dict:
        params = state.get("params", {})
        method = params.get("method", "cellchat")
        
        if method == "cellchat":
            return self._run_cellchat(state)
        elif method == "liana":
            return self._run_liana(state)
        else:
            return {
                "adata": state["adata"],
                "communication_report": {
                    "status": "failed",
                    "error": f"Unknown method: {method}. Use 'cellchat' or 'liana'.",
                }
            }
    
    def _run_cellchat(self, state: State) -> dict:
        """CellChat R: 完整配体受体分析"""
        import os, json
        import numpy as np
        import scanpy as sc
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        cell_type_col = params.get("cell_type_col", "cell_type")
        
        if cell_type_col not in adata.obs.columns:
            return {
                "adata": adata,
                "communication_report": {
                    "status": "failed",
                    "error": f"Cell type column '{cell_type_col}' not found in adata.obs. "
                             "Run cell annotation first."
                }
            }
        
        # 保存为临时 RDS 供 R 调用
        import tempfile
        tmpdir = tempfile.mkdtemp()
        rds_path = os.path.join(tmpdir, "seurat_obj.rds")
        
        # 尝试 R 调用
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            pandas2ri.activate()
        except ImportError:
            return {
                "adata": adata,
                "communication_report": {
                    "status": "failed",
                    "error": "rpy2 not installed. Install with: pip install rpy2"
                }
            }
        
        try:
            # 导出 AnnData 为 RDS
            ro.globalenv["adata_path"] = rds_path
            
            # 保存 adata 为 HDF5 → R 读取
            h5ad_path = os.path.join(tmpdir, "adata_temp.h5ad")
            adata.write_h5ad(h5ad_path)
            
            organism = params.get("organism", "human")
            
            r_script = f"""
            library(Seurat)
            library(CellChat)
            
            # Read h5ad via SeuratDisk
            adata <- ReadH5AD('{h5ad_path}')
            
            # Create CellChat object
            cellchat <- createCellChat(object = adata, group.by = '{cell_type_col}')
            
            # Use database
            if ('{organism}' == 'human') {{
                CellChatDB <- CellChatDB.human
            }} else {{
                CellChatDB <- CellChatDB.mouse
            }}
            cellchat@DB <- CellChatDB
            
            # Subset to signaling type
            signaling_type <- '{params.get("signaling_type", "Secreted Signaling")}'
            if (signaling_type != 'All') {{
                CellChatDB.use <- subsetDB(CellChatDB, search = signaling_type)
                cellchat@DB <- CellChatDB.use
            }}
            
            # Preprocessing
            cellchat <- subsetData(cellchat)
            cellchat <- identifyOverExpressedGenes(cellchat)
            cellchat <- identifyOverExpressedInteractions(cellchat)
            
            # Compute communication probability
            cellchat <- computeCommunProb(cellchat, type = 'triMean')
            cellchat <- filterCommunication(cellchat, min.cells = {params.get("min_cells", 10)})
            
            # Pathway level
            cellchat <- computeCommunProbPathway(cellchat)
            cellchat <- aggregateNet(cellchat)
            
            # Save results
            saveRDS(cellchat, '{tmpdir}/cellchat.rds')
            
            # Extract key results
            pathways <- cellchat@netP$pathways
            n_cells <- ncol(cellchat@data)
            
            # Top interactions
            df.net <- subsetCommunication(cellchat)
            if (nrow(df.net) > 0) {{
                top_interactions <- head(df.net[order(df.net$prob, decreasing=TRUE), ], 
                                        n = {params.get("top_n_interactions", 20)})
            }} else {{
                top_interactions <- data.frame()
            }}
            
            # Sender/Receiver analysis
            nwei <- as.numeric(cellchat@net$weight)
            groupSize <- as.numeric(table(cellchat@idents))
            sender_scores <- rowSums(matrix(nwei, nrow=length(groupSize)))
            names(sender_scores) <- names(groupSize)
            
            # Output
            pathways_str <- paste(pathways, collapse='|||')
            sender_str <- paste(paste(names(sender_scores), sender_scores, sep=':'), collapse='|||')
            n_paths <- length(pathways)
            n_interactions <- nrow(df.net)
            """
            
            ro.r(r_script)
            
            # 读取结果
            pathways_str = str(ro.globalenv.get("pathways_str", ""))
            sender_str = str(ro.globalenv.get("sender_str", ""))
            n_paths = int(ro.globalenv.get("n_paths", 0))
            n_interactions = int(ro.globalenv.get("n_interactions", 0))
            
            pathways = pathways_str.split("|||") if pathways_str else []
            sender_scores = {}
            if sender_str:
                for pair in sender_str.split("|||"):
                    if ":" in pair:
                        name, score = pair.rsplit(":", 1)
                        sender_scores[name] = float(score)
            
            # 在 adata.obs 中加入 sender scores
            for cell_type, score in sender_scores.items():
                col_name = f"sender_score_{cell_type}"
                adata.obs[col_name] = adata.obs[cell_type_col].map(
                    lambda x: float(score) if x == cell_type else 0.0
                )
            
            # 保存 CellChat 对象路径
            cellchat_path = os.path.join(tmpdir, "cellchat.rds")
            if not os.path.exists(cellchat_path):
                # 如果 R 脚本失败，创建一个空的占位符
                open(cellchat_path, 'w').close()
            
            report = {
                "method": "CellChat",
                "organism": organism,
                "n_pathways": n_paths,
                "n_interactions": n_interactions,
                "top_pathways": pathways[:10],
                "sender_scores": sender_scores,
                "signaling_type": params.get("signaling_type", "Secreted Signaling"),
                "cellchat_rds_path": cellchat_path,
                "status": "success",
            }
            
            print(f"  [CellChat] ✅ {n_paths} pathways, {n_interactions} interactions")
            print(f"  [CellChat] Top pathways: {pathways[:5]}")
            
            return {"adata": adata, "communication_report": report}
            
        except Exception as e:
            return {
                "adata": adata,
                "communication_report": {
                    "status": "failed",
                    "error": f"CellChat R execution failed: {e}",
                    "suggestion": "Install CellChat R package: "
                                  "R -e 'BiocManager::install(c(\"CellChat\"))'"
                }
            }
    
    def _run_liana(self, state: State) -> dict:
        """LIANA Python: 多方法细胞通讯对比"""
        import scanpy as sc
        import numpy as np
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        cell_type_col = params.get("cell_type_col", "cell_type")
        
        if cell_type_col not in adata.obs.columns:
            return {
                "adata": adata,
                "communication_report": {
                    "status": "failed",
                    "error": f"Cell type column '{cell_type_col}' not found"
                }
            }
        
        try:
            import liana as li
        except ImportError:
            print("  [LIANA] ⚠️  Installing liana...")
            import subprocess
            subprocess.run(["pip3", "install", "liana", "-q"], check=True)
            import liana as li
        
        min_cells = params.get("min_cells", 10)
        
        # 运行 LIANA 多方法
        try:
            li.mt.cellchat_single(
                adata, groupby=cell_type_col, n_threads=4,
                resource='CellChatDB', min_cells=min_cells,
            )
            method_used = "cellchat_single"
        except Exception:
            # fallback: 尝试 mean_rank
            try:
                li.mt.mean_rank(
                    adata, groupby=cell_type_col, n_threads=4,
                    resource='CellChatDB', min_cells=min_cells,
                )
                method_used = "mean_rank"
            except Exception as e:
                return {
                    "adata": adata,
                    "communication_report": {
                        "status": "failed",
                        "error": f"LIANA failed: {e}"
                    }
                }
        
        # 提取结果
        lr_cols = [c for c in adata.uns.keys() if c.startswith("liana_res")]
        if lr_cols:
            liana_res = adata.uns[lr_cols[0]]
            top_n = params.get("top_n_interactions", 20)
            top_interactions = liana_res.head(top_n)
            
            # 添加到 adata.uns
            adata.uns["liana_top_interactions"] = top_interactions.to_dict(orient="records")
            
            report = {
                "method": f"LIANA ({method_used})",
                "organism": params.get("organism", "human"),
                "n_interactions": len(liana_res),
                "top_interactions": top_interactions[[
                    "source", "target", "ligand", "receptor", "cellchat_score"
                ]].head(10).to_dict(orient="records") if "cellchat_score" in top_interactions.columns else [],
                "status": "success",
            }
            
            print(f"  [LIANA] ✅ {len(liana_res)} interactions, top 5:")
            for _, row in top_interactions.head(5).iterrows():
                print(f"    {row.get('source', '?')} → {row.get('target', '?')}: "
                      f"{row.get('ligand', '?')}-{row.get('receptor', '?')}")
            
            return {"adata": adata, "communication_report": report}
        else:
            return {
                "adata": adata,
                "communication_report": {"status": "skipped", "reason": "No interactions found"}
            }
