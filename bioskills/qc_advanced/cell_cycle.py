"""
L0: CellCycle Skill — 细胞周期评分与回归

bioSkills 原版: 相关方法（bioSkills 内嵌于 preprocessing）
翻译重构: 独立 L0 Skill

功能：
- Cell Cycle Scoring（Scanpy）：S/G2M 基因集评分
- Regress Out Cell Cycle（Scanpy）：从 PCA 中移除细胞周期影响
- 适用于：去除细胞周期偏倚，使其他生物学信号更清晰

输出:
  - adata.obs: 增加 S_score, G2M_score, phase 列
  - cell_cycle_report: {S_score_mean, G2M_score_mean, dominant_phase}
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State


@register
class CellCycleSkill(AbstractSkill):
    """
    Score cell cycle phase and optionally regress out cell cycle effects.
    
    细胞周期基因集（bioSkills 标准）：
    - S genes: MCM5, PCLAF, ... (DNA replication)
    - G2M genes: HMGB2, CKAP2L, ... (G2/M transition)
    
    当细胞周期信号干扰生物学发现时（如比较肿瘤 vs 正常），
    使用 regress_out=True 从数据中移除细胞周期影响。
    """
    
    name = "cell_cycle"
    description = (
        "Score cell cycle phase (S/G2M) and regress out cell cycle effects. "
        "Use when cell cycle confounding masks biological variation. "
        "bioSkills pattern: cell cycle scoring in scRNA-seq preprocessing"
    )
    input_contract = ["adata"]
    output_contract = ["adata", "cell_cycle_report"]
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA]
    
    tunable_parameters = {
        "action": {"type": "str", "default": "score"},
        "species": {"type": "str", "default": "human"},  # human | mouse
        "regress_out": {"type": "bool", "default": False},
        "s_genes": {"type": "list", "default": []},  # 空=自动
        "g2m_genes": {"type": "list", "default": []},  # 空=自动
        "n_comps": {"type": "int", "default": 50},
    }
    
    def _run(self, state: State) -> dict:
        import scanpy as sc
        import numpy as np
        
        adata = state["adata"].copy()
        params = state.get("params", {})
        
        # 自动基因集（bioSkills 内置）
        # scanpy.queries.biology not available in all scanpy versions
        # Use built-in gene lists instead
        scq = None
        try:
            if params.get("species", "human") == "human":
                s_genes = params.get("s_genes") or scq.cell_cycle_genes("S")
                g2m_genes = params.get("g2m_genes") or scq.cell_cycle_genes("G2M")
            else:
                s_genes = params.get("s_genes") or []
                g2m_genes = params.get("g2m_genes") or []
        except Exception:
            # 备用标准基因集
            s_genes = params.get("s_genes") or [
                "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM3", "MCM4",
                "RFC2", "RPA2", "NASP", "RAD51", "GMNN", "CDC6", "CCND2",
            ]
            g2m_genes = params.get("g2m_genes") or [
                "HMGB2", "HMGB3", "HMGN1", "HMGN2", "HMGN3", "HMGN4",
                "TUBB", "TUBB4B", "CKAP2L", "CKAP2", "AURKB", "AURKA",
                "CDC20", "CENPF", "KIF2C", "KIF20B", "PLK1", "RANGAP1",
            ]
        
        # 过滤可用基因
        s_genes = [g for g in s_genes if g in adata.var_names]
        g2m_genes = [g for g in g2m_genes if g in adata.var_names]
        
        print(f"  [CellCycle] S genes: {len(s_genes)}, G2M genes: {len(g2m_genes)}")
        
        # Score
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
        
        s_mean = float(adata.obs["S_score"].mean())
        g2m_mean = float(adata.obs["G2M_score"].mean())
        
        phase_counts = adata.obs["phase"].value_counts()
        
        print(f"  [CellCycle] ✅ S_mean={s_mean:.3f}, G2M_mean={g2m_mean:.3f}")
        print(f"  [CellCycle] Phases: {dict(phase_counts)}")
        
        report = {
            "S_score_mean": round(s_mean, 4),
            "G2M_score_mean": round(g2m_mean, 4),
            "phase_distribution": phase_counts.to_dict(),
            "n_S_genes": len(s_genes),
            "n_G2M_genes": len(g2m_genes),
            "status": "success",
        }
        
        # Regress out（可选）
        if params.get("regress_out", False):
            print(f"  [CellCycle] Regressing out cell cycle from PCA...")
            sc.pp.scale(adata)
            sc.tl.pca(adata, n_comps=params.get("n_comps", 50))
            sc.pp.regress_out(adata, ["S_score", "G2M_score"])
            sc.pp.scale(adata)
            sc.tl.pca(adata, n_comps=params.get("n_comps", 50))
            report["regressed"] = True
            print(f"  [CellCycle] ✅ Cell cycle regressed from PCA")
        
        return {"adata": adata, "cell_cycle_report": report}
