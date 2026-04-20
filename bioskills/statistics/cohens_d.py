"""
L0: Cohen's d Skill — 效应量计算（核心统计）

从 adata.obs 读取 GSVA 评分，按 tumor/normal 分组计算 Cohen's d。
完全通用：任意基因集 × 任意细胞类型。
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu

@register
class CohensDSkill(AbstractSkill):
    """
    计算两组间的 Cohen's d 效应量 + Mann-Whitney U 检验。
    
    分组来源（优先级）：
    1. state['params']['groups'] = {"tumor": [barcodes], "normal": [barcodes]}
    2. adata.obs['group']（从 Phase 4 节点设置）
    3. adata.obs['sample_type'] / 'condition'
    4. barcode 前缀启发式（tumor/normal）
    
    输出：
      - effect_size (Cohen's d)
      - p_value (Mann-Whitney U, one-sided)
      - direction (increased/decreased)
      - confidence (high/medium/low)
      - groups_report
    """
    
    name = "cohens_d"
    description = (
        "Compute Cohen's d effect size and Mann-Whitney U test "
        "between two groups (tumor vs normal). "
        "Reads groups from adata.obs or params."
    )
    input_contract = ["adata"]
    output_contract = ["effect_size_result"]
    stage = Stage.STATISTICS
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.SPATIAL, Modality.GENERIC]
    
    tunable_parameters = {
        "score_column": {"type": "str", "default": ""},   # 空=自动
        "group_a": {"type": "str", "default": "tumor"},
        "group_b": {"type": "str", "default": "normal"},
        "alternative": {"type": "str", "default": "greater"},  # one-sided
        "confidence_threshold_high": {"type": "float", "default": 0.05},
        "effect_size_threshold_high": {"type": "float", "default": 0.5},
        "effect_size_threshold_medium": {"type": "float", "default": 0.2},
    }
    
    def _run(self, state: State) -> dict:
        import numpy as np
        from scipy.stats import mannwhitneyu
        
        adata = state["adata"]
        params = state.get("params", {})
        
        # ── Step 1: 找到 GSVA 评分列 ────────────
        score_col = params.get("score_column", "")
        
        if not score_col:
            # 自动找第一个 gsva_ 列
            gsva_cols = [c for c in adata.obs.columns if c.startswith("gsva_")]
            if gsva_cols:
                score_col = gsva_cols[0]
            else:
                return {
                    "effect_size_result": {
                        "status": "failed",
                        "error": "No GSVA score columns found in adata.obs. "
                                 "Run gsva skill first.",
                        "suggestion": "Add gsva skill to pipeline before cohens_d"
                    }
                }
        
        print(f"  [CohensD] Using score column: {score_col}")
        
        # ── Step 2: 确定分组 ────────────────────
        group_a_name = params.get("group_a", "tumor")
        group_b_name = params.get("group_b", "normal")
        
        groups = self._infer_groups(adata, params)
        
        if groups is None or not groups.get(group_a_name):
            return {
                "effect_size_result": {
                    "status": "failed",
                    "error": f"Cannot infer group '{group_a_name}' from data",
                    "suggestion": "Provide state['params']['groups'] with tumor/normal barcodes"
                }
            }
        
        a_barcodes = groups.get(group_a_name, [])
        b_barcodes = groups.get(group_b_name, [])
        
        # 取交集
        a_barcodes = [bc for bc in a_barcodes if bc in adata.obs_names]
        b_barcodes = [bc for bc in b_barcodes if bc in adata.obs_names]
        
        if not a_barcodes or not b_barcodes:
            return {
                "effect_size_result": {
                    "status": "failed",
                    "error": f"Empty group after intersection: "
                             f"{group_a_name}={len(a_barcodes)}, "
                             f"{group_b_name}={len(b_barcodes)}",
                }
            }
        
        scores_a = adata.obs.loc[a_barcodes, score_col].dropna().values.astype(float)
        scores_b = adata.obs.loc[b_barcodes, score_col].dropna().values.astype(float)
        
        if len(scores_a) < 3 or len(scores_b) < 3:
            return {
                "effect_size_result": {
                    "status": "failed",
                    "error": f"Too few samples: {group_a_name}={len(scores_a)}, "
                             f"{group_b_name}={len(scores_b)} (need ≥3 each)"
                }
            }
        
        # ── Step 3: Cohen's d ────────────────────
        mean_a = float(np.mean(scores_a))
        mean_b = float(np.mean(scores_b))
        std_a = float(np.std(scores_a, ddof=1))
        std_b = float(np.std(scores_b, ddof=1))
        
        n_a, n_b = len(scores_a), len(scores_b)
        pooled_sd = np.sqrt(((n_a - 1) * std_a**2 + (n_b - 1) * std_b**2) / (n_a + n_b - 2))
        
        if pooled_sd == 0:
            cohens_d = 0.0
        else:
            cohens_d = float((mean_a - mean_b) / pooled_sd)
        
        # ── Step 4: Mann-Whitney U ────────────────
        try:
            stat, p_two = mannwhitneyu(scores_a, scores_b, alternative='two-sided')
            p_value = float(min(p_two, 1.0))
        except Exception:
            stat, p_value = None, None
        
        # ── Step 5: 方向和置信度 ─────────────────
        direction = "increased" if cohens_d > 0 else "decreased"
        
        eff_h = params.get("effect_size_threshold_high", 0.5)
        eff_m = params.get("effect_size_threshold_medium", 0.2)
        p_thr = params.get("confidence_threshold_high", 0.05)
        
        if abs(cohens_d) >= eff_h and (p_value is not None and p_value < p_thr):
            confidence = "high"
        elif abs(cohens_d) >= eff_m:
            confidence = "medium"
        else:
            confidence = "low"
        
        result = {
            "status": "success",
            "score_column": score_col,
            "effect_size": round(cohens_d, 4),
            "p_value": round(p_value, 6) if p_value is not None else None,
            "direction": direction,
            "confidence": confidence,
            f"n_{group_a_name}": int(n_a),
            f"n_{group_b_name}": int(n_b),
            f"mean_{group_a_name}": round(mean_a, 4),
            f"mean_{group_b_name}": round(mean_b, 4),
            f"std_{group_a_name}": round(std_a, 4),
            f"std_{group_b_name}": round(std_b, 4),
            "groups_report": {
                f"{group_a_name}_barcodes": a_barcodes[:100],  # 截断避免太大
                f"{group_b_name}_barcodes": b_barcodes[:100],
                "total_a": len(a_barcodes),
                "total_b": len(b_barcodes),
            },
        }
        
        print(f"  [CohensD] ✅ d={cohens_d:.3f}, p={p_value:.4f if p_value else 'N/A'}, "
              f"direction={direction}, confidence={confidence}")
        
        return {"effect_size_result": result}
    
    def _infer_groups(self, adata, params) -> dict:
        """推断 tumor vs normal 分组"""
        # 1. 显式传入
        if "groups" in params:
            return params["groups"]
        
        # 2. adata.obs
        for col in ["group", "condition", "sample_type", "tumor_normal"]:
            if col in adata.obs.columns:
                groups = {}
                for cat in adata.obs[col].cat.categories:
                    barcodes = adata.obs_names[adata.obs[col] == cat].tolist()
                    groups[str(cat).lower()] = barcodes
                if groups:
                    return groups
        
        # 3. 启发式：从 barcode
        groups = {"tumor": [], "normal": []}
        keywords_tumor = ["tumor", "tumour", "cancer", "malignant", "primary"]
        keywords_normal = ["normal", "ctrl", "control", "nml", "healthy", "adjacent"]
        
        for bc in adata.obs_names:
            prefix = bc.rsplit("-", 1)[0].lower()
            if any(k in prefix for k in keywords_tumor):
                groups["tumor"].append(bc)
            elif any(k in prefix for k in keywords_normal):
                groups["normal"].append(bc)
        
        # 默认全部 tumor（保守策略）
        if not groups["tumor"] and not groups["normal"]:
            groups["tumor"] = adata.obs_names.tolist()
        
        return groups if (groups["tumor"] or groups["normal"]) else None
