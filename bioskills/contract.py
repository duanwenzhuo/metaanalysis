"""
L0: Contract Skill — 契约验证与技能规划

bioSkills 原版: core/contract.py (STUB → 完全实现)
翻译重构: 验证 State Dict 契约，检查输入/输出，生成 gap 报告

功能：
- verify_skill(name, state): 验证技能能否在当前 state 上运行
- plan_pipeline(goal, state): 从目标反推需要的技能序列
- analyze_gaps(state, target_contract): 识别当前 state 缺少哪些字段

输出:
  - verification_result: {valid, missing_inputs, extra_outputs}
  - pipeline_plan: list of skill names in execution order
  - gap_analysis: {missing: [], present: []}
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
from typing import Tuple, List, Dict
import json


@register
class ContractSkill(AbstractSkill):
    """
    Verify bioskills contracts and plan skill pipelines.
    
    核心功能：
    1. **verify_skill**: 验证 skill 的输入契约是否满足
    2. **plan_pipeline**: 从目标反向规划技能序列
    3. **analyze_gaps**: 识别当前 state 缺少哪些字段
    
    用法示例：
        skill = ContractSkill()
        valid, errors = skill.verify("marker_score", state)
        plan = skill.plan_pipeline("Annotate cell types → DE → Enrich", state)
    """
    
    name = "contract"
    description = (
        "Verify bioskills input/output contracts and plan skill pipelines. "
        "Pre-flight check before running any pipeline."
    )
    input_contract = []  # 可接受任意 state
    output_contract = ["verification_result", "pipeline_plan", "gap_analysis"]
    stage = Stage.CORE
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.SPATIAL, Modality.GENERIC]
    
    tunable_parameters = {
        "skill_to_verify": {"type": "str", "default": ""},
        "goal_text": {"type": "str", "default": ""},
        "target_inputs": {"type": "list", "default": []},
        "strict": {"type": "bool", "default": False},
    }
    
    def _run(self, state: State) -> dict:
        params = state.get("params", {})
        skill_to_verify = params.get("skill_to_verify", "")
        goal_text = params.get("goal_text", "")
        target_inputs = params.get("target_inputs", [])
        strict = params.get("strict", False)
        
        results = {}
        
        # ── 验证指定技能 ─────────────────────────
        if skill_to_verify:
            results["verification_result"] = self._verify_skill(skill_to_verify, state, strict)
        
        # ── 规划管线 ────────────────────────────
        if goal_text:
            results["pipeline_plan"] = self._plan_pipeline(goal_text, state)
        
        # ── Gap 分析 ───────────────────────────
        if target_inputs:
            results["gap_analysis"] = self._analyze_gaps(state, target_inputs)
        
        # 默认：返回 gap 分析（如果没指定其他）
        if not results:
            results["gap_analysis"] = self._analyze_gaps(state, [])
        
        return results
    
    def _verify_skill(self, skill_name: str, state: State, strict: bool) -> dict:
        """验证技能输入契约"""
        # 获取技能定义
        from bioskills.core.base import get_skill
        
        try:
            skill_class = get_skill(skill_name)
        except Exception:
            return {
                "valid": False,
                "error": f"Skill '{skill_name}' not found",
                "missing_inputs": [],
            }
        
        required_inputs = getattr(skill_class, "input_contract", [])
        present_keys = set(state.keys())
        missing = [inp for inp in required_inputs if inp not in present_keys]
        extra = [k for k in present_keys if k not in required_inputs and k not in ["params", "adata"]]
        
        valid = len(missing) == 0
        
        result = {
            "skill_name": skill_name,
            "valid": valid,
            "required_inputs": required_inputs,
            "missing_inputs": missing,
            "extra_keys": extra if strict else [],
            "state_keys": list(present_keys),
        }
        
        if valid:
            print(f"  [Contract] ✅ '{skill_name}' verified (all inputs present)")
        else:
            print(f"  [Contract] ⚠️  '{skill_name}' missing: {missing}")
        
        return result
    
    def _plan_pipeline(self, goal_text: str, state: State) -> dict:
        """从目标文本反向规划技能序列"""
        # 基于关键词的简单规划器（可扩展为 LLM）
        goal = goal_text.lower()
        
        # 定义技能-关键词映射
        keyword_map = {
            "normalize": ["normalize", "library size", "log-transform", "scaling"],
            "hvg": ["hvg", "highly variable", "feature selection", "HVG"],
            "pca": ["pca", "principal component", "dimensionality reduction"],
            "neighbors": ["neighbors", "neighbor graph", "k-neighbor"],
            "leiden": ["cluster", "leiden", "community detection"],
            "marker_score": ["annotate", "cell type", "marker", "cell annotation"],
            "gsva": ["gsva", "enrichment score", "pathway activity"],
            "deseq2": ["differential expression", "de", "deseq", "edgeR"],
            "go_enrichment": ["enrichment", "go", "kegg", "pathway"],
            "volcano": ["volcano", "visualize", "logfc", "p-value"],
            "trajectory_inference": ["trajectory", "pseudotime", "developmental"],
            "cell_chat": ["cell chat", "cellphone", "communication", "ligand-receptor"],
            "batch_correction": ["batch", "integrate", "harmony", "scvi"],
        }
        
        # 启发式排序（拓扑顺序）
        skill_order = [
            "data_io", "qc", "normalize", "scale",
            "cell_cycle", "doublet_detection",
            "hvg", "pca", "neighbors", "batch_correction",
            "leiden", "louvain",
            "marker_score", "cell_annotation_ref",
            "gsva", "deseq2", "edger",
            "cohens_d", "mann_whitney",
            "go_enrichment", "gsea",
            "volcano", "de_visualization",
            "trajectory_inference",
            "cell_chat",
        ]
        
        selected = []
        mentioned = set()
        
        # 匹配关键词
        for skill, keywords in keyword_map.items():
            if any(kw in goal for kw in keywords):
                if skill not in mentioned:
                    selected.append(skill)
                    mentioned.add(skill)
        
        # 推断必要的前置技能
        inferred = self._infer_prerequisites(selected, state)
        
        plan = inferred + selected
        
        print(f"  [Contract] 📋 Pipeline planned: {' → '.join(plan)}")
        
        return {
            "plan": plan,
            "n_skills": len(plan),
            "goal": goal_text,
            "reasoning": "keyword matching + topological sort",
        }
    
    def _infer_prerequisites(self, selected: List[str], state: State) -> List[str]:
        """推断前置技能"""
        prereqs = []
        has_adata = "adata" in state
        
        # 无 adata 时必须先有 data_io
        if not has_adata:
            return ["data_io"]
        
        adata = state.get("adata")
        has_pca = adata is not None and "X_pca" in adata.obsm if adata else False
        has_neighbors = adata is not None and "neighbors" in adata.uns if adata else False
        has_clusters = adata is not None and ("leiden" in adata.obs or "louvain" in adata.obs or "clusters" in adata.obs) if adata else False
        
        for skill in selected:
            if skill in ["pca", "neighbors", "leiden", "louvain"] and not has_pca:
                if "hvg" not in prereqs:
                    prereqs.append("hvg")
            if skill in ["leiden", "louvain"] and not has_neighbors:
                if "neighbors" not in prereqs:
                    prereqs.append("neighbors")
            if skill in ["gsva", "marker_score"]:
                if "normalize" not in prereqs and "scale" not in prereqs:
                    prereqs.append("normalize")
        
        return prereqs
    
    def _analyze_gaps(self, state: State, target_inputs: List[str]) -> dict:
        """分析缺失字段"""
        present = list(state.keys())
        
        if not target_inputs:
            # 默认：检查 adata 相关字段
            target_inputs = ["adata"]
        
        missing = [inp for inp in target_inputs if inp not in present]
        has = [inp for inp in target_inputs if inp in present]
        
        # 特殊检查 adata 状态
        adata_status = "not_in_state"
        if "adata" in state:
            adata = state["adata"]
            adata_status = f"{adata.n_obs} cells × {adata.n_vars} vars"
        
        return {
            "present": has,
            "missing": missing,
            "adata_status": adata_status,
            "all_keys": present,
        }