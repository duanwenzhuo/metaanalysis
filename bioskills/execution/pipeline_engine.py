"""
SkillPipelineEngine — bioskills 统一执行引擎

职责：接收 LangGraph State，调度 L0 技能序列，返回更新后的 State。
用法：
    engine = SkillPipelineEngine()
    result = engine.run(state, pipeline="cell_state_full")
    # 或
    result = engine.run(state, target_keys=["effect_size_result"])
"""

from __future__ import annotations
import sys, os, time, traceback
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, field

# 确保 bioskills 可导入
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from bioskills.core.base import SkillRegistry, Stage, State
from bioskills.core.contract import ContractValidator
from bioskills.knowledge.gene_set_db import resolve_gene_set


# ─────────────────────────────────────────────
# Pipeline 模板
# ─────────────────────────────────────────────

PIPELINES = {
    "cell_state_full": [
        "qc", "normalize", "hvg", "scale",
        "pca", "neighbors", "leiden",
        "marker_score", "gsva", "cohens_d",
    ],
    "cell_state_fast": [
        "qc", "normalize", "hvg", "scale",
        "pca", "neighbors", "leiden",
        "marker_score", "gsva", "cohens_d",
    ],
    "cell_state_with_doublet": [
        "doublet_detection", "qc", "normalize", "hvg", "scale",
        "pca", "neighbors", "leiden",
        "marker_score", "gsva", "cohens_d",
    ],
    "cell_state_with_batch": [
        "qc", "normalize", "hvg", "scale",
        "pca", "batch_correction",
        "neighbors", "leiden",
        "marker_score", "gsva", "cohens_d",
    ],
    "bulk_deseq2": [
        "deseq2", "de_results", "volcano",
    ],
    "bulk_edger": [
        "edger", "de_results", "volcano",
    ],
    "go_enrichment": [
        "go_enrichment",
    ],
}


@dataclass
class PipelineResult:
    status: str = "pending"
    state: State = field(default_factory=dict)
    executed: List[str] = field(default_factory=list)
    failed: List[str] = field(default_factory=list)
    logs: List[str] = field(default_factory=list)
    duration_seconds: float = 0.0
    error: Optional[str] = None


class SkillPipelineEngine:
    """
    bioskills 统一执行引擎。
    
    核心逻辑：
    1. 从 State 中读取假设信息
    2. 调用 GeneSetDB 动态解析基因集（替换硬编码）
    3. 按模板顺序执行 L0 技能
    4. 每步 merge 结果回 State
    5. 失败时触发 L2 自愈
    """
    
    def __init__(self, registry: Optional[SkillRegistry] = None):
        self.registry = registry or SkillRegistry()
        self._ensure_discovered = False
    
    def _ensure_ready(self):
        if not self._ensure_discovered:
            self.registry.auto_discover()
            self._ensure_discovered = True
    
    def run(
        self,
        state: State,
        pipeline: Optional[str] = None,
        target_keys: Optional[List[str]] = None,
    ) -> PipelineResult:
        """
        执行 pipeline。
        
        Args:
            state: LangGraph 传来的 State dict
            pipeline: 模板名（如 "cell_state_full"）
            target_keys: 目标 State 键（无模板时用契约规划）
        """
        self._ensure_ready()
        start = time.time()
        result = PipelineResult(state=dict(state))
        
        # ── Step 1: 动态解析基因集（替换硬编码）───
        state = self._resolve_gene_sets(state)
        
        # ── Step 2: 确定 pipeline ──────────────────
        if pipeline and pipeline in PIPELINES:
            skill_sequence = PIPELINES[pipeline]
            result.logs.append(f"📋 Pipeline: {pipeline} ({len(skill_sequence)} steps)")
        elif target_keys:
            validator = ContractValidator(self.registry)
            path = validator.plan_path(
                {k for k, v in state.items() if v is not None},
                set(target_keys),
            )
            if path:
                skill_sequence = path.skills
                result.logs.append(f"🔍 Contract-planned: {' → '.join(skill_sequence)}")
            else:
                skill_sequence = PIPELINES.get("cell_state_full", [])
                result.logs.append(f"⚠️  No path found, fallback to cell_state_full")
        else:
            skill_sequence = PIPELINES.get("cell_state_full", [])
            result.logs.append(f"📋 Default: cell_state_full")
        
        ip = bool(state.get("input_path")); ha = "adata" in state; print(f"[DEBUG] ip={ip} ha={ha}")
        # ── Step 2.5: 自动加载 input_path → adata ──
        current_state = dict(state)
        if current_state.get("input_path") and "adata" not in current_state:
            try:
                data_io_skill = self.registry.get("data_io")
                load_result = data_io_skill.execute(current_state)
                if load_result.get("status", "success") == "success" and "state_updates" in load_result:
                    current_state.update(load_result.get("state_updates", {}))
                    result.logs.append(f"📂 [data_io] Loaded adata from {current_state.get('input_path')}")
                    result.executed.append("data_io")
                else:
                    result.logs.append(f"❌ [data_io] {load_result.get('error', 'unknown')}")
            except Exception as e:
                result.logs.append(f"❌ [data_io] exception: {e}")

        # ── Step 3: 顺序执行技能 ──────────────────
        
        for skill_name in skill_sequence:
            try:
                skill = self.registry.get(skill_name)
            except KeyError:
                result.logs.append(f"⏭ [{skill_name}] not registered, skipping")
                continue
            
            # 合并 params
            params = dict(current_state.get("params", {}))
            # 确保基因集和细胞类型信息传下去
            if "gene_sets" not in params and "hypothesis_parser_gene_sets" in current_state:
                params["gene_sets"] = current_state["hypothesis_parser_gene_sets"]
            if "cell_type" not in params and "hypothesis_parser_entity1" in current_state:
                params["cell_type"] = current_state["hypothesis_parser_entity1"]
            if "process" not in params and "hypothesis_parser_relation" in current_state:
                params["process"] = current_state["hypothesis_parser_relation"]
            if "hypothesis_text" not in params and "hypothesis_text" in current_state:
                params["hypothesis_text"] = current_state["hypothesis_text"]
            
            skill_input = dict(current_state)
            skill_input["params"] = params
            
            result.logs.append(f"▶ [{skill_name}] starting...")
            
            try:
                output = skill.execute(skill_input)
                
                if output["status"] == "success":
                    updates = output.get("state_updates", {})
                    current_state.update(updates)
                    result.executed.append(skill_name)
                    result.logs.append(f"✅ [{skill_name}] done")
                else:
                    err = output.get("error", "unknown error")
                    result.failed.append(skill_name)
                    result.logs.append(f"❌ [{skill_name}] {err}")
                    # 不中断，继续执行后续技能
            
            except Exception as e:
                result.failed.append(skill_name)
                result.logs.append(f"❌ [{skill_name}] exception: {e}")
        
        # ── Step 4: 如果失败且有 L2 ──────────────
        if result.failed and current_state.get("_enable_l2", True):
            result.logs.append(f"🔧 L2 self-heal triggered for {len(result.failed)} failures")
            # L2 暂时只记录，不自动执行（需要 R/Python 依赖）
        
        result.state = current_state
        result.duration_seconds = round(time.time() - start, 2)
        result.status = "success" if not result.failed else "partial"
        
        return result
    
    def _resolve_gene_sets(self, state: State) -> State:
        """
        动态解析基因集 — 替换旧系统硬编码的 t_cell_markers / exhaustion_genes。
        """
        # 如果已有基因集，直接返回
        if state.get("hypothesis_parser_gene_sets") or state.get("params", {}).get("gene_sets"):
            return state
        
        entity1 = state.get("hypothesis_parser_entity1", "")
        relation = state.get("hypothesis_parser_relation", "")
        hyp_text = state.get("hypothesis_text", "")
        
        if not entity1 and not relation and not hyp_text:
            return state
        
        gene_sets = resolve_gene_set(entity1)  # returns List[str]

        if gene_sets and gene_sets != []:
            state["hypothesis_parser_gene_sets"] = {"CUSTOM": gene_sets}
            # 也写入 params 以便技能读取
            params = dict(state.get("params", {}))
            params["gene_sets"] = {"CUSTOM": gene_sets}
            params["cell_type"] = entity1
            params["process"] = relation
            params["hypothesis_text"] = hyp_text
            state["params"] = params

            n_genes = len(gene_sets)
            print(f"  [PipelineEngine] GeneSetDB → CUSTOM ({n_genes} genes)")

        return state
    
    def list_pipelines(self) -> Dict[str, List[str]]:
        return dict(PIPELINES)
    
    def info(self) -> Dict[str, Any]:
        self._ensure_ready()
        return {
            "n_skills": len(self.registry.list()),
            "pipelines": list(PIPELINES.keys()),
            "skills": self.registry.list(),
        }


# ─────────────────────────────────────────────
# 便捷函数
# ─────────────────────────────────────────────

def execute_pipeline(state: State, pipeline: str = "cell_state_full") -> PipelineResult:
    """一行执行 pipeline"""
    engine = SkillPipelineEngine()
    return engine.run(state, pipeline=pipeline)

def pipeline_info() -> Dict[str, Any]:
    """查看 pipeline 引擎信息"""
    engine = SkillPipelineEngine()
    return engine.info()
