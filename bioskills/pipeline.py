"""
L1: Pipeline Skill — 多技能编排执行器

bioSkills 原版: core/pipeline.py (STUB → 完全实现)
翻译重构: 按顺序执行多个技能，管理 State 传递，错误恢复

功能：
- run_pipeline(skills, state): 执行技能序列
- 支持 dry-run（只验证不执行）
- 每个技能的结果自动合并到 state
- 中途失败时报告失败点

输出:
  - final_state: 更新后的完整 state
  - execution_log: 每个技能的输入/输出/错误
  - pipeline_report: {n_skills, n_success, n_failed, elapsed_time}
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
from typing import List, Dict, Optional
import time, traceback


@register
class PipelineSkill(AbstractSkill):
    """
    Orchestrate multiple bioskills in sequence with state management.
    
    用法:
        pipeline = PipelineSkill()
        result = pipeline.run_pipeline(
            skills=["qc", "normalize", "hvg", "pca", "neighbors", "leiden"],
            state={"adata": adata, "params": {"min_genes": 200}}
        )
        
        final_state = result["final_state"]
        report = result["pipeline_report"]
    """
    
    name = "pipeline"
    description = (
        "Run a sequence of bioskills with automatic state propagation. "
        "Use for multi-step pipelines (QC → Normalize → HVG → PCA → Cluster → Annotate). "
        "bioSkills: pipeline orchestration"
    )
    input_contract = []  # 接受任意输入（从 state 读取 adata）
    output_contract = ["final_state", "execution_log", "pipeline_report"]
    stage = Stage.CORE
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.SPATIAL, Modality.GENERIC]
    
    tunable_parameters = {
        "skills": {"type": "list", "default": []},
        "dry_run": {"type": "bool", "default": False},
        "stop_on_error": {"type": "bool", "default": True},
        "report_every": {"type": "int", "default": 1},
    }
    
    def _run(self, state: State) -> dict:
        from bioskills.core.base import get_skill
        params = state.get("params", {})
        
        skills = params.get("skills", [])
        dry_run = params.get("dry_run", False)
        stop_on_error = params.get("stop_on_error", True)
        report_every = params.get("report_every", 1)
        
        if not skills:
            return {
                "final_state": state,
                "execution_log": [],
                "pipeline_report": {
                    "status": "skipped",
                    "reason": "No skills specified in params['skills']"
                }
            }
        
        print(f"  [Pipeline] 🚀 Starting pipeline: {' → '.join(skills)}")
        
        t0 = time.time()
        current_state = dict(state)  # 副本
        execution_log = []
        n_success = 0
        n_failed = 0
        
        for i, skill_name in enumerate(skills):
            idx = i + 1
            
            try:
                skill_class = get_skill(skill_name)
                skill_instance = skill_class()
            except Exception as e:
                log_entry = {
                    "skill": skill_name,
                    "index": idx,
                    "status": "not_found",
                    "error": str(e),
                    "elapsed_ms": 0,
                }
                execution_log.append(log_entry)
                n_failed += 1
                
                if stop_on_error:
                    break
                else:
                    continue
            
            # 执行技能
            t_skill = time.time()
            skill_params = current_state.get("params", {})
            
            if dry_run:
                status = "dry_run"
                skill_output = {}
            else:
                try:
                    skill_output = skill_instance.run(current_state)
                    status = "success"
                    n_success += 1
                except Exception as e:
                    status = "failed"
                    error_msg = f"{type(e).__name__}: {e}"
                    log_entry = {
                        "skill": skill_name,
                        "index": idx,
                        "status": status,
                        "error": error_msg,
                        "traceback": traceback.format_exc(),
                        "elapsed_ms": round((time.time() - t_skill) * 1000, 1),
                    }
                    execution_log.append(log_entry)
                    n_failed += 1
                    
                    if stop_on_error:
                        print(f"  [Pipeline] ❌ Skill {idx}/{len(skills)} '{skill_name}' failed: {error_msg}")
                        break
                    else:
                        continue
            
            elapsed_ms = round((time.time() - t_skill) * 1000, 1)
            
            log_entry = {
                "skill": skill_name,
                "index": idx,
                "status": status,
                "elapsed_ms": elapsed_ms,
                "output_keys": list(skill_output.keys()) if skill_output else [],
            }
            execution_log.append(log_entry)
            
            if idx % report_every == 0 or idx == len(skills):
                print(f"  [Pipeline]   {idx}/{len(skills)} ✅ {skill_name} ({elapsed_ms:.0f}ms)")
            
            # 合并输出到 state（跳过 adata 本身，避免复制开销）
            for key, value in skill_output.items():
                if key != "adata":
                    current_state[key] = value
            
            # 更新 adata（如果输出包含）
            if "adata" in skill_output:
                current_state["adata"] = skill_output["adata"]
        
        elapsed_total = round((time.time() - t0) * 1000, 1)
        
        report = {
            "status": "completed" if n_failed == 0 else "partial" if n_success > 0 else "failed",
            "n_skills_total": len(skills),
            "n_success": n_success,
            "n_failed": n_failed,
            "elapsed_ms": elapsed_total,
            "pipeline": skills,
            "failed_skills": [e["skill"] for e in execution_log if e["status"] == "failed"],
        }
        
        if report["status"] == "completed":
            print(f"  [Pipeline] ✅ Pipeline completed: {n_success}/{len(skills)} skills "
                  f"in {elapsed_total:.0f}ms")
        else:
            print(f"  [Pipeline] ⚠️  Pipeline partial: {n_success} ok, {n_failed} failed "
                  f"({report['failed_skills']})")
        
        return {
            "final_state": current_state,
            "execution_log": execution_log,
            "pipeline_report": report,
        }