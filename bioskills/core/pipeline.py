from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
from typing import List, Dict, Optional, Any
from dataclasses import dataclass, field
import time as _time


@dataclass
class PipelineResult:
    status: str = "pending"
    final_state: Optional[Dict[str, Any]] = None
    executed_skills: List[str] = field(default_factory=list)
    failed_skills: List[str] = field(default_factory=list)
    logs: List[str] = field(default_factory=list)
    artifacts: Dict[str, Any] = field(default_factory=dict)
    duration_seconds: float = 0.0
    error: Optional[str] = None
    gaps: List[Dict] = field(default_factory=list)
    target_keys: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def is_success(self) -> bool:
        return self.status == "success"

    def summary(self) -> str:
        n = len(self.executed_skills)
        f = len(self.failed_skills)
        return f"PipelineResult(status={self.status}, executed={n}, failed={f}, duration={self.duration_seconds:.2f}s)"

    def __repr__(self) -> str:
        return self.summary()


class PipelineOrchestrator:
    def __init__(self, registry=None, enable_self_heal=True, max_retries=3):
        from bioskills.core.base import SkillRegistry
        self.registry = registry or SkillRegistry()
        self.enable_self_heal = enable_self_heal
        self.max_retries = max_retries

    def execute(self, state: dict, target_keys: List[str], force_contract_aware=False):
        start = _time.time()
        logs = []
        executed = []
        failed = []
        current_state = dict(state)
        remaining = set(target_keys)
        max_iters = len(self.registry.list()) * 2
        for _ in range(max_iters):
            if not remaining:
                break
            made_progress = False
            for target in sorted(remaining):
                for sn in self.registry.find_producers(target):
                    if sn in plan if "plan" in dir() else []:
                        continue
                    skill_cls = self.registry.get_class(sn)
                    if all(k in current_state for k in skill_cls.input_contract):
                        logs.append(f"-> {sn} ({target})")
                        try:
                            result = skill_cls().execute(current_state)
                            if result.get("status") != "failed":
                                for k, v in result.items():
                                    if k not in ("logs",):
                                        current_state[k] = v
                                executed.append(sn)
                                remaining.discard(target)
                                made_progress = True
                            else:
                                err = result.get("error") or "unknown"
                                logs.append(f"  FAIL {sn}: {err}")
                                failed.append(sn)
                        except Exception as e:
                            logs.append(f"  FAIL {sn}: {type(e).__name__}: {e}")
                            failed.append(sn)
                        break
            if not made_progress:
                break
        return PipelineResult(
            status="success" if not remaining else "failed",
            final_state=current_state, executed_skills=executed,
            failed_skills=failed, logs=logs,
            duration_seconds=_time.time() - start, target_keys=target_keys,
        )


@register
class PipelineSkill(AbstractSkill):
    name = "pipeline"
    description = "Execute a predefined bioskills pipeline (cell_state_full, etc.)"
    input_contract = ["adata", "params"]
    output_contract = ["adata"]
    stage = Stage.DATA
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.SPATIAL, Modality.GENERIC]

    def _run(self, state: dict) -> dict:
        from bioskills.execution import execute_pipeline
        name = state.get("params", {}).get("pipeline", "cell_state_full")
        result = execute_pipeline(state, pipeline=name)
        return {
            "pipeline_result": result,
            "adata": result.final_state.get("adata") if result.final_state else None,
        }
