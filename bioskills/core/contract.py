"""L0: Contract Skill — 契约验证与技能规划"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
from typing import Tuple, List, Dict, Set, Optional, Any


class ContractError(Exception):
    def __init__(self, message: str, missing_keys=None):
        super().__init__(message)
        self.missing_keys = missing_keys or []


class ContractValidator:
    def __init__(self, registry=None):
        from bioskills.core.base import SkillRegistry
        self.registry = registry or SkillRegistry()

    def verify(self, skill_name: str, state: dict) -> dict:
        try:
            skill_cls = self.registry.get_class(skill_name)
        except KeyError:
            return {"valid": False, "error": f"Unknown skill: {skill_name}"}
        missing = [k for k in skill_cls.input_contract if k not in state or state[k] is None]
        return {"valid": len(missing) == 0, "missing_inputs": missing, "extra_outputs": []}

    def analyze_gaps(self, state: dict, target_contract: List[str]) -> dict:
        present = [k for k in target_contract if k in state and state[k] is not None]
        missing = [k for k in target_contract if k not in state or state[k] is None]
        return {"missing": missing, "present": present}

    def plan_pipeline(self, target_keys: List[str], state: dict) -> List[str]:
        plan = []
        covered = set(k for k, v in state.items() if v is not None)
        remaining = set(target_keys)
        max_iters = len(self.registry.list()) + 5
        for _ in range(max_iters):
            if not remaining:
                break
            made_progress = False
            for key in sorted(remaining):
                for sn in self.registry.find_producers(key):
                    if sn in plan:
                        continue
                    skill_cls = self.registry.get_class(sn)
                    if all(k in covered for k in skill_cls.input_contract):
                        plan.append(sn)
                        for out in skill_cls.output_contract:
                            covered.add(out)
                        remaining.discard(key)
                        made_progress = True
                        break
            if not made_progress:
                break
        return plan


@register
class ContractSkill(AbstractSkill):
    name = "contract"
    description = "Verify skill contracts, plan pipelines, analyze capability gaps"
    input_contract = ["params"]
    output_contract = ["verification_result", "pipeline_plan", "gap_analysis"]
    stage = Stage.DATA
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.SPATIAL, Modality.GENERIC]

    def _run(self, state: dict) -> dict:
        params = state.get("params", {})
        action = params.get("action", "verify")
        validator = ContractValidator(self.registry)
        if action == "verify":
            return {"verification_result": validator.verify(params.get("skill_name"), state),
                    "pipeline_plan": [], "gap_analysis": {}}
        elif action == "plan":
            return {"verification_result": {}, "pipeline_plan": validator.plan_pipeline(params.get("target_keys", []), state), "gap_analysis": {}}
        elif action == "analyze":
            return {"verification_result": {}, "pipeline_plan": [], "gap_analysis": validator.analyze_gaps(state, params.get("target_contract", []))}
        return {"verification_result": {}, "pipeline_plan": [], "gap_analysis": {}}
