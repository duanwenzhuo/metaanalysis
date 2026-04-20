"""
L2 进化层 — 自愈与技能自动合成

L2SelfHealingPipeline:
    当 L1 无法规划出有效路径时：
    1. GapAnalyzer 识别能力缺口
    2. SkillSynthesizer 生成新技能代码（可选：LLM）
    3. 热加载并重试
"""

from bioskills.l2.self_heal import (
    L2SelfHealingPipeline,
    GapAnalyzer,
    SkillSynthesizer,
    CONTRACTS_DIR,
)

__all__ = [
    "L2SelfHealingPipeline",
    "GapAnalyzer",
    "SkillSynthesizer",
]
