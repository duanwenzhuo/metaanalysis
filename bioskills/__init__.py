"""
BioSkills — 三层自进化生物信息学技能系统

L0: 原子技能库（不可再分的最小功能单元）
L1: 复合层（动态管线和契约感知规划）
L2: 进化层（缺口分析 + 自动合成 + 自愈重试）

使用示例:
    from bioskills import SkillRegistry, AbstractSkill, register
    
    # 列出所有技能
    registry = SkillRegistry()
    print(registry.list())
    
    # 列出某类技能
    print(registry.list_by_stage(Stage.PREPROCESSING))
    
    # 获取技能实例
    qc = registry.get("qc")
    
    # 执行技能
    result = qc.execute({"adata": my_adata, "params": {"min_genes": 200}})
"""

from bioskills.core.base import (
    AbstractSkill,
    register,
    SkillRegistry,
    Stage,
    Modality,
    State,
    SkillInput,
    SkillOutput,
    StatusLiteral,
    ContractError,
    SkillExecutionError,
)
from bioskills.core.contract import ContractValidator

# 执行引擎
try:
    from bioskills.execution.pipeline_engine import (
        SkillPipelineEngine,
        execute_pipeline,
        pipeline_info,
    )
except ImportError:
    SkillPipelineEngine = None
    execute_pipeline = None
    pipeline_info = None

# Knowledge 层工具（可选导入）
try:
    from bioskills.knowledge.gene_set_db import (
        GeneSetDB,
        resolve_gene_set,
        CELL_MARKER_DB,
        BIOLOGICAL_PROCESS_DB,
    )
except ImportError:
    GeneSetDB = None
    resolve_gene_set = lambda *a, **k: {}
    CELL_MARKER_DB = {}
    BIOLOGICAL_PROCESS_DB = {}

__version__ = "1.0.0"
__all__ = [
    # 核心类
    "AbstractSkill",
    "register",
    "SkillRegistry",
    "Stage",
    "Modality",
    "State",
    "SkillInput",
    "SkillOutput",
    "StatusLiteral",
    "ContractError",
    "SkillExecutionError",
    "ContractValidator",
    # 执行引擎
    "SkillPipelineEngine",
    "execute_pipeline",
    "pipeline_info",
    # Knowledge 工具
    "GeneSetDB",
    "resolve_gene_set",
    "CELL_MARKER_DB",
    "BIOLOGICAL_PROCESS_DB",
]
