"""
Bio-Analysis State Schema (TypedDict)

整个多智能体系统的共享状态定义。
Python 3.9 兼容版（无 | 语法，无 NotRequired）。

设计原则：
- 所有 Phase 共用的字段放这里
- 每个 Agent 只修改自己负责的字段
- 不允许写入未知字段
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import (
    TypedDict, Any, Optional, Dict, List, Set, Union,
    Literal, Callable, Sequence,
)


# ─────────────────────────────────────────────
# 核心数据类型
# ─────────────────────────────────────────────

@dataclass
class DatasetInfo:
    """单个 GEO 数据集的信息"""
    gse_id: str
    disease: str
    cells: int
    platform: str
    n_tumor: int = 0
    n_normal: int = 0
    tier: str = "P1"           # P0=必须, P1=推荐, P2=可选
    priority_score: float = 0.5
    status: str = "pending"     # pending / downloaded / analyzed / failed
    error: Optional[str] = None
    file_path: Optional[str] = None


@dataclass
class LiteratureSource:
    """文献来源"""
    pmid: str
    title: str
    journal: str
    year: int
    citation_count: int = 0
    geo_ids: List[str] = field(default_factory=list)
    key_findings: List[str] = field(default_factory=list)


# ─────────────────────────────────────────────
# PhaseStatus 字面量（Python 3.9 兼容）
# ─────────────────────────────────────────────

PhaseStatusLiteral = Literal["pending", "running", "done", "failed", "no_data"]


# ─────────────────────────────────────────────
# BioAnalysisState — 贯穿整个流程的共享状态
# ─────────────────────────────────────────────

@dataclass
class DatasetResult:
    """单个数据集的 Meta 分析结果"""
    gse_id: str
    effect_size: float          # Cohen's d
    p_value: float            # 统计显著性
    direction: int             # 1=正效应, -1=负效应, 0=中性
    n_samples: int
    disease: str
    platform: str
    genes_used: List[str]
    n_high: int = 0
    n_low: int = 0
    method: str = "GSVA"
    status: str = "success"
    error: Optional[str] = None


class BioAnalysisState(TypedDict, total=False):
    """
    LangGraph StateSchema。
    贯穿 Phase 1-6 的所有共享状态。
    
    字段命名规则：
    - hypothesis_parser_* = Phase 1 输出
    - literature_miner_* = Phase 2 输出
    - geo_retriever_* = Phase 3 输出
    - effect_size_analyst_* = Phase 4 输出
    - meta_analyst_* = Phase 5 输出
    - report_curator_* = Phase 6 输出
    """

    # ── 元信息 ──
    hypothesis_id: str
    hypothesis_text: str
    created_at: str

    # ── Phase 1 ──
    hypothesis_parser_type: str
    hypothesis_parser_entity1: str
    hypothesis_parser_entity2: str
    hypothesis_parser_relation: str
    hypothesis_parser_strategy: str
    hypothesis_parser_target_genes: List[str]
    hypothesis_parser_status: PhaseStatusLiteral
    hypothesis_parser_error: Optional[str]

    # ── Phase 2 ──
    literature_miner_supporting_papers: List[dict]
    literature_miner_gene_sets: Dict[str, List[str]]
    literature_miner_recommended_datasets: List[dict]
    literature_miner_status: PhaseStatusLiteral
    literature_miner_error: Optional[str]

    # ── Phase 3 ──
    geo_retriever_datasets: List[dict]
    geo_retriever_n_downloaded: int
    geo_retriever_n_failed: int
    geo_retriever_failed_datasets: List[dict]
    geo_retriever_status: PhaseStatusLiteral
    geo_retriever_error: Optional[str]

    # ── Phase 4 ──
    effect_size_analyst_results: Dict[str, dict]       # {gse_id: DatasetResult}
    effect_size_analyst_n_analyzed: int
    effect_size_analyst_n_total: int
    effect_size_analyst_errors: List[dict]
    effect_size_analyst_status: PhaseStatusLiteral
    effect_size_analyst_error: Optional[str]

    # ── Phase 5 ──
    meta_analyst_reliability_score: float
    meta_analyst_support: float
    meta_analyst_consistency: float
    meta_analyst_significance: float
    meta_analyst_heterogeneity: float
    meta_analyst_combined_d: float
    meta_analyst_ci_lower: float
    meta_analyst_ci_upper: float
    meta_analyst_combined_p: float
    meta_analyst_i_squared: float
    meta_analyst_loo_robust: str
    meta_analyst_grade: str
    meta_analyst_eggertest_p: float
    meta_analyst_status: PhaseStatusLiteral
    meta_analyst_error: Optional[str]

    # ── Phase 6 ──
    report_curator_report_path: Optional[str]
    report_curator_status: PhaseStatusLiteral
    report_curator_error: Optional[str]

    # ── 全局字段 ──
    current_phase: str
    errors: List[str]
    checkpoints: List[str]
    updated_at: Optional[str]


# ─────────────────────────────────────────────
# 便捷构造函数
# ─────────────────────────────────────────────

def new_state(hypothesis_text: str) -> BioAnalysisState:
    """创建初始状态（Python 3.9 兼容）"""
    import uuid
    from datetime import datetime, timezone

    hyp_id = "h_{}_{}".format(
        datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S"),
        uuid.uuid4().hex[:6],
    )

    state: BioAnalysisState = {
        "hypothesis_id": hyp_id,
        "hypothesis_text": hypothesis_text,
        "created_at": datetime.now(timezone.utc).isoformat(),

        # 默认 Phase 状态
        "hypothesis_parser_status": "pending",
        "literature_miner_status": "pending",
        "geo_retriever_status": "pending",
        "effect_size_analyst_status": "pending",
        "meta_analyst_status": "pending",
        "report_curator_status": "pending",

        # 全局
        "current_phase": "phase1",
        "errors": [],
        "checkpoints": [],

        # Phase 1 defaults
        "hypothesis_parser_type": "",
        "hypothesis_parser_entity1": "",
        "hypothesis_parser_entity2": "",
        "hypothesis_parser_relation": "",
        "hypothesis_parser_strategy": "literature_first",
        "hypothesis_parser_target_genes": [],

        # Phase 2 defaults
        "literature_miner_gene_sets": {},
        "literature_miner_supporting_papers": [],
        "literature_miner_recommended_datasets": [],

        # Phase 3 defaults
        "geo_retriever_datasets": [],
        "geo_retriever_n_downloaded": 0,
        "geo_retriever_n_failed": 0,
        "geo_retriever_failed_datasets": [],

        # Phase 4 results
        "effect_size_analyst_results": {},
        "effect_size_analyst_n_analyzed": 0,
        "effect_size_analyst_n_total": 0,
        "effect_size_analyst_errors": [],

        # Phase 5 defaults
        "meta_analyst_reliability_score": 0.0,
        "meta_analyst_support": 0.0,
        "meta_analyst_consistency": 0.0,
        "meta_analyst_significance": 0.0,
        "meta_analyst_heterogeneity": 0.0,
        "meta_analyst_combined_d": 0.0,
        "meta_analyst_ci_lower": 0.0,
        "meta_analyst_ci_upper": 0.0,
        "meta_analyst_combined_p": 1.0,
        "meta_analyst_i_squared": 0.0,
        "meta_analyst_loo_robust": "unknown",
        "meta_analyst_grade": "INSUFFICIENT",
        "meta_analyst_eggertest_p": 1.0,

        # Phase 6 defaults
        "report_curator_report_path": "",
        "effect_size_analyst_n_total": 0,
        "effect_size_analyst_errors": [],

        # Phase 3
        "geo_retriever_datasets": [],
        "geo_retriever_n_downloaded": 0,
        "geo_retriever_n_failed": 0,
        "geo_retriever_failed_datasets": [],

        # Phase 5 默认
        "meta_analyst_reliability_score": 0.0,
        "meta_analyst_support": 0.0,
        "meta_analyst_consistency": 0.0,
        "meta_analyst_significance": 0.0,
        "meta_analyst_heterogeneity": 0.0,
        "meta_analyst_combined_d": 0.0,
        "meta_analyst_ci_lower": 0.0,
        "meta_analyst_ci_upper": 0.0,
        "meta_analyst_combined_p": 1.0,
        "meta_analyst_i_squared": 0.0,
        "meta_analyst_eggertest_p": 0.0,
        "meta_analyst_grade": "⚪ N/A",
        "meta_analyst_loo_robust": "N/A",
    }

    return state


def validate_state(state: BioAnalysisState) -> List[str]:
    """验证状态合法性，返回错误列表（空=合法）"""
    errors = []

    required = ["hypothesis_id", "hypothesis_text", "current_phase"]
    for field_name in required:
        if field_name not in state or not state.get(field_name):
            errors.append("Missing required field: {}".format(field_name))

    # Agent 状态值校验
    valid_statuses = {"pending", "running", "done", "failed", "no_data"}
    for phase in ["phase1", "phase2", "phase3", "phase4", "phase5", "phase6"]:
        status_field = "{}_status".format(phase)
        if status_field in state:
            val = state.get(status_field)
            if val not in valid_statuses:
                errors.append("Invalid {}: {}".format(status_field, val))

    return errors


# ─────────────────────────────────────────────
# Phase 状态转移规则
# ─────────────────────────────────────────────

VALID_TRANSITIONS = {
    "phase1": ["phase2"],
    "phase2": ["phase3"],
    "phase3": ["phase4"],
    "phase4": ["phase5"],
    "phase5": ["phase6"],
    "phase6": [],  # 终态
}


def next_phase(phase: str) -> Optional[str]:
    """获取下一个 Phase（若无则返回 None）"""
    phases = list(VALID_TRANSITIONS.keys())
    try:
        idx = phases.index(phase)
        if idx + 1 < len(phases):
            return phases[idx + 1]
    except ValueError:
        pass
    return None


def phase_status_label(status: str) -> str:
    """Phase 状态的中文标签"""
    labels = {
        "pending": "⏳ 待执行",
        "running": "🔄 执行中",
        "done": "✅ 完成",
        "failed": "❌ 失败",
    }
    return labels.get(status, status)
