"""core/types.py — 核心数据类型定义（Python 3.9 兼容）。"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import TypedDict, Any, Optional, Dict, List, Literal


# ─────────────────────────────────────────────
# 核心数据类型（dataclass）
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


# ─────────────────────────────────────────────
# PhaseStatus 字面量（Python 3.9 兼容）
# ─────────────────────────────────────────────

PhaseStatusLiteral = Literal["pending", "running", "done", "failed"]


# ─────────────────────────────────────────────
# BioAnalysisState — 贯穿整个流程的共享状态
# ─────────────────────────────────────────────

class BioAnalysisState(TypedDict, total=False):
    """
    LangGraph StateSchema。
    贯穿 HypothesisParser-6 的所有共享状态。
    """

    # ── 元信息 ──
    hypothesis_id: str
    hypothesis_text: str
    created_at: str

    # ── HypothesisParser ──
    hypothesis_parser_type: str
    hypothesis_parser_entity1: str
    hypothesis_parser_entity2: str
    hypothesis_parser_relation: str
    hypothesis_parser_strategy: str
    hypothesis_parser_target_genes: List[str]
    hypothesis_parser_status: PhaseStatusLiteral
    hypothesis_parser_error: Optional[str]

    # ── LiteratureMiner ──
    literature_miner_supporting_papers: List[dict]
    literature_miner_gene_sets: Dict[str, List[str]]
    literature_miner_recommended_datasets: List[dict]
    literature_miner_status: PhaseStatusLiteral
    literature_miner_error: Optional[str]

    # ── GEORetriever ──
    geo_retriever_datasets: List[dict]
    geo_retriever_n_downloaded: int
    geo_retriever_n_failed: int
    geo_retriever_failed_datasets: List[dict]
    geo_retriever_status: PhaseStatusLiteral
    geo_retriever_error: Optional[str]

    # ── EffectSizeAnalyst ──
    effect_size_analyst_results: Dict[str, dict]
    effect_size_analyst_n_analyzed: int
    effect_size_analyst_n_total: int
    effect_size_analyst_errors: List[dict]
    effect_size_analyst_status: PhaseStatusLiteral
    effect_size_analyst_error: Optional[str]

    # ── MetaAnalyst ──
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

    # ── ReportCurator ──
    report_curator_report_path: Optional[str]
    report_curator_status: PhaseStatusLiteral
    report_curator_error: Optional[str]

    # ── 全局字段 ──
    current_phase: str
    errors: List[str]
    checkpoints: List[str]
    updated_at: Optional[str]
