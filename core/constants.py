"""core/constants.py — 全局常量定义。"""

import os
from pathlib import Path

# ── 工作空间路径 ──
WORKSPACE_DIR = Path(os.environ.get(
    "BIO_ANALYSIS_WORKSPACE",
    Path(__file__).parent.parent
))

SCRIPTS_DIR = WORKSPACE_DIR / "scripts"
DATA_DIR = WORKSPACE_DIR / "data"
CACHE_DIR = DATA_DIR / "cache"
HYPOTHESES_DIR = DATA_DIR / "hypotheses"

# ── 并行执行 ──
DEFAULT_MAX_WORKERS = 6

# ── GEO ──
GEO_BASE_URL = "https://www.ncbi.nlm.nih.gov/geo"
GEO_DOWNLOAD_BASE = "https://www.ncbi.nlm.nih.gov/geo/download"
DEFAULT_PLATFORM = "GPL570"

# ── 可靠性阈值 ──
RELIABILITY_THRESHOLDS = {
    "min_datasets": 3,
    "min_combined_n": 10,
    "effect_size_strong": 0.8,
    "effect_size_moderate": 0.5,
    "effect_size_small": 0.2,
    "p_value_significant": 0.05,
    "p_value_highly_significant": 0.01,
    "i_squared_low": 25.0,
    "i_squared_moderate": 50.0,
    "i_squared_high": 75.0,
}

# ── GRADE 评级 ──
GRADE_CUTOFFS = {
    "🟢 High": 0.75,
    "🟡 Moderate": 0.50,
    "🟠 Low": 0.25,
}

# ── Phase ──
PHASE_ORDER = ["phase1", "phase2", "phase3", "phase4", "phase5", "phase6"]

# ── 数据缓存 TTL（秒） ──
CACHE_TTL_PUBMED = 86400      # 24h
CACHE_TTL_ANNOTATION = 604800 # 7d
CACHE_TTL_DOWNLOAD = 259200   # 3d

# ── API Keys（从环境变量读取） ──
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY", "")
PUBMED_EMAIL = os.environ.get("PUBMED_EMAIL", "bioanalysis@local")

# ── Registry 版本 ──
REGISTRY_VERSION = 3

# ── 基因集 ──
IMMUNE_CHECKPOINT_GENES = [
    "PDCD1", "PD-L1", "PD-L2", "HAVCR2", "LAG3", "TIGIT",
    "CTLA4", "BTLA", "CD160", "LAG3", "HAVCR2",
]
T_CELL_EXHAUSTION_MARKERS = [
    "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4",
    "TOX", "ENTPD1", "BATF", "IRF4",
]
MACROPHAGE_M2_MARKERS = [
    "CD163", "MRC1", "ARG1", "IL10", "TGFB1",
    "VEGFA", "CCL2", "CCL17", "CCL22",
]

# ── 细胞标记物 ──
CELL_TYPE_MARKERS = {
    "T_cells": ["CD3D", "CD3E", "CD8A", "CD4"],
    "Macrophages": ["CD68", "CD14", "AIF1"],
    "NK_cells": ["NKG7", "GNLY", "KLRD1"],
    "B_cells": ["CD19", "MS4A1", "CD79A"],
    "Fibroblasts": ["FAP", "PDGFRA", "ACTA2", "COL1A1"],
    "Tumor_cells": ["EPCAM", "KRT18", "KRT19"],
}
