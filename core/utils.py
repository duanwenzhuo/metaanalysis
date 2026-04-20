"""core/utils.py — 通用工具函数。"""

from __future__ import annotations
import json, hashlib, os, sys
from pathlib import Path
from typing import Any, Optional
from datetime import datetime


# ─────────────────────────────────────────────
# JSON 工具
# ─────────────────────────────────────────────

class NumpyEncoder(json.JSONEncoder):
    """支持 numpy 类型的 JSON 编码器"""
    def default(self, obj):
        try:
            import numpy as np
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.bool_):
                return bool(obj)
        except ImportError:
            pass
        return super().default(obj)


def safe_json_dump(obj: Any, fp, **kwargs) -> None:
    """安全的 JSON 写入（自动使用 NumpyEncoder）"""
    kwargs.setdefault("ensure_ascii", False)
    kwargs.setdefault("indent", 2)
    encoder = NumpyEncoder
    json.dump(obj, fp, cls=encoder, **kwargs)


def safe_json_load(fp_or_path) -> dict:
    """安全的 JSON 读取"""
    if isinstance(fp_or_path, (str, Path)):
        with open(fp_or_path, encoding="utf-8") as f:
            return json.load(f)
    return json.load(fp_or_path)


# ─────────────────────────────────────────────
# 文件系统工具
# ─────────────────────────────────────────────

def ensure_dir(path: str | Path, mode: int = 0o755) -> Path:
    """确保目录存在（不存在则创建）"""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True, mode=mode)
    return p


def compute_sha1(file_path: str | Path, chunk_size: int = 8192) -> str:
    """计算文件 SHA1"""
    h = hashlib.sha1()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def get_file_age_hours(path: str | Path) -> float:
    """获取文件年龄（小时）"""
    mtime = os.path.getmtime(path)
    return (datetime.now().timestamp() - mtime) / 3600.0


def is_cache_stale(path: str | Path, ttl_seconds: float) -> bool:
    """判断缓存是否过期"""
    if not os.path.exists(path):
        return True
    return get_file_age_hours(path) * 3600 > ttl_seconds


# ─────────────────────────────────────────────
# 字符串工具
# ─────────────────────────────────────────────

def truncate(text: str, max_len: int = 200, suffix: str = "...") -> str:
    """截断字符串"""
    if len(text) <= max_len:
        return text
    return text[: max_len - len(suffix)] + suffix


def sanitize_filename(name: str) -> str:
    """将任意字符串转换为合法的文件名"""
    import re
    name = re.sub(r"[^\w\-\.]", "_", name)
    name = re.sub(r"_+", "_", name).strip("._")
    return name[: 200] or "unnamed"


# ─────────────────────────────────────────────
# 时间工具
# ─────────────────────────────────────────────

def now_iso() -> str:
    """返回当前 UTC ISO 时间字符串"""
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat()


def parse_year(year_str: str) -> int:
    """解析年份字符串"""
    import re
    m = re.search(r"\d{4}", str(year_str))
    return int(m.group()) if m else 0


# ─────────────────────────────────────────────
# 统计工具
# ─────────────────────────────────────────────

def cohens_d(mean1: float, mean2: float, std1: float, std2: float,
             n1: int, n2: int) -> float:
    """计算 Cohen's d 效应量"""
    import math
    pooled_std = math.sqrt(
        ((n1 - 1) * std1 ** 2 + (n2 - 1) * std2 ** 2) / (n1 + n2 - 2)
    )
    if pooled_std == 0:
        return 0.0
    return (mean1 - mean2) / pooled_std


def z_to_p(z: float) -> float:
    """Z 值转 P 值（双尾）"""
    import math
    return 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2))))


# ─────────────────────────────────────────────
# GEO / 数据验证工具
# ─────────────────────────────────────────────

def validate_gse_id(gse_id: str) -> bool:
    """验证 GSE ID 格式"""
    import re
    if not gse_id:
        return False
    # GSE + 数字，常见长度 6-8 位
    return bool(re.match(r"^GSE\d{5,8}$", gse_id, re.IGNORECASE))


def format_effect_size(d: float) -> str:
    """格式化效应量显示"""
    if abs(d) < 0.2:
        return f"{d:.3f} (negligible)"
    elif abs(d) < 0.5:
        return f"{d:.3f} (small)"
    elif abs(d) < 0.8:
        return f"{d:.3f} (medium)"
    else:
        return f"{d:.3f} (large)"


def format_confidence_interval(lower: float, upper: float, decimals: int = 2) -> str:
    """格式化置信区间"""
    return f"[{lower:.{decimals}f}, {upper:.{decimals}f}]"


def grade_from_score(score: float) -> str:
    """根据可靠性评分返回等级"""
    if score >= 0.8:
        return "HIGH"
    elif score >= 0.6:
        return "MODERATE"
    elif score >= 0.4:
        return "LOW"
    else:
        return "INSUFFICIENT"


def emoji_for_grade(grade: str) -> str:
    """等级对应的 emoji"""
    return {"HIGH": "🟢", "MODERATE": "🟡", "LOW": "🟠", "INSUFFICIENT": "⚪"}.get(grade, "⚪")
