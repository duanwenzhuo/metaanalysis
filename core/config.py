"""core/config.py — 配置管理（环境变量 + 配置文件）。"""

from __future__ import annotations
import os, json
from pathlib import Path
from typing import Optional, Dict, Any
from core.exceptions import ConfigurationError


class Config:
    """
    全局配置单例。
    优先从环境变量读取，其次从配置文件读取，最后使用默认值。
    """

    _instance: Optional["Config"] = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._init()
        return cls._instance

    def _init(self):
        from core.constants import WORKSPACE_DIR, DEFAULT_MAX_WORKERS
        self.WORKSPACE_DIR = WORKSPACE_DIR
        self.DATA_DIR = WORKSPACE_DIR / "data"
        self.MAX_WORKERS = int(
            os.environ.get("BIO_MAX_WORKERS", DEFAULT_MAX_WORKERS)
        )

        # API Keys
        self.NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
        self.PUBMED_EMAIL = os.environ.get("PUBMED_EMAIL", "bioanalysis@local")
        self.OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY", "")

        # 路径配置
        self.CACHE_DIR = self.DATA_DIR / "cache"
        self.HYPOTHESES_DIR = self.DATA_DIR / "hypotheses"

        # 阈值
        self.MIN_DATASETS = int(os.environ.get("BIO_MIN_DATASETS", "3"))
        self.EFFECT_SIZE_STRONG = float(
            os.environ.get("BIO_EFFECT_SIZE_STRONG", "0.8")
        )

        # 调试
        self.DEBUG = os.environ.get("BIO_DEBUG", "0") == "1"

    def get(self, key: str, default: Any = None) -> Any:
        return getattr(self, key, default)

    def to_dict(self) -> Dict[str, Any]:
        """返回所有配置（不含敏感信息）"""
        return {
            k: v for k, v in self.__dict__.items()
            if not k.startswith("_") and "KEY" not in k
        }

    @classmethod
    def load_from_file(cls, config_path: str | Path) -> None:
        """从 JSON 配置文件加载配置"""
        p = Path(config_path)
        if not p.exists():
            raise ConfigurationError(f"Config file not found: {p}")
        with open(p, encoding="utf-8") as f:
            data = json.load(f)
        # 写入环境变量（后续 Config() 会自动读取）
        for k, v in data.items():
            os.environ.setdefault(f"BIO_{k.upper()}", str(v))


# 便捷单例访问
def get_config() -> Config:
    return Config()
