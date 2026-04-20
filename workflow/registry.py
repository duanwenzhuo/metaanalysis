"""
Bio-Analysis State Registry

负责状态持久化、跨 Agent 共享、完整性验证。
所有 Agent 通过 Registry 读写状态，禁止直接操作文件系统。

设计原则：
- 内存缓存优先（避免重复读文件）
- 写时持久化（每个 Phase 完成后自动保存）
- 读时检查版本（防止状态损坏）
- 完整性校验（不允许写入未知字段）
"""

from __future__ import annotations
import json, shutil, os
from pathlib import Path
from datetime import datetime
from typing import Any, Optional

from workflow.state import BioAnalysisState, validate_state


WORKSPACE_DIR = Path(__file__).parent.parent  # bio-analysis/
DATA_DIR = WORKSPACE_DIR / "data" / "hypotheses"  # 新路径
REGISTRY_VERSION = 3  # 版本号，用于迁移


class StateRegistry:
    """
    共享状态管理器。
    
    用法：
        registry = StateRegistry("h_001_t_cell_exhaustion")
        
        # 读取（优先缓存）
        state = registry.get()
        
        # 更新（自动持久化）
        registry.set(state, phase="phase2")
        
        # 检查某 Phase 是否完成
        if registry.is_done("phase4"):
            results = state["phase4_results"]
    """

    def __init__(self, hypothesis_id: str):
        self.hypothesis_id = hypothesis_id
        self.hyp_dir = DATA_DIR / hypothesis_id
        self._cache: Optional[BioAnalysisState] = None
        self._dirty = False
        
        # 确保目录存在
        self.hyp_dir.mkdir(parents=True, exist_ok=True)
        (self.hyp_dir / "results").mkdir(exist_ok=True)

    # ── 路径 ─────────────────────────────────────

    @property
    def state_file(self) -> Path:
        return self.hyp_dir / "state.json"

    @property
    def log_file(self) -> Path:
        return self.hyp_dir / "_registry_log.txt"

    # ── 核心读写 ─────────────────────────────────

    def get(self, use_cache: bool = True) -> BioAnalysisState:
        """
        获取当前状态。
        优先返回内存缓存（避免重复读文件）。
        """
        if use_cache and self._cache is not None:
            return self._cache
        
        if not self.state_file.exists():
            raise FileNotFoundError(
                f"State not found for {self.hypothesis_id}. "
                f"Run Phase 1 first."
            )
        
        with open(self.state_file, encoding="utf-8") as f:
            raw = json.load(f)
        
        # 版本迁移（v1/v2 → v3）
        raw = self._migrate(raw)
        
        # 完整性验证
        errors = validate_state(raw)
        if errors:
            self._log(f"WARNING: state validation errors: {errors}")
        
        self._cache = raw
        return self._cache

    def set(
        self,
        state: BioAnalysisState,
        phase: Optional[str] = None,
        validate: bool = True,
    ) -> None:
        """
        更新状态并持久化。
        
        Args:
            state: 新的完整状态
            phase: 当前完成的 Phase（自动更新 checkpoints）
            validate: 是否校验状态合法性
        """
        if validate:
            errors = validate_state(state)
            if errors:
                raise ValueError(f"Invalid state: {errors}")
        
        # 自动更新时间戳
        state["updated_at"] = datetime.now().isoformat()
        
        # 如果 phase 完成，加入 checkpoints
        if phase:
            if "checkpoints" not in state:
                state["checkpoints"] = []
            if phase not in state["checkpoints"]:
                state["checkpoints"].append(phase)
            state["current_phase"] = phase
        
        # 写入文件
        output = {
            **state,
            "_version": REGISTRY_VERSION,
            "_hypothesis_id": self.hypothesis_id,
        }
        
        # 先写临时文件，再 rename（原子写入）
        tmp = self.state_file.with_suffix(".tmp")
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(output, f, ensure_ascii=False, indent=2)
        tmp.rename(self.state_file)
        
        # 更新缓存
        self._cache = state
        self._dirty = False
        
        self._log(f"Persisted state: phase={phase or 'update'}, "
                   f"hypothesis_id={self.hypothesis_id}")

    # ── 便捷方法 ─────────────────────────────────

    def is_done(self, phase: str) -> bool:
        """检查某 Phase 是否已完成"""
        try:
            state = self.get()
            return state.get("checkpoints", []) and phase in state.get("checkpoints", [])
        except FileNotFoundError:
            return False

    def is_failed(self, phase: str) -> bool:
        """检查某 Phase 是否失败"""
        try:
            state = self.get()
            return state.get(f"{phase}_status") == "failed"
        except FileNotFoundError:
            return False

    def get_phase_status(self, phase: str) -> str:
        """获取某 Phase 的状态"""
        try:
            return self.get().get(f"{phase}_status", "unknown")
        except FileNotFoundError:
            return "unknown"

    def get_checkpoint(self) -> Optional[str]:
        """获取最后一个完成的 Phase"""
        try:
            checkpoints = self.get().get("checkpoints", [])
            return checkpoints[-1] if checkpoints else None
        except FileNotFoundError:
            return None

    # ── Phase 输出文件便捷方法 ─────────────────────

    def read_phase2(self) -> Optional[dict]:
        path = self.hyp_dir / "phase2_gene_set.json"
        return json.loads(open(path).read()) if path.exists() else None

    def read_manifest(self) -> Optional[dict]:
        path = self.hyp_dir / "datasets_manifest.json"
        return json.loads(open(path).read()) if path.exists() else None

    def read_phase5(self) -> Optional[dict]:
        path = self.hyp_dir / "phase5_meta_result.json"
        return json.loads(open(path).read()) if path.exists() else None

    def list_phase4_results(self) -> list[dict]:
        """列出所有 Phase 4 结果文件"""
        results = []
        for f in sorted((self.hyp_dir / "results").glob("*_result.json")):
            try:
                d = json.loads(open(f).read())
                # 跳过元文件
                if "individual_results" in d and "combined_p_value" in d:
                    continue
                results.append(d)
            except Exception:
                pass
        return results

    # ── Phase 4 单数据集结果 ───────────────────────

    def write_dataset_result(self, gse_id: str, result: dict) -> None:
        """写入单个数据集的分析结果"""
        path = self.hyp_dir / "results" / f"{gse_id}_result.json"
        with open(path, "w", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
        
        # 立即更新 registry state
        state = self.get()
        if "phase4_results" not in state:
            state["phase4_results"] = {}
        state["phase4_results"][gse_id] = result
        self.set(state, phase=None)  # 不改变 current_phase

    # ── 初始化 ─────────────────────────────────

    @classmethod
    def init(cls, hypothesis_text: str) -> tuple[StateRegistry, str]:
        """
        创建新假设，返回 (registry, hypothesis_id)。
        
        等价于 Phase 1：解析假设 + 初始化状态。
        """
        from workflow.state import new_state, next_phase
        
        # 生成状态
        state = new_state(hypothesis_text)
        hyp_id = state["hypothesis_id"]
        
        # 创建 Registry 并保存
        registry = cls(hyp_id)
        registry.set(state, phase="hypothesis_parser")
        
        # 同时写 hypothesis.json（兼容旧脚本）
        hyp_path = registry.hyp_dir / "hypothesis.json"
        with open(hyp_path, "w", encoding="utf-8") as f:
            json.dump({
                "id": hyp_id,
                "hypothesis": hypothesis_text,
                "type": "Cell State",
                "analysis_strategy": "GSVA / AUCell scoring",
                "minimum_datasets": 12,
            }, f, ensure_ascii=False, indent=2)
        
        return registry, hyp_id

    @classmethod
    def load(cls, hypothesis_id: str) -> StateRegistry:
        """加载已有假设"""
        registry = cls(hypothesis_id)
        registry.get()  # 验证存在
        return registry

    # ── 内部工具 ─────────────────────────────────

    def _log(self, msg: str) -> None:
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] [Registry] {msg}"
        print(line)
        self.hyp_dir.mkdir(parents=True, exist_ok=True)
        with open(self.log_file, "a", encoding="utf-8") as f:
            f.write(line + "\n")

    def _migrate(self, raw: dict) -> dict:
        """迁移旧版本状态到当前版本"""
        version = raw.get("_version", 1)
        
        # v1 → v2: 添加 phase4_results 作为 dict
        if version < 2:
            raw["phase4_results"] = raw.get("phase4_results", {})
            raw["phase4_errors"] = raw.get("phase4_errors", [])
            version = 2
        
        # v2 → v3: 添加 checkpoints 字段
        if version < 3:
            if "checkpoints" not in raw:
                raw["checkpoints"] = []
            # 从 current_phase 推断 checkpoints
            if not raw["checkpoints"] and raw.get("current_phase"):
                phase_order = ["phase1","phase2","phase3","phase4","phase5","phase6"]
                try:
                    idx = phase_order.index(raw["current_phase"])
                    raw["checkpoints"] = phase_order[:idx+1]
                except ValueError:
                    pass
            version = 3
        
        raw["_version"] = REGISTRY_VERSION
        return raw

    def __repr__(self) -> str:
        try:
            s = self.get()
            cp = s.get("current_phase", "?")
            cps = s.get("checkpoints", [])
            return (f"StateRegistry({self.hypothesis_id}, "
                    f"current_phase={cp}, checkpoints={cps})")
        except Exception:
            return f"StateRegistry({self.hypothesis_id}, not-loaded)"

    @classmethod
    def list_all(cls) -> list:
        """列出 hypotheses 目录下所有假设（轻量，不加载完整 state）"""
        from core.constants import HYPOTHESES_DIR
        import json
        results = []
        if not HYPOTHESES_DIR.exists():
            return results
        for hyp_dir in HYPOTHESES_DIR.iterdir():
            if not hyp_dir.is_dir() or hyp_dir.name.startswith("."):
                continue
            state_file = hyp_dir / "state.json"
            if state_file.exists():
                try:
                    data = json.loads(state_file.read_text(encoding="utf-8"))
                    results.append({
                        "hypothesis_id": data.get("hypothesis_id", hyp_dir.name),
                        "hypothesis_text": data.get("hypothesis_text", ""),
                        "current_phase": data.get("current_phase", "?"),
                        "checkpoints": data.get("checkpoints", []),
                        "phase1_type": data.get("phase1_type", ""),
                        "phase5_reliability_score": data.get("phase5_reliability_score"),
                    })
                except Exception:
                    pass
        return results
