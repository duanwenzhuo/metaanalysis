"""
L0 原子技能层 — 核心基类与全局注册表

设计原则：
1. 每个 Skill 必须显式声明契约（input_contract / output_contract）
2. 执行前自动验证契约，执行后检查输出
3. 不抛异常：用 status=failed + error 字段
4. SkillRegistry 全局单例，支持热插拔
5. 状态传递：execute(state) → state'（不可变语义）
"""

from __future__ import annotations
import time
import hashlib
import json
import inspect
import importlib
from pathlib import Path
from typing import (
    TypedDict, Any, Optional, Dict, List, Set, Callable,
    Type, get_type_hints, get_origin, get_args,
    Literal,
)
from dataclasses import dataclass, field
from enum import Enum
import threading


# ─────────────────────────────────────────────
# Stage — 技能阶段枚举
# ─────────────────────────────────────────────

class Stage(Enum):
    """技能所属的处理阶段"""
    DATA = "data"                    # 数据读取/下载
    PREPROCESSING = "preprocessing"   # QC / 归一化 / HVG
    DIMENSION = "dimension"           # PCA / UMAP / Harmony
    CLUSTER = "cluster"               # Leiden / Louvain
    ANNOTATION = "annotation"         # Marker score / CellTypist
    ENRICHMENT = "enrichment"        # GSVA / AUCell / ssGSEA
    STATISTICS = "statistics"         # 效应量 / 检验 / Meta
    VISUALIZATION = "visualization"   # 绘图
    REPORTING = "reporting"           # 报告生成


class Modality(Enum):
    """数据模态"""
    SCRNA = "scrna"                   # 单细胞 RNA-seq
    BULKRNA = "bulk_rna"              # 批量 RNA-seq
    SPATIAL = "spatial"               # 空间转录组
    PROTEOMICS = "proteomics"         # 蛋白质组
    METABOLOMICS = "metabolomics"     # 代谢组
    METAGENOMICS = "metagenomics"     # 宏基因组
    MULTIOME = "multiome"             # 多组学
    GENERIC = "generic"                # 通用（任何模态）


StatusLiteral = Literal["pending", "success", "failed", "skipped"]


# ─────────────────────────────────────────────
# SkillInput / SkillOutput — 标准 I/O 类型
# ─────────────────────────────────────────────

class SkillInput(TypedDict, total=False):
    """标准 Skill 输入。技能从 State 中读取这些键。"""
    adata: Any                       # AnnData 对象（可选）
    state: "State"                   # 共享状态字典（必须）
    params: Dict[str, Any]           # 运行时参数（可选，有默认值）
    cache_path: Optional[str]        # 缓存路径（可选）


class SkillOutput(TypedDict, total=False):
    """标准 Skill 输出。技能返回新的 State 片段。"""
    state_updates: Dict[str, Any]    # 要 merge 回 State 的更新
    status: StatusLiteral
    error: Optional[str]
    logs: List[str]
    artifacts: Dict[str, Any]        # 产物路径或数据 {name: value}
    duration_seconds: float
    new_keys: List[str]              # 本次新增的 State 键


# ─────────────────────────────────────────────
# AbstractSkill — 所有原子技能的基类
# ─────────────────────────────────────────────

class AbstractSkill:
    """
    原子技能基类。
    
    子类必须定义：
    - name: str          唯一标识符
    - description: str   功能描述（用于 LLM 规划）
    - input_contract: List[str]   执行前必须存在于 State 的键
    - output_contract: List[str]   执行后将会添加到 State 的键
    
    可选定义：
    - modality: List[Modality]    适用模态
    - tunable_parameters: Dict     可调参数声明
    
    执行流程：
        execute(state) → validate → _run → check_output → return state_updates
    """
    
    name: str = "abstract_skill"
    description: str = "Abstract base skill"
    
    # 契约声明（子类必须覆盖）
    input_contract: List[str] = []   # 必须的 State 键
    output_contract: List[str] = []  # 执行后新增的 State 键
    optional_inputs: List[str] = []  # 可选的 State 键
    
    # 模态声明（默认通用）
    modality: List[Modality] = [Modality.GENERIC]
    
    # 可调参数声明
    tunable_parameters: Dict[str, Dict[str, Any]] = {}
    
    # 超时（秒）
    timeout_seconds: int = 300
    
    # ── 标准执行入口 ────────────────────────────
    
    def execute(self, state: State) -> SkillOutput:
        """
        标准执行入口：try/catch 自动包装。
        所有技能调用此方法，不直接调用 _run()。
        """
        start = time.time()
        output: SkillOutput = {
            "state_updates": {},
            "status": "pending",
            "logs": [],
            "artifacts": {},
            "duration_seconds": 0.0,
            "new_keys": [],
        }
        
        try:
            # 1. 契约验证
            missing = self._validate_input(state)
            if missing:
                raise ContractError(
                    f"[{self.name}] Missing required input keys: {missing}"
                )
            
            # 2. 执行核心逻辑
            self._log(output, f"[{self.name}] ▶ Starting...")
            result = self._run(state)
            
            # 3. 检查输出契约
            self._validate_output(result)
            
            output["state_updates"] = result
            output["status"] = "success"
            self._log(output, f"[{self.name}] ✅ Done")
            
        except ContractError as e:
            output["status"] = "failed"
            output["error"] = str(e)
            self._log(output, f"[{self.name}] ❌ ContractError: {e}")
            
        except Exception as e:
            output["status"] = "failed"
            output["error"] = f"{type(e).__name__}: {e}"
            self._log(output, f"[{self.name}] ❌ Exception: {e}")
        
        output["duration_seconds"] = round(time.time() - start, 3)
        return output
    
    # ── 子类实现 ────────────────────────────────
    
    def _run(self, state: State) -> Dict[str, Any]:
        """
        子类实现核心逻辑。
        返回 dict，将被 merge 到 State。
        """
        raise NotImplementedError(f"{self.name}._run() must be implemented")
    
    # ── 契约验证 ───────────────────────────────
    
    def _validate_input(self, state: State) -> List[str]:
        """返回缺失的必需输入键列表（空=全部满足）"""
        missing = []
        for key in self.input_contract:
            if key not in state:
                missing.append(key)
        return missing
    
    def _validate_output(self, result: Dict[str, Any]) -> None:
        """检查输出是否包含声明的 output_contract"""
        for key in self.output_contract:
            if key not in result:
                raise ContractError(
                    f"[{self.name}] Output contract violated: "
                    f"'{key}' declared but not produced"
                )
    
    # ── 工具方法 ────────────────────────────────
    
    def _log(self, output: SkillOutput, msg: str) -> None:
        """追加日志"""
        output["logs"].append(msg)
        print(msg)
    
    def _get_params(self, state: State) -> Dict[str, Any]:
        """从 state['params'] 或空字典安全获取参数"""
        return state.get("params", {})
    
    def _has_key(self, state: State, key: str) -> bool:
        """检查 State 中是否存在某键（含嵌套）"""
        if key in state:
            val = state[key]
            return val is not None and val != {}
        return False
    
    def _hash_inputs(self, state: State, keys: List[str]) -> str:
        """生成输入数据的哈希（用于缓存键）"""
        parts = []
        for k in keys:
            if k in state:
                val = state[k]
                try:
                    parts.append(str(hash(str(val)))[:16])
                except Exception:
                    parts.append(k)
        return hashlib.md5("_".join(parts).encode()).hexdigest()[:12]
    
    def compatible_with(self, state: State) -> bool:
        """检查当前 State 是否与此技能的模态兼容"""
        if Modality.GENERIC in self.modality:
            return True
        current_modality = state.get("_modality", Modality.GENERIC)
        return current_modality in self.modality


# ─────────────────────────────────────────────
# SkillExecutionError / ContractError
# ─────────────────────────────────────────────

class SkillExecutionError(Exception):
    """技能执行失败"""
    pass

class ContractError(Exception):
    """契约验证失败"""
    pass


# ─────────────────────────────────────────────
# State 类型别名
# ─────────────────────────────────────────────

State = Dict[str, Any]


# ─────────────────────────────────────────────
# SkillRegistry — 全局注册表（线程安全单例）
# ─────────────────────────────────────────────

class SkillRegistry:
    """
    全局技能注册表（单例，线程安全）。
    
    功能：
    - register(skill_cls): 注册技能类
    - get(name): 获取技能实例
    - list(): 列出所有技能
    - list_by_stage(stage): 按阶段过滤
    - list_by_modality(modality): 按模态过滤
    - find_producers(key): 查找能产出某 State 键的技能
    - find_consumers(key): 查找消费某 State 键的技能
    - auto_discover(): 自动发现 bioskills/ 下的所有技能
    """
    
    _instance: Optional["SkillRegistry"] = None
    _lock = threading.Lock()
    
    def __new__(cls) -> "SkillRegistry":
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance._init()
        return cls._instance
    
    def _init(self) -> None:
        self._skills: Dict[str, Type[AbstractSkill]] = {}
        self._by_stage: Dict[Stage, List[str]] = {s: [] for s in Stage}
        self._by_modality: Dict[Modality, List[str]] = {m: [] for m in Modality}
        self._reverse_input: Dict[str, List[str]] = {}  # key → [skill_names]
        self._reverse_output: Dict[str, List[str]] = {}  # key → [skill_names]
        self._initialized: bool = False
    
    @classmethod
    def reset(cls) -> None:
        """重置注册表（用于测试）"""
        with cls._lock:
            cls._instance = None
    
    # ── 注册 ──────────────────────────────────
    
    def register(self, cls: Type[AbstractSkill]) -> Type[AbstractSkill]:
        """注册技能类（支持装饰器 @register，幂等）"""
        name = cls.name
        
        if name in self._skills:
            # 幂等：已注册则跳过（装饰器在模块加载时已注册）
            return cls
        
        self._skills[name] = cls
        
        # 按 Stage 索引
        if hasattr(cls, "stage"):
            stage = getattr(cls, "stage")
            if isinstance(stage, Stage):
                if name not in self._by_stage[stage]:
                    self._by_stage[stage].append(name)
        
        # 按 Modality 索引
        modalities = getattr(cls, "modality", [Modality.GENERIC])
        for m in modalities:
            if name not in self._by_modality[m]:
                self._by_modality[m].append(name)
        
        # 逆向索引：key → consumers（谁需要这个 key）
        for key in getattr(cls, "input_contract", []):
            self._reverse_input.setdefault(key, []).append(name)
        
        # 逆向索引：key → producers（谁能产出这个 key）
        for key in getattr(cls, "output_contract", []):
            self._reverse_output.setdefault(key, []).append(name)
        
        print(f"[SkillRegistry] ✅ Registered: {name}")
        return cls
    
    # ── 查询 ──────────────────────────────────
    
    def get(self, name: str) -> AbstractSkill:
        """获取技能实例（每次返回新实例，保证无状态）"""
        if name not in self._skills:
            available = list(self._skills.keys())
            raise KeyError(
                f"Skill '{name}' not found. Available: {available}"
            )
        return self._skills[name]()
    
    def get_class(self, name: str) -> Type[AbstractSkill]:
        """获取技能类（不实例化）"""
        if name not in self._skills:
            raise KeyError(f"Skill '{name}' not found")
        return self._skills[name]
    
    def list(self) -> List[str]:
        """列出所有技能名称"""
        return list(self._skills.keys())
    
    def list_by_stage(self, stage: Stage) -> List[str]:
        """按阶段列出技能"""
        return list(self._by_stage.get(stage, []))
    
    def list_by_modality(self, modality: Modality) -> List[str]:
        """按模态列出技能"""
        return list(self._by_modality.get(modality, []))
    
    def find_producers(self, key: str) -> List[str]:
        """查找能产出某 State 键的技能"""
        return list(self._reverse_output.get(key, []))
    
    def find_consumers(self, key: str) -> List[str]:
        """查找需要某 State 键的技能"""
        return list(self._reverse_input.get(key, []))
    
    def find_by_contract(
        self,
        needs: List[str],
        produces: Optional[List[str]] = None,
    ) -> List[str]:
        """
        按契约查找技能。
        
        Args:
            needs: 需要的 State 键
            produces: 需要产出的 State 键（可选）
        
        Returns:
            匹配的技能名称列表
        """
        candidates = set(self._skills.keys())
        
        # 必须满足所有 needs
        for key in needs:
            consumers = set(self._reverse_input.get(key, []))
            candidates &= consumers
        
        # 如果指定了 produces，必须匹配
        if produces:
            results = []
            for name in candidates:
                skill = self._skills[name]
                if all(k in skill.output_contract for k in produces):
                    results.append(name)
            return results
        
        return list(candidates)
    
    # ── 自动发现 ────────────────────────────────
    
    def auto_discover(self, base_dir: Optional[Path] = None) -> int:
        """
        自动发现并注册 bioskills/ 下的所有技能。
        
        扫描规则：
        - bioskills/<category>/<name>.py
        - 类名以 "Skill" 结尾或函数名以 "skill_" 开头
        """
        if base_dir is None:
            base_dir = Path(__file__).parent.parent
            # 在 dev 环境中找到 bioskills 目录
            for p in [base_dir / "bioskills", base_dir.parent / "bioskills"]:
                if p.exists():
                    base_dir = p.parent
                    break
        
        bioskills_dir = Path(__file__).parent.parent / "bioskills"
        if not bioskills_dir.exists():
            bioskills_dir = Path(__file__).parent.parent
        
        count = 0
        for category_dir in bioskills_dir.iterdir():
            if not category_dir.is_dir() or category_dir.name.startswith("_"):
                continue
            
            for py_file in category_dir.glob("*.py"):
                if py_file.name.startswith("_"):
                    continue
                
                try:
                    module_name = (
                        f"bioskills.{category_dir.name}.{py_file.stem}"
                    )
                    # 尝试从已加载的模块导入
                    try:
                        mod = importlib.import_module(module_name)
                    except ImportError:
                        # 可能路径不对，尝试相对导入
                        continue
                    
                    # 查找 Skill 类
                    for attr_name in dir(mod):
                        if attr_name.startswith("_"):
                            continue
                        attr = getattr(mod, attr_name)
                        if (
                            isinstance(attr, type)
                            and issubclass(attr, AbstractSkill)
                            and attr is not AbstractSkill
                        ):
                            self.register(attr)
                            count += 1
                
                except Exception as e:
                    print(f"[SkillRegistry] ⚠️  Failed to load {py_file}: {e}")
        
        self._initialized = True
        print(f"[SkillRegistry] 🔍 Auto-discovered {count} skills")
        return count
    
    def summary(self) -> Dict[str, Any]:
        """返回注册表摘要（用于调试）"""
        return {
            "total": len(self._skills),
            "by_stage": {s.value: len(v) for s, v in self._by_stage.items()},
            "by_modality": {m.value: len(v) for m, v in self._by_modality.items()},
            "skills": sorted(self._skills.keys()),
        }


# ─────────────────────────────────────────────
# register 装饰器
# ─────────────────────────────────────────────

def register(cls: Type[AbstractSkill]) -> Type[AbstractSkill]:
    """
    技能注册装饰器。
    
    用法:
        @register
        class MySkill(AbstractSkill):
            name = "my_skill"
            ...
    
    等价于:
        class MySkill(AbstractSkill):
            ...
        SkillRegistry().register(MySkill)
    """
    registry = SkillRegistry()
    registry.register(cls)
    return cls


# ─────────────────────────────────────────────
# 默认导出 convenience 函数
# ─────────────────────────────────────────────

def list_skills() -> List[str]:
    """列出所有注册的技能"""
    return SkillRegistry().list()

def get_skill(name: str) -> AbstractSkill:
    """获取技能实例"""
    return SkillRegistry().get(name)
