from __future__ import annotations

DEFAULT_SKILL_TEMPLATE = '''"""L2 Synthesized Skill: {skill_name}
{description}

⚠️  AUTO-GENERATED CODE — Review before production use
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class {class_name}(AbstractSkill):
    name = "{skill_name}"
    description = "{description}"
    input_contract = {input_contract}
    output_contract = {output_contract}
    stage = Stage.DATA
    modality = [Modality.GENERIC]

    def _run(self, state: State) -> dict:
        print(f"[{self.name}] ⚠️  Placeholder — implement actual logic")
        return {{"{missing_key}": None}}
'''

GSVA_TEMPLATE = '''"""L2 Synthesized Skill: {skill_name}
GSVA / Gene Set Enrichment Analysis Skill
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
import numpy as np

@register
class {class_name}(AbstractSkill):
    name = "{skill_name}"
    description = "Gene Set Variation Analysis (GSVA) for gene set scoring"
    input_contract = ["adata", "params"]
    output_contract = ["gsva_scores", "adata"]
    stage = Stage.ENRICHMENT
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.GENERIC]

    def _run(self, state: State) -> dict:
        import scanpy as sc
        import gseapy as gp
        adata = state["adata"]
        params = state.get("params", {{}})
        gene_sets = params.get("gene_sets", {{"CUSTOM": []}})
        expr = adata.to_df().T
        available_gs = {{}}
        for gs_name, genes in gene_sets.items():
            avail = [g for g in genes if g in expr.index]
            if len(avail) >= 5:
                available_gs[gs_name] = avail
        if not available_gs:
            return {{"gsva_scores": {{}}, "adata": adata}}
        res = gp.ssgsea(
            data=expr, gene_sets=available_gs,
            sample_norm_method=params.get("sample_norm", "rank"),
            min_size=5, max_size=500, no_plot=True, outdir=None,
        )
        scores = res.resultsOnSamples.T
        for gs_name, row in scores.iterrows():
            adata.obs[f"gsva_{{gs_name}}"] = adata.obs_names.map(row.to_dict()).values
        return {{"gsva_scores": scores.to_dict(), "adata": adata}}
'''

SCALE_TEMPLATE = '''"""L2 Synthesized Skill: {skill_name}
Scale Skill
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class {class_name}(AbstractSkill):
    name = "{skill_name}"
    description = "Z-score scaling per gene"
    input_contract = ["adata"]
    output_contract = ["adata"]
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.GENERIC]

    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {{}})
        sc.pp.scale(adata, max_value=params.get("max_value", 10.0))
        return {{"adata": adata}}
'''
NORMALIZE_TEMPLATE = '''"""L2 Synthesized Skill: {skill_name}
Normalization Skill
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class {class_name}(AbstractSkill):
    name = "{skill_name}"
    description = "Expression matrix normalization (total count + log1p)"
    input_contract = ["adata"]
    output_contract = ["adata"]
    stage = Stage.PREPROCESSING
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.GENERIC]

    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {{}})
        target_sum = params.get("normalize_target_sum", 1e4)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        return {{"adata": adata}}
'''



CLUSTER_TEMPLATE = '''"""L2 Synthesized Skill: {skill_name}
Clustering Skill
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class {class_name}(AbstractSkill):
    name = "{skill_name}"
    description = "Cell clustering (Leiden/Louvain)"
    input_contract = ["adata"]
    output_contract = ["clusters", "adata"]
    stage = Stage.CLUSTER
    modality = [Modality.SCRNA, Modality.SPATIAL, Modality.GENERIC]

    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"]
        params = state.get("params", {{}})
        resolution = params.get("leiden_resolution", 0.5)
        if "neighbors" not in adata.uns:
            sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.leiden(adata, resolution=resolution, key_added="clusters")
        return {{"clusters": adata.obs["clusters"].to_dict(), "adata": adata}}
'''

AUCELL_TEMPLATE = '''"""L2 Synthesized Skill: {skill_name}
AUCell Skill
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class {class_name}(AbstractSkill):
    name = "{skill_name}"
    description = "AUCell gene set enrichment scoring"
    input_contract = ["adata", "params"]
    output_contract = ["gsva_scores", "adata"]
    stage = Stage.ENRICHMENT
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.GENERIC]

    def _run(self, state: State) -> dict:
        import gseapy as gp
        adata = state["adata"]
        params = state.get("params", {{}})
        gene_sets = params.get("gene_sets", {{"CUSTOM": []}})
        expr = adata.to_df().T
        available_gs = {{}}
        for gs_name, genes in gene_sets.items():
            avail = [g for g in genes if g in expr.index]
            if len(avail) >= 5:
                available_gs[gs_name] = avail
        if not available_gs:
            return {{"gsva_scores": {{}}, "adata": adata}}
        res = gp.aucell(expr, gene_sets=available_gs, min_size=5, max_size=500)
        for gs_name, row in res.iterrows():
            adata.obs[f"gsva_{{gs_name}}"] = adata.obs_names.map(row.to_dict()).values
        return {{"gsva_scores": res.to_dict(), "adata": adata}}
'''

EMBEDDING_TEMPLATE = '''"""L2 Synthesized Skill: {skill_name}
Embedding Skill
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class {class_name}(AbstractSkill):
    name = "{skill_name}"
    description = "Dimensionality reduction embedding"
    input_contract = ["adata"]
    output_contract = ["adata"]
    stage = Stage.DIMENSION
    modality = [Modality.SCRNA, Modality.SPATIAL, Modality.GENERIC]

    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"].copy()
        params = state.get("params", {{}})
        n_pcs = params.get("n_pcs", 50)
        if "X_pca" not in adata.obsm:
            sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars-1))
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]))
        sc.tl.umap(adata)
        return {{"adata": adata}}
'''

ANNOTATION_TEMPLATE = '''"""L2 Synthesized Skill: {skill_name}
Cell Annotation Skill
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class {class_name}(AbstractSkill):
    name = "{skill_name}"
    description = "Cell type annotation by marker genes"
    input_contract = ["adata"]
    output_contract = ["adata"]
    stage = Stage.ANNOTATION
    modality = [Modality.SCRNA, Modality.SPATIAL, Modality.GENERIC]

    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"]
        params = state.get("params", {{}})
        markers = params.get("markers", [])
        cell_type = params.get("cell_type", "Unknown")
        if not markers:
            return {{"adata": adata}}
        available = [g for g in markers if g in adata.var_names]
        score_name = f"score_{{cell_type.lower().replace(\" \", \"_\")}}"
        sc.tl.score_genes(adata, available, score_name=score_name, copy=False)
        return {{"adata": adata}}
'''
DIFFERENTIAL_EXPRESSION_TEMPLATE = '''"""L2 Synthesized Skill: {skill_name}
Differential Expression Skill
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class {class_name}(AbstractSkill):
    name = "{skill_name}"
    description = "Differential gene expression analysis"
    input_contract = ["adata"]
    output_contract = ["de_results", "adata"]
    stage = Stage.STATISTICS
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.GENERIC]

    def _run(self, state: State) -> dict:
        import scanpy as sc
        import pandas as pd
        adata = state["adata"]
        params = state.get("params", {{}})

        group_col = params.get("group_col", "clusters")
        if group_col not in adata.obs.columns:
            for alt in ["cell_type", "annotation", "batch"]:
                if alt in adata.obs.columns:
                    group_col = alt
                    break

        if group_col not in adata.obs.columns:
            return {{"de_results": {{"error": "no_group_column"}}, "adata": adata}}

        method = params.get("method", "wilcoxon")
        sc.tl.rank_genes_groups(adata, groupby=group_col, method=method)

        de_results = {{}}
        groups = adata.uns["rank_genes_groups"]["names"].dtype.names
        for g in groups:
            genes = adata.uns["rank_genes_groups"]["names"][g]
            scores = adata.uns["rank_genes_groups"]["scores"][g]
            pvals = adata.uns["rank_genes_groups"]["pvals"][g]
            de_results[g] = pd.DataFrame({{
                "gene": genes,
                "score": scores,
                "pval": pvals,
            }})

        return {{"de_results": de_results, "adata": adata}}
'''



MARKER_FIND_TEMPLATE = '''"""L2 Synthesized Skill: {skill_name}
Marker Gene Finding Skill
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State

@register
class {class_name}(AbstractSkill):
    name = "{skill_name}"
    description = "Find marker genes for cell clusters"
    input_contract = ["adata"]
    output_contract = ["markers"]
    stage = Stage.ANNOTATION
    modality = [Modality.SCRNA, Modality.GENERIC]

    def _run(self, state: State) -> dict:
        import scanpy as sc
        adata = state["adata"]
        params = state.get("params", {{}})
        cluster_col = params.get("cluster_col", "clusters")
        if cluster_col not in adata.obs.columns:
            return {{"markers": {{}}}}
        sc.tl.rank_genes_groups(adata, groupby=cluster_col, method="wilcoxon")
        marker_dict = {{}}
        for g in adata.uns["rank_genes_groups"]["names"]:
            for name in g.dtype.names or []:
                if name not in marker_dict:
                    marker_dict[name] = []
                marker_dict[name].append(g[name])
        return {{"markers": marker_dict}}
'''


"""
L2 自愈层 — 能力缺口分析与自动合成

三层工作流：
    GapAnalyzer → SkillSynthesizer → SelfHealingPipeline

使用方式:
    l2 = L2SelfHealingPipeline(registry=registry, failed_result=p1_result)
    result = l2.execute(state)
"""

# from __future__ import annotations  # moved to line 1
import os, sys, time, json, re
from pathlib import Path
from typing import List, Dict, Any, Optional, Set, Tuple
from dataclasses import dataclass, field
from collections import defaultdict

from bioskills.core.pipeline import PipelineResult as PIPResult
from bioskills.core.base import (
    SkillRegistry, AbstractSkill, State,
    StatusLiteral,
)
from bioskills.core.contract import (
    ContractError, ContractValidator,
)


# ─────────────────────────────────────────────
# 常量
# ─────────────────────────────────────────────

# 合成技能存放目录
SKILLS_DIR = Path(__file__).parent.parent
SYNTHESIZED_DIR = SKILLS_DIR / "_synthesized"
CONTRACTS_DIR = SYNTHESIZED_DIR / "contracts"
CODE_TEMPLATES_DIR = SYNTHESIZED_DIR / "templates"

SYNTHESIZED_DIR.mkdir(parents=True, exist_ok=True)
CONTRACTS_DIR.mkdir(parents=True, exist_ok=True)
CODE_TEMPLATES_DIR.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────
# GapAnalyzer — 缺口分析
# ─────────────────────────────────────────────

@dataclass
class SkillGap:
    """单个技能缺口"""
    missing_key: str               # 缺失的输出键
    needed_by: List[str]            # 需要这个键的技能
    suggested_inputs: List[str]    # 推断的输入键
    suggested_description: str     # LLM 生成描述
    gap_severity: str               # critical / warning




class GapAnalyzer:
    """
    能力缺口分析器。
    
    分析当前注册表中无法通过已有技能产出的 State 键，
    为 SkillSynthesizer 提供生成依据。
    """
    
    def __init__(self, registry: Optional[SkillRegistry] = None):
        self.registry = registry or SkillRegistry()
        self.validator = ContractValidator(self.registry)
    
    def analyze(
        self,
        start_keys: Set[str],
        goal_keys: Set[str],
        failed_skills: Optional[List[str]] = None,
    ) -> List[SkillGap]:
        """
        分析能力缺口。
        
        Args:
            start_keys: 当前 State 中已就绪的键
            goal_keys: 目标 State 键
            failed_skills: L1 层失败的技能列表
        
        Returns:
            List[SkillGap] — 所有无法填补的缺口
        """
        gaps: List[SkillGap] = []
        covered = start_keys.copy()
        remaining_goals = goal_keys.copy()
        
        # 迭代扩展覆盖范围，直到无法再扩展
        iteration = 0
        max_iter = len(self.registry.list()) + 5
        
        while remaining_goals and iteration < max_iter:
            iteration += 1
            made_progress = False
            
            for goal_key in list(remaining_goals):
                producers = self.registry.find_producers(goal_key)
                
                if not producers:
                    # 找不到生产者 → 缺口
                    consumers = self.registry.find_consumers(goal_key)
                    inputs = self._infer_inputs(goal_key, consumers)
                    
                    gaps.append(SkillGap(
                        missing_key=goal_key,
                        needed_by=consumers[:3],
                        suggested_inputs=inputs,
                        suggested_description=self._describe_gap(goal_key, inputs),
                        gap_severity="critical" if goal_key in goal_keys else "warning",
                    ))
                    remaining_goals.discard(goal_key)
                
                else:
                    # 找到了生产者，将其输出加入覆盖
                    for sn in producers:
                        if failed_skills and sn in failed_skills:
                            continue  # 已知失败，跳过
                        skill_cls = self.registry.get_class(sn)
                        for out_k in skill_cls.output_contract:
                            if out_k not in covered:
                                covered.add(out_k)
                                remaining_goals.discard(out_k)
                                made_progress = True
            
            if not made_progress and remaining_goals:
                # 无法继续扩展 → 剩余都是缺口
                for goal_key in remaining_goals:
                    consumers = self.registry.find_consumers(goal_key)
                    gaps.append(SkillGap(
                        missing_key=goal_key,
                        needed_by=consumers[:3],
                        suggested_inputs=[],
                        suggested_description=f"需要能产出 '{goal_key}' 的技能",
                        gap_severity="critical",
                    ))
                break
        
        return gaps
    
    def _infer_inputs(self, key: str, consumers: List[str]) -> List[str]:
        """根据消费者推断缺失技能的输入"""
        inputs_sets = []
        for sn in consumers:
            try:
                skill_cls = self.registry.get_class(sn)
                inputs_sets.append(set(skill_cls.input_contract))
            except Exception:
                continue
        
        if not inputs_sets:
            return []
        
        # 取所有消费者的输入交集（最保守）
        inferred = set.intersection(*inputs_sets) if inputs_sets else set()
        # 过滤掉 key 自身
        return sorted(inferred - {key})
    
    def _describe_gap(self, key: str, inputs: List[str]) -> str:
        """为缺口生成自然语言描述"""
        input_str = ", ".join(inputs) if inputs else "?"
        return (
            f"需要一个新的技能：将 {input_str} 转换为 '{key}'。"
            f"参考命名：custom_{key}_skill"
        )
    
    def generate_report(self, gaps: List[SkillGap]) -> str:
        """生成缺口报告"""
        lines = ["=" * 60, "🔍 L2 Gap Analysis Report", "=" * 60]
        
        critical = [g for g in gaps if g.gap_severity == "critical"]
        warning = [g for g in gaps if g.gap_severity == "warning"]
        
        lines.append(f"\n⛔ Critical ({len(critical)}):")
        for g in critical:
            lines.append(f"  • {g.missing_key}")
            lines.append(f"    ← needed by: {', '.join(g.needed_by)}")
            lines.append(f"    ← inputs: {g.suggested_inputs}")
        
        if warning:
            lines.append(f"\n⚠️  Warning ({len(warning)}):")
            for g in warning:
                lines.append(f"  • {g.missing_key}")
        
        lines.append("=" * 60)
        return "\n".join(lines)


# ─────────────────────────────────────────────
# SkillSynthesizer — 技能自动合成
# ─────────────────────────────────────────────

@dataclass
class SynthesizedSkill:
    """合成的技能元信息"""
    skill_name: str
    class_name: str
    file_path: Path
    input_contract: List[str]
    output_contract: List[str]
    description: str
    code: str
    status: str = "pending"    # pending / compiled / loaded


class SkillSynthesizer:
    """
    技能代码自动合成器。
    
    两种模式：
    1. template_mode: 基于代码模板生成（无 LLM API 依赖）
    2. llm_mode: 调用 LLM 生成（更智能，需要 API）
    
    生成流程：
    1. 分析缺口信息
    2. 选择代码模板
    3. 填充占位符
    4. 写入 .py 文件
    5. 验证语法
    6. 热加载到注册表
    """
    
    def __init__(
        self,
        registry: Optional[SkillRegistry] = None,
        mode: str = "template",   # "template" | "llm"
        llm_config: Optional[Dict] = None,
    ):
        self.registry = registry or SkillRegistry()
        self.mode = mode
        self.llm_config = llm_config or {}
    
    def synthesize(self, gap: SkillGap) -> Optional[SynthesizedSkill]:
        """
        根据缺口合成新技能。
        
        策略选择：
        - 如果 gap 有标准模板（如 gsva_*, normalize_*）→ 用模板
        - 否则用通用占位符模板
        """
        # 生成技能名
        skill_name = self._make_skill_name(gap.missing_key)
        class_name = "".join(
            p.capitalize() for p in re.split(r"[_\-]", gap.missing_key)
        ) + "Skill"
        
        # 选择模板
        template = self._select_template(gap)
        
        # 填充代码
        code = self._fill_template(
            template=template,
            skill_name=skill_name,
            class_name=class_name,
            gap=gap,
        )
        
        # 写入文件
        file_path = SYNTHESIZED_DIR / f"{skill_name}.py"
        file_path.write_text(code, encoding="utf-8")
        
        return SynthesizedSkill(
            skill_name=skill_name,
            class_name=class_name,
            file_path=file_path,
            input_contract=gap.suggested_inputs,
            output_contract=[gap.missing_key],
            description=gap.suggested_description,
            code=code,
        )
    
    def synthesize_batch(
        self, gaps: List[SkillGap]
    ) -> List[SynthesizedSkill]:
        """批量合成多个缺口"""
        results = []
        for gap in gaps:
            result = self.synthesize(gap)
            if result:
                results.append(result)
        return results
    
    def load_synthesized(self, skill: SynthesizedSkill) -> bool:
        """
        热加载合成的技能到注册表。
        
        流程：
        1. 验证语法（compile）
        2. importlib 动态导入
        3. 注册到 SkillRegistry
        """
        try:
            # 语法验证
            compile(skill.code, str(skill.file_path), "exec")
            
            # 动态导入
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                skill.skill_name, skill.file_path
            )
            if spec is None or spec.loader is None:
                return False
            
            module = importlib.util.module_from_spec(spec)
            sys.modules[skill.skill_name] = module
            spec.loader.exec_module(module)
            
            # 查找 Skill 类
            for attr_name in dir(module):
                if attr_name.startswith("_"):
                    continue
                attr = getattr(module, attr_name)
                if (
                    isinstance(attr, type)
                    and issubclass(attr, AbstractSkill)
                    and attr is not AbstractSkill
                ):
                    self.registry.register(attr)
                    skill.status = "loaded"
                    print(f"[SkillSynthesizer] ✅ Loaded: {attr_name}")
                    return True
            
            return False
            
        except SyntaxError as e:
            print(f"[SkillSynthesizer] ❌ SyntaxError in {skill.skill_name}: {e}")
            skill.status = "failed"
            return False
        except Exception as e:
            print(f"[SkillSynthesizer] ❌ Failed to load {skill.skill_name}: {e}")
            skill.status = "failed"
            return False
    
    def _make_skill_name(self, key: str) -> str:
        """从输出键生成技能名"""
        safe = re.sub(r"[^a-z0-9_]", "_", key.lower())
        safe = re.sub(r"_+", "_", safe).strip("_")
        return f"custom_{safe}"
    
    def _select_template(self, gap: SkillGap) -> str:
        """根据缺口选择代码模板"""
        key = gap.missing_key.lower()
        
        templates = {
            "gsva": GSVA_TEMPLATE,
            "aucell": AUCELL_TEMPLATE,
            "normalize": NORMALIZE_TEMPLATE,
            "scale": SCALE_TEMPLATE,
            "cluster": CLUSTER_TEMPLATE,
            "embed": EMBEDDING_TEMPLATE,
            "annot": ANNOTATION_TEMPLATE,
            "de": DIFFERENTIAL_EXPRESSION_TEMPLATE,
            "marker": MARKER_FIND_TEMPLATE,
            "default": DEFAULT_SKILL_TEMPLATE,
        }
        
        for keyword, tmpl in templates.items():
            if keyword in key:
                return tmpl
        return DEFAULT_SKILL_TEMPLATE
    
    def _fill_template(
        self,
        template: str,
        skill_name: str,
        class_name: str,
        gap: SkillGap,
    ) -> str:
        """填充模板"""
        inputs_str = ", ".join(f'"{k}"' for k in gap.suggested_inputs)
        outputs_str = ", ".join(f'"{k}"' for k in [gap.missing_key])
        description = gap.suggested_description.replace('"', '\\"')
        
        return template.format(
            skill_name=skill_name,
            class_name=class_name,
            input_contract=f"[{inputs_str}]" if inputs_str else "[]",
            output_contract=f"[{outputs_str}]",
            description=description,
            missing_key=gap.missing_key,
        )


# ─────────────────────────────────────────────
# L2SelfHealingPipeline — 自愈执行器
# ─────────────────────────────────────────────

class L2SelfHealingPipeline:
    """
    L2 自愈层主执行器。
    
    工作流程：
    1. 接收 L1 失败结果，提取失败原因
    2. 调用 GapAnalyzer 分析缺口
    3. 调用 SkillSynthesizer 生成新技能
    4. 热加载新技能
    5. 调用 L1 重试
    6. 最多重试 MAX_RETRIES 次
    """
    
    MAX_RETRIES: int = 3
    
    def __init__(
        self,
        registry: Optional[SkillRegistry] = None,
        failed_result: Optional[Any] = None,
        llm_mode: bool = False,
    ):
        self.registry = registry or SkillRegistry()
        self.failed_result = failed_result
        self.llm_mode = llm_mode
        self.gap_analyzer = GapAnalyzer(self.registry)
        self.synthesizer = SkillSynthesizer(
            self.registry, mode="llm" if llm_mode else "template"
        )
    
    def execute(self, state: State) -> PIPResult:
        """
        执行 L2 自愈流程。
        
        Returns:
            PipelineResult — 成功时返回新结果，失败时返回含 gaps 的结果
        """
        logs: List[str] = []
        all_artifacts: Dict[str, Any] = {}
        start_time = time.time()
        
        # 提取失败原因
        if self.failed_result:
            failed_skills = self.failed_result.get("failed_skills", [])
            gaps = self.failed_result.get("gaps", [])
        else:
            failed_skills = []
            gaps = []
        
        logs.append(f"🔧 L2 Self-Healing triggered (retries={self.MAX_RETRIES})")
        
        # ── Step 1: 缺口分析 ─────────────────────
        if not gaps:
            start_keys = {k for k, v in state.items() if v is not None}
            goal_keys = self.failed_result.get("target_keys", set()) if self.failed_result else set()
            
            gaps = self.gap_analyzer.analyze(
                start_keys=start_keys,
                goal_keys=goal_keys,
                failed_skills=failed_skills,
            )
        
        if not gaps:
            logs.append("✅ No gaps found — original failure was not a capability gap")
            return PIPResult(
                status="failed",
                final_state=state,
                executed_skills=[],
                failed_skills=failed_skills,
                logs=logs,
                artifacts=all_artifacts,
                duration_seconds=time.time() - start_time,
                error="L2: No capability gaps found",
            )
        
        # 输出缺口报告
        report = self.gap_analyzer.generate_report(gaps)
        for line in report.split("\n"):
            logs.append(line)
        
        # ── Step 2: 合成技能 ─────────────────────
        logs.append(f"\n🧬 Synthesizing {len(gaps)} skills...")
        synthesized = self.synthesizer.synthesize_batch(gaps)
        
        loaded = []
        for skill in synthesized:
            ok = self.synthesizer.load_synthesized(skill)
            if ok:
                loaded.append(skill.skill_name)
                logs.append(f"  ✅ {skill.skill_name}")
            else:
                logs.append(f"  ❌ {skill.skill_name}")
        
        if not loaded:
            logs.append("❌ L2: No skills could be synthesized")
            return PIPResult(
                status="failed",
                final_state=state,
                executed_skills=[],
                failed_skills=failed_skills,
                logs=logs,
                artifacts=all_artifacts,
                duration_seconds=time.time() - start_time,
                error="L2: Synthesis failed for all gaps",
                gaps=[{"missing_key": g.missing_key} for g in gaps],
            )
        
        # ── Step 3: 重试 L1 ─────────────────────
        logs.append(f"\n🔄 Retrying with {len(loaded)} new skills...")
        
        from bioskills.core.pipeline import (
            L1ContractAwarePipeline, PipelineOrchestrator,
        )
        
        orchestrator = PipelineOrchestrator(
            registry=self.registry,
            enable_self_heal=False,  # 防止无限递归
        )
        
        target_keys = self.failed_result.get("target_keys", []) if self.failed_result else []
        retry_result = orchestrator.execute(
            state,
            target_keys=target_keys,
            force_contract_aware=True,
        )
        
        logs.extend(retry_result.logs)
        all_artifacts.update(retry_result.artifacts)
        
        if retry_result.status == "success":
            logs.append("✅ L2 Self-Healing SUCCESS")
            return PIPResult(
                status="success",
                final_state=retry_result.final_state,
                executed_skills=loaded + retry_result.executed_skills,
                failed_skills=[],
                logs=logs,
                artifacts=all_artifacts,
                duration_seconds=time.time() - start_time,
            )
        else:
            logs.append(f"❌ L2 Self-Healing FAILED after retry")
            return PIPResult(
                status="failed",
                final_state=retry_result.final_state,
                executed_skills=loaded,
                failed_skills=failed_skills + retry_result.failed_skills,
                logs=logs,
                artifacts=all_artifacts,
                duration_seconds=time.time() - start_time,
                error=f"L2: Still failed after {len(loaded)} synthesized skills",
                gaps=[{"missing_key": g.missing_key} for g in gaps],
            )


# ─────────────────────────────────────────────
# 代码模板库
# ─────────────────────────────────────────────



# ─────────────────────────────────────────────
# 代码模板库
# ─────────────────────────────────────────────

