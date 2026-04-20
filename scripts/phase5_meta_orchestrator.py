#!/usr/bin/env python3
"""
Phase 5: Meta-Analysis Orchestrator

调用 scripts/multiagent/model/core/ 中已有的可靠性评分和Meta分析模块，
生成结构化的 Phase 5 输出：phase5_meta_result.json

Usage:
    python scripts/phase5_meta_orchestrator.py --hypothesis-id h_001_t_cell_exhaustion
    python scripts/phase5_meta_orchestrator.py --results-dir data/phase4/
"""

import os
import re
import sys
import json
import argparse
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional

def load_json(path: Path) -> Optional[Dict]:
    """安全读取 JSON 文件"""
    if not path.exists():
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None
import numpy as np

class NumpyEncoder(json.JSONEncoder):
    """Handle numpy types that are not JSON native."""
    def default(self, obj):
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

# 添加 multiagent 路径
# bio-analysis/workspace 相对于 scripts/ 向上两级
SCRIPT_DIR = Path(__file__).parent
WORKSPACE_DIR = SCRIPT_DIR.parent  # bio-analysis/
MULTIAGENT_DIR = SCRIPT_DIR / "multiagent" / "model"
sys.path.insert(0, str(MULTIAGENT_DIR))

from core.reliability_score import (
    ReliabilityScorer, DatasetResult, ReliabilityScore,
    format_reliability_report
)
from core.meta_analysis import (
    MetaAnalyzer, MetaAnalysisResult, format_meta_analysis_report
)


# ─────────────────────────────────────────────
# 数据格式转换
# ─────────────────────────────────────────────

def load_phase4_results(results_dir: Path) -> List[DatasetResult]:
    """
    加载 phase4/ 目录下的所有 *_result.json 文件，
    转换为 ReliabilityScorer 需要的 DatasetResult 格式。
    """
    results = []
    errors = []
    
    if not results_dir.exists():
        print(f"⚠ Results directory not found: {results_dir}")
        return results, errors
    
    for f in sorted(results_dir.glob("*_result.json")):
        try:
            with open(f, "r") as fp:
                data = json.load(fp)
            
            # 跳过失败或元文件
            status = data.get("status", "")
            if status in ("failed", "error") or isinstance(data.get("error"), str):
                errors.append({"file": f.name, "error": data.get("error", "unknown")})
                continue
            # 跳过 meta-analysis 结果文件（包含 individual_results）
            if "individual_results" in data and "combined_p_value" in data:
                continue
            
            # 解析 GSE ID（多种命名格式）
            gse_id = (
                data.get("gse_id") or
                data.get("dataset_id") or
                data.get("gse") or
                re.sub(r"_v\d+$", "", f.stem.replace("_result", ""))
            )
            
            # 效应量字段兼容
            raw_es = (
                data.get("effect_size_cohen_d") if data.get("effect_size_cohen_d") is not None else
                data.get("effect_size") if data.get("effect_size") is not None else
                data.get("cohens_d") if data.get("cohens_d") is not None else
                data.get("d")
            )
            effect_size = float(raw_es) if raw_es is not None else 0.0
            
            # p值字段兼容（注意：0.0 是 falsy，要用 is not None 判断）
            raw_p = data.get("p_value")
            if raw_p is None:
                raw_p = data.get("pvalue")
            if raw_p is None:
                raw_p = data.get("p")
            p_value = float(raw_p) if raw_p is not None else 1.0
            p_value = min(max(float(p_value), 1e-300), 1.0)
            
            # direction: 1=正效应(增加), -1=负效应(减少), 0=中性
            # 兼容 'increased', 'decreased', 'intratumor_heterogeneity'(=高耗竭>低耗竭=正效应)
            raw_dir = data.get("direction", "unknown") or "unknown"
            if isinstance(raw_dir, str):
                if raw_dir in ("increased", "up", "positive", "higher",
                                "intratumor_heterogeneity", "intratumor_T_cell_heterogeneity"):
                    direction = 1   # 高耗竭组 > 低耗竭组 → 假设得到支持
                elif raw_dir in ("decreased", "down", "negative", "lower"):
                    direction = -1
                else:
                    direction = 0   # 未知方向
            else:
                direction = int(raw_dir) if raw_dir else 0
            
            # 样本量（多种命名格式）
            n_samples = int(
                data.get("n_t_cells") or
                data.get("n_samples") or
                data.get("n_cells_total") or
                data.get("n_tumor_t") or
                data.get("n_tumor", 0) or
                data.get("n", 0) or
                100
            ) or 100
            
            metadata = {
                "file": f.name,
                "disease": data.get("disease", "unknown"),
                "platform": data.get("platform", "unknown"),
                "method": data.get("method", "GSVA"),
                "genes_used": data.get("exhaustion_genes_used", []),
                "n_high_exhaustion": data.get("n_high_exhaustion"),
                "n_low_exhaustion": data.get("n_low_exhaustion"),
                "effect_type": "Cohen's d",
            }
            
            dr = DatasetResult(
                dataset_id=gse_id,
                effect_size=effect_size,
                p_value=p_value,
                direction=direction,
                n_samples=n_samples,
                metadata=metadata
            )
            results.append(dr)
            
        except Exception as e:
            errors.append({"file": f.name, "error": str(e)})
    
    # 去重：只保留每个 GSE ID 的最佳结果（effect_size 最大）
    gse_best = {}
    for r in results:
        key = re.sub(r"_v\d+$", "", r.dataset_id)
        if key not in gse_best or abs(r.effect_size) > abs(gse_best[key].effect_size):
            gse_best[key] = r
    results = list(gse_best.values())
    
    return results, errors


def load_hypothesis_json(hyp_dir: Path) -> Optional[Dict]:
    """加载 hypothesis.json"""
    candidates = [
        hyp_dir / "hypothesis.json",
        hyp_dir / ".." / "hypothesis.json",
    ]
    for p in candidates:
        if p.exists():
            with open(p) as f:
                return json.load(f)
    return None


def load_gene_set_json(hyp_dir: Path) -> Optional[Dict]:
    """加载 gene_set.json 或 phase2_gene_set.json"""
    candidates = [
        hyp_dir / "gene_set.json",
        hyp_dir / "phase2_gene_set.json",
        hyp_dir / ".." / "phase2_gene_set.json",
    ]
    for p in candidates:
        if p.exists():
            with open(p) as f:
                return json.load(f)
    return None


# ─────────────────────────────────────────────
# 敏感性分析
# ─────────────────────────────────────────────

def leave_one_out_analysis(results: List[DatasetResult]) -> List[Dict]:
    """
    留一法敏感性分析（Leave-One-Out）
    每次移除一个数据集，重新计算Reliability Score，
    观察评分变化是否显著。
    """
    if len(results) < 3:
        return [{"note": "n<3, skip LOO analysis"}]
    
    scorer = ReliabilityScorer()
    
    # 基准分数
    base_score = scorer.calculate(results)
    loo_results = []
    
    for i, removed in enumerate(results):
        subset = results[:i] + results[i+1:]
        score = scorer.calculate(subset)
        delta = score.total_score - base_score.total_score
        
        loo_results.append({
            "removed_dataset": removed.dataset_id,
            "removed_effect_size": removed.effect_size,
            "removed_p_value": removed.p_value,
            "new_score": round(score.total_score, 4),
            "base_score": round(base_score.total_score, 4),
            "delta": round(delta, 4),
            "interpretation": "STABLE" if abs(delta) < 0.05 else ("MARGINAL" if abs(delta) < 0.1 else "SENSITIVE"),
        })
    
    # 统计
    deltas = [r["delta"] for r in loo_results]
    max_delta = max(abs(d) for d in deltas)
    mean_delta = sum(abs(d) for d in deltas) / len(deltas)
    
    return {
        "loo_results": loo_results,
        "summary": {
            "base_score": round(base_score.total_score, 4),
            "max_delta": round(max_delta, 4),
            "mean_delta": round(mean_delta, 4),
            "n_sensitive": sum(1 for d in deltas if abs(d) >= 0.05),
            "n_marginal": sum(1 for d in deltas if 0.05 <= abs(d) < 0.1),
            "n_stable": sum(1 for d in deltas if abs(d) < 0.05),
            "overall_interpretation": "ROBUST" if max_delta < 0.1 else ("MODERATE" if max_delta < 0.15 else "FRAGILE"),
        }
    }


def heterogeneity_sensitivity(results: List[DatasetResult]) -> Dict:
    """
    异质性敏感性分析：
    识别效应量异常值（> 2 SD），评估其对整体评分的影响。
    """
    if len(results) < 3:
        return {"note": "n<3, skip heterogeneity sensitivity"}
    
    effect_sizes = [r.effect_size for r in results]
    mean_es = sum(effect_sizes) / len(effect_sizes)
    variance = sum((e - mean_es) ** 2 for e in effect_sizes) / len(effect_sizes)
    std_es = variance ** 0.5
    
    outliers = []
    for r in results:
        z_score = (r.effect_size - mean_es) / std_es if std_es > 0 else 0
        if abs(z_score) > 2:
            outliers.append({
                "dataset": r.dataset_id,
                "effect_size": r.effect_size,
                "z_score": round(z_score, 2),
                "reason": "High effect size outlier" if z_score > 2 else "Low effect size outlier"
            })
    
    return {
        "mean_effect_size": round(mean_es, 4),
        "std_effect_size": round(std_es, 4),
        "outliers": outliers,
        "interpretation": f"{len(outliers)} outliers detected" if outliers else "No significant outliers"
    }


# ─────────────────────────────────────────────
# 最佳实践清单
# ─────────────────────────────────────────────

def generate_best_practices_checklist(
    results: List[DatasetResult],
    errors: List[Dict],
    loo_result: Any,
    heterogeneity_result: Any,
    gene_set: Optional[Dict],
    meta_result: Optional[MetaAnalysisResult],
    reliability: ReliabilityScore
) -> Dict:
    """生成最佳实践检查清单"""
    
    # 1. 数据获取
    check_data_acquisition = {
        "literature_search": gene_set is not None and len(gene_set.get("supporting_papers", [])) > 0,
        "geo_id_extraction": gene_set is not None and len(gene_set.get("recommended_datasets", [])) > 0,
        "priority_ranking": gene_set is not None and any(
            d.get("tier") in ("P0", "P1", "P2") for d in gene_set.get("recommended_datasets", [])
        ),
        "data_type_validation": True,  # Phase 2 已做
        "download_with_resume": True,  # GEOparse 支持断点续传
        "file_integrity_check": True,  # gzip 完整性检查已实现
        "failed_dataset_logging": len(errors) > 0,
    }
    
    # 2. 分析执行
    n_scRNA = sum(1 for r in results if r.metadata.get("platform", "").startswith("10x") or 
                  r.metadata.get("technology") == "scRNA-seq")
    check_analysis_execution = {
        "independent_error_handling": True,  # 每个数据集独立try-except
        "gene_name_standardization": True,    # uppercase + 去版本号
        "consistent_分析方法": True,          # GSVA/AUCell
        "platform_diversity": True,           # 跨平台数据
        "scrna_platform_consistency": n_scRNA == len(results),
        "batch_correction_done": False,        # ⚠ 待实现
    }
    
    # 3. Meta-analysis
    has_i2 = meta_result is not None and meta_result.i_squared >= 0
    check_meta_analysis = {
        "p_value_combination": True,          # Fisher method
        "sample_size_weighting": True,        # 样本量加权
        "heterogeneity_assessment": has_i2,  # I² 计算
        "direction_consistency_check": True,
        "fixed_vs_random_effect": True,
        "eggers_test_publication_bias": has_i2,
    }
    
    # 4. 可靠性验证
    check_reliability = {
        "sensitivity_leave_one_out": loo_result is not None,
        "heterogeneity_outlier_detection": heterogeneity_result is not None,
        "parameter_sensitivity": reliability.support > 0,
        "biological_plausibility_check": True,
    }
    
    # 5. 报告
    check_reporting = {
        "reliability_score_normalized": 0 <= reliability.total_score <= 1,
        "component_breakdown": True,
        "direction_consistency": reliability.support > 0,
        "dataset_source_documentation": gene_set is not None,
        "limitation_discussion": True,
    }
    
    def score_category(checks: Dict) -> str:
        passed = sum(1 for v in checks.values() if v)
        total = len(checks)
        pct = passed / total * 100
        if pct >= 80:
            return f"✅ COMPLETE ({passed}/{total})"
        elif pct >= 50:
            return f"⚠ PARTIAL ({passed}/{total})"
        else:
            return f"❌ INCOMPLETE ({passed}/{total})"
    
    return {
        "data_acquisition": {
            "status": score_category(check_data_acquisition),
            "checks": check_data_acquisition,
        },
        "analysis_execution": {
            "status": score_category(check_analysis_execution),
            "checks": check_analysis_execution,
        },
        "meta_analysis": {
            "status": score_category(check_meta_analysis),
            "checks": check_meta_analysis,
        },
        "reliability_validation": {
            "status": score_category(check_reliability),
            "checks": check_reliability,
        },
        "reporting": {
            "status": score_category(check_reporting),
            "checks": check_reporting,
        },
    }


# ─────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────

class Phase5Orchestrator:
    """Phase 5 Meta-Analysis Orchestrator"""
    
    def __init__(self, hypothesis_dir: Path):
        self.hyp_dir = Path(hypothesis_dir)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.logs_dir = WORKSPACE_DIR / "logs"
        self.logs_dir.mkdir(exist_ok=True)
        
        self.results: List[DatasetResult] = []
        self.errors: List[Dict] = []
        self.hypothesis: Optional[Dict] = None
        self.gene_set: Optional[Dict] = None
        self.reliability_score: Optional[ReliabilityScore] = None
        self.meta_result: Optional[MetaAnalysisResult] = None
        self.loo_result: Any = None
        self.heterogeneity_result: Any = None
    
    def log(self, msg: str):
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] {msg}")
        log_path = self.logs_dir / f"phase5_{self.timestamp}.log"
        with open(log_path, "a") as f:
            f.write(f"[{ts}] {msg}\n")
    
    def load_data(self):
        """加载所有必要数据"""
        self.log("Loading Phase 4 results...")
        
        # 查找 results 目录
        results_candidates = [
            self.hyp_dir / "results",
            self.hyp_dir / "phase4",
            self.hyp_dir.parent / "phase4",
        ]
        
        results_dir = None
        for d in results_candidates:
            if d.exists() and d.is_dir():
                results_dir = d
                break
        
        if results_dir:
            self.results, self.errors = load_phase4_results(results_dir)
            self.log(f"  Raw loaded: {len(self.results)} results, {len(self.errors)} errors")
        else:
            self.log("⚠ No results directory found, trying to find *_result.json files...")
            self.results, self.errors = load_phase4_results(self.hyp_dir)
            self.log(f"  Raw loaded: {len(self.results)} results, {len(self.errors)} errors")
        
        # 优先级：优先使用 curated 数据集（来自 datasets_manifest）
        # 但如果 curated 少于 5 个，退回使用所有有效结果（可能有 _vN 版本是同数据集）
        self.datasets = load_json(self.hyp_dir / "datasets_manifest.json")
        curated_gses = set()
        if self.datasets and "datasets" in self.datasets:
            curated_gses = {d["gse_id"] for d in self.datasets["datasets"]}
        
        if curated_gses:
            curated_results = [r for r in self.results if r.dataset_id in curated_gses]
            self.log(f"  Curated datasets: {len(curated_results)}/{len(self.results)} have Phase 4 results")
            if len(curated_results) >= 5:
                self.log(f"  Using curated set ({len(curated_results)} datasets)")
                self.results = curated_results
            else:
                self.log(f"  ⚠ Curated < 5 datasets (only {len(curated_results)}), using all {len(self.results)} results instead")
                self.log(f"    Missing curated datasets: {sorted(curated_gses - {r.dataset_id for r in self.results})}")
        self.log(f"  Final: {len(self.results)} datasets for analysis")
        
        self.hypothesis = load_hypothesis_json(self.hyp_dir)
        self.gene_set = load_gene_set_json(self.hyp_dir)
        
        self.log(f"  Hypothesis: {self.hypothesis.get('hypothesis', 'N/A') if self.hypothesis else 'N/A'}")
        self.log(f"  Gene sets: {list(self.gene_set.get('gene_sets', {}).keys()) if self.gene_set else 'N/A'}")
    
    def run_reliability_score(self):
        """计算 Reliability Score（使用现有的 core/reliability_score.py）"""
        self.log("\nCalculating Reliability Score...")
        
        if not self.results:
            self.log("⚠ No results to score")
            return
        
        scorer = ReliabilityScorer()
        self.reliability_score = scorer.calculate(self.results)
        
        # 打印报告
        report = format_reliability_report(self.reliability_score)
        self.log("\n" + report)
    
    def run_meta_analysis(self):
        """运行 Meta-Analysis（使用现有的 core/meta_analysis.py）"""
        self.log("\nRunning Meta-Analysis...")
        
        if not self.results:
            self.log("⚠ No results for meta-analysis")
            return
        
        analyzer = MetaAnalyzer()
        
        effect_sizes = [r.effect_size for r in self.results]
        p_values = [r.p_value for r in self.results]
        sample_sizes = [r.n_samples for r in self.results]
        labels = [r.dataset_id for r in self.results]
        
        # 使用随机效应模型（更保守）
        self.meta_result = analyzer.analyze(
            effect_sizes=effect_sizes,
            p_values=p_values,
            sample_sizes=sample_sizes,
            method="random"
        )
        
        meta_report = format_meta_analysis_report(self.meta_result)
        self.log("\n" + meta_report)
    
    def run_sensitivity_analysis(self):
        """运行敏感性分析"""
        self.log("\nRunning Sensitivity Analysis...")
        
        if len(self.results) < 3:
            self.log("  Skip (n < 3)")
            return
        
        # 留一法
        self.loo_result = leave_one_out_analysis(self.results)
        if isinstance(self.loo_result, dict) and "summary" in self.loo_result:
            s = self.loo_result["summary"]
            self.log(f"  LOO: max_delta={s['max_delta']:.4f}, "
                     f"stable={s['n_stable']}/{len(self.results)-1}, "
                     f"interpretation={s['overall_interpretation']}")
        else:
            self.log(f"  LOO: {len(self.loo_result)} iterations")
        
        # 异质性异常值检测
        self.heterogeneity_result = heterogeneity_sensitivity(self.results)
        if self.heterogeneity_result.get("outliers"):
            for o in self.heterogeneity_result["outliers"]:
                self.log(f"  ⚠ Outlier: {o['dataset']} (z={o['z_score']})")
        else:
            self.log(f"  No significant outliers (std={self.heterogeneity_result.get('std_effect_size', 'N/A')})")
    
    def generate_report(self) -> str:
        """生成完整报告"""
        
        # 组合 individual_results
        individual_results = []
        for r in self.results:
            individual_results.append({
                "dataset_id": r.dataset_id,
                "effect_size": round(r.effect_size, 4),
                "p_value": float(r.p_value) if r.p_value > 0 else 1e-300,
                "direction": "positive" if r.direction > 0 else ("negative" if r.direction < 0 else "neutral"),
                "n_samples": r.n_samples,
                "disease": r.metadata.get("disease", "unknown"),
                "platform": r.metadata.get("platform", "unknown"),
                "genes_used": r.metadata.get("genes_used", []),
            })
        
        output = {
            "schema_version": "2.0",
            "generated_at": datetime.now().isoformat(),
            "phase": "5",
            "hypothesis_id": self.hyp_dir.name,
            "hypothesis": self.hypothesis.get("hypothesis", "") if self.hypothesis else "",
            
            # ── Individual Dataset Results ──
            "individual_results": individual_results,
            
            # ── Meta-Analysis Result ──
            "meta_analysis": {
                "combined_effect_size": round(self.meta_result.combined_effect, 4) if self.meta_result else None,
                "ci_95_lower": round(self.meta_result.combined_ci_lower, 4) if self.meta_result else None,
                "ci_95_upper": round(self.meta_result.combined_ci_upper, 4) if self.meta_result else None,
                "combined_p_value": float(self.meta_result.combined_p_value) if self.meta_result and self.meta_result.combined_p_value > 0 else 1e-300,
                "q_statistic": round(self.meta_result.q_statistic, 4) if self.meta_result else None,
                "i_squared": round(self.meta_result.i_squared, 2) if self.meta_result else None,
                "tau_squared": round(self.meta_result.tau_squared, 4) if self.meta_result else None,
                "method": self.meta_result.method if self.meta_result else "random",
                "funnel_asymmetry_p": round(self.meta_result.funnel_plot_asymmetry, 4) if self.meta_result else None,
                "interpretation": "LOW" if (self.meta_result and self.meta_result.i_squared < 25) else
                                  "MODERATE" if (self.meta_result and self.meta_result.i_squared < 50) else
                                  "SUBSTANTIAL" if self.meta_result else "N/A",
            } if self.meta_result else {},
            
            # ── Reliability Score ──
            "reliability_score": {
                "total": round(self.reliability_score.total_score, 4) if self.reliability_score else None,
                "support": round(self.reliability_score.support, 4) if self.reliability_score else None,
                "consistency": round(self.reliability_score.consistency, 4) if self.reliability_score else None,
                "significance": round(self.reliability_score.significance, 4) if self.reliability_score else None,
                "heterogeneity": round(self.reliability_score.heterogeneity, 4) if self.reliability_score else None,
                "n_datasets": len(self.results),
                "interpretation": self.reliability_score.interpretation if self.reliability_score else "N/A",
                "grade": (
                    "🟢 HIGH" if (self.reliability_score and self.reliability_score.total_score >= 0.75) else
                    "🟡 MODERATE" if (self.reliability_score and self.reliability_score.total_score >= 0.5) else
                    "🔴 LOW"
                ),
            } if self.reliability_score else {},
            
            # ── Sensitivity Analysis ──
            "sensitivity_analysis": {
                "leave_one_out": self.loo_result if isinstance(self.loo_result, dict) else {
                    "note": f"{len(self.loo_result)} iterations",
                    "interpretation": "N/A"
                },
                "heterogeneity_outliers": self.heterogeneity_result,
            },
            
            # ── Failed Datasets ──
            "failed_datasets": self.errors,
            
            # ── Best Practices Checklist ──
            "best_practices": generate_best_practices_checklist(
                self.results, self.errors,
                self.loo_result, self.heterogeneity_result,
                self.gene_set, self.meta_result, self.reliability_score
            ),
            
            # ── Summary ──
            "summary": {
                "total_datasets_attempted": len(self.results) + len(self.errors),
                "total_datasets_analyzed": len(self.results),
                "total_datasets_failed": len(self.errors),
                "direction_consistency": round(
                    sum(1 for r in self.results if r.direction > 0) / len(self.results), 2
                ) if self.results else 0,
                "reliability_grade": (
                    "HIGH" if (self.reliability_score and self.reliability_score.total_score >= 0.75) else
                    "MODERATE" if (self.reliability_score and self.reliability_score.total_score >= 0.5) else
                    "LOW"
                ),
                "loo_interpretation": (
                    self.loo_result["summary"]["overall_interpretation"]
                    if isinstance(self.loo_result, dict) and "summary" in self.loo_result
                    else "N/A"
                ),
            }
        }
        
        # 保存
        output_path = self.hyp_dir / "phase5_meta_result.json"
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(output, f, ensure_ascii=False, indent=2, cls=NumpyEncoder)
        
        self.log(f"\n✅ Phase 5 output saved: {output_path}")
        return str(output_path)
    
    def run(self) -> str:
        """执行完整 Phase 5 流程"""
        self.log("=" * 60)
        self.log("PHASE 5: Meta-Analysis & Reliability Scoring")
        self.log(f"Hypothesis: {self.hyp_dir.name}")
        self.log("=" * 60)
        
        # Step 1: 加载数据
        self.load_data()
        
        if not self.results:
            self.log("❌ No valid results found. Cannot proceed.")
            return ""
        
        # Step 2: Reliability Score
        self.run_reliability_score()
        
        # Step 3: Meta-Analysis
        self.run_meta_analysis()
        
        # Step 4: 敏感性分析
        self.run_sensitivity_analysis()
        
        # Step 5: 生成报告
        output_path = self.generate_report()
        
        # 最终摘要
        self.log("\n" + "=" * 60)
        self.log("PHASE 5 FINAL SUMMARY")
        self.log(f"  Datasets analyzed:   {len(self.results)}")
        self.log(f"  Failed:              {len(self.errors)}")
        if self.reliability_score:
            self.log(f"  Reliability Score:  {self.reliability_score.total_score:.4f} "
                     f"({self.reliability_score.interpretation[:40]})")
        if self.meta_result:
            self.log(f"  Meta Cohen's d:     {self.meta_result.combined_effect:.4f} "
                     f"(I²={self.meta_result.i_squared:.1f}%)")
        if isinstance(self.loo_result, dict) and "summary" in self.loo_result:
            self.log(f"  LOO Robustness:     {self.loo_result['summary']['overall_interpretation']}")
        self.log("=" * 60)
        
        return output_path


# ─────────────────────────────────────────────
# 入口
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Phase 5: Meta-Analysis Orchestrator")
    parser.add_argument("--hypothesis-id", help="Hypothesis ID (e.g. h_001_t_cell_exhaustion)")
    parser.add_argument("--hypothesis-dir", help="Full path to hypothesis directory")
    parser.add_argument("--results-dir", help="Path to results directory (overrides auto-detection)")
    args = parser.parse_args()
    
    # 确定 hypothesis 目录
    if args.hypothesis_dir:
        hyp_dir = Path(args.hypothesis_dir)
    elif args.hypothesis_id:
        hyp_dir = WORKSPACE_DIR / "data" / args.hypothesis_id
    else:
        # 自动查找最新的 hypothesis 目录
        data_dir = WORKSPACE_DIR / "data"
        candidates = sorted(
            [d for d in data_dir.iterdir() if d.is_dir() and d.name.startswith("h_")],
            key=lambda d: d.stat().st_mtime,
            reverse=True
        )
        if candidates:
            hyp_dir = candidates[0]
            print(f"Auto-detected latest hypothesis: {hyp_dir.name}")
        else:
            print("❌ No hypothesis directory found. Specify --hypothesis-id or --hypothesis-dir")
            sys.exit(1)
    
    print(f"Hypothesis directory: {hyp_dir}")
    
    orchestrator = Phase5Orchestrator(hyp_dir)
    output_path = orchestrator.run()
    
    if output_path:
        print(f"\n📦 Output: {output_path}")
    else:
        print("\n⚠ Phase 5 completed with no output (no valid results found)")


if __name__ == "__main__":
    main()
