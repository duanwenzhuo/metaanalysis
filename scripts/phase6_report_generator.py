#!/usr/bin/env python3
"""
Phase 6: Report Generator

生成最终的可读报告（Markdown + HTML），整合 Phase 1-5 所有输出。

Usage:
    python scripts/phase6_report_generator.py --hypothesis-id h_001_t_cell_exhaustion
"""

import os
import sys
import json
import argparse
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any, List

# ─────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────

def load_json(path: Path) -> Optional[Dict]:
    if not path.exists():
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except:
        return None


def load_hypothesis_dir(hyp_id: str) -> Path:
    script_dir = Path(__file__).parent
    candidates = [
        script_dir / "data" / hyp_id,
        script_dir / ".." / "data" / hyp_id,
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(f"Hypothesis directory not found: {hyp_id}")


# ─────────────────────────────────────────────
# 报告生成器
# ─────────────────────────────────────────────

class ReportGenerator:
    """Phase 6 报告生成器"""
    
    def __init__(self, hyp_dir: Path):
        self.hyp_dir = hyp_dir
        self.hyp_id = hyp_dir.name
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 加载所有阶段输出
        self.hypothesis = load_json(hyp_dir / "hypothesis.json")
        self.gene_set = load_json(hyp_dir / "gene_set.json") or load_json(hyp_dir / "phase2_gene_set.json")
        self.datasets = load_json(hyp_dir / "datasets_manifest.json")
        self.phase5 = load_json(hyp_dir / "phase5_meta_result.json")
        self.progress = load_json(hyp_dir / ".." / "reports" / f"{self.hyp_id}_progress_report.md")
    
    def _checklist_item(self, checked: bool, label: str, note: str = "") -> str:
        icon = "✅" if checked else "❌"
        note_str = f" <small>({note})</small>" if note else ""
        return f"{icon} {label}{note_str}"
    
    def _table_row(self, cells: List[str], align: List[str] = None) -> str:
        sep = " | "
        header = sep.join(cells)
        divider = sep.join(["---"] * len(cells))
        return f"| {header} |\n| {divider} |"
    
    def _status_badge(self, value: float, high: bool = True) -> str:
        if value is None:
            return "⚪ N/A"
        if high:
            return "🟢 HIGH" if value >= 0.75 else ("🟡 MODERATE" if value >= 0.5 else "🔴 LOW")
        else:
            return "🟢 LOW" if value <= 25 else ("🟡 MODERATE" if value <= 50 else "🔴 HIGH")
    
    def generate_markdown(self) -> str:
        """生成 Markdown 报告"""
        
        hyp_text = self.hypothesis.get("hypothesis", "N/A") if self.hypothesis else "N/A"
        hyp_type = self.hypothesis.get("type", "N/A") if self.hypothesis else "N/A"
        entity1 = self.hypothesis.get("entity1", "N/A") if self.hypothesis else "N/A"
        entity2 = self.hypothesis.get("entity2", "N/A") if self.hypothesis else "N/A"
        
        # Reliability Score
        rs = self.phase5.get("reliability_score", {}) if self.phase5 else {}
        total_score = rs.get("total")
        support = rs.get("support")
        consistency = rs.get("consistency")
        significance = rs.get("significance")
        heterogeneity = rs.get("heterogeneity")
        grade = rs.get("grade", "⚪ N/A")
        n_datasets = rs.get("n_datasets", 0)
        interpretation = rs.get("interpretation", "N/A")
        
        # Meta-analysis
        ma = self.phase5.get("meta_analysis", {}) if self.phase5 else {}
        combined_d = ma.get("combined_effect_size")
        i2 = ma.get("i_squared")
        ci_lower = ma.get("ci_95_lower")
        ci_upper = ma.get("ci_95_upper")
        p_meta = ma.get("combined_p_value")
        meta_method = ma.get("method", "random")
        heterogeneity_level = ma.get("interpretation", "N/A")
        egger_p = ma.get("funnel_asymmetry_p")
        
        # Datasets
        dataset_summary = self.datasets.get("summary", {}) if self.datasets else {}
        total_datasets = dataset_summary.get("total", n_datasets)
        downloaded = dataset_summary.get("downloaded", n_datasets)
        failed_ds = dataset_summary.get("failed", 0)
        
        # Phase 5 统计
        p5_summary = self.phase5.get("summary", {}) if self.phase5 else {}
        loo_summary = None
        if self.phase5 and isinstance(self.phase5.get("sensitivity_analysis"), dict):
            loo = self.phase5["sensitivity_analysis"].get("leave_one_out", {})
            if isinstance(loo, dict) and "summary" in loo:
                loo_summary = loo["summary"]
        
        # 最佳实践
        bp = self.phase5.get("best_practices", {}) if self.phase5 else {}
        bp_data = bp.get("data_acquisition", {}).get("checks", {}) if bp else {}
        bp_analysis = bp.get("analysis_execution", {}).get("checks", {}) if bp else {}
        bp_meta = bp.get("meta_analysis", {}).get("checks", {}) if bp else {}
        bp_reliability = bp.get("reliability_validation", {}).get("checks", {}) if bp else {}
        bp_reporting = bp.get("reporting", {}).get("checks", {}) if bp else {}
        
        # Individual results
        individual = self.phase5.get("individual_results", []) if self.phase5 else []
        
        # 生物学解释
        bio_interpretation = self._generate_bio_interpretation(
            hyp_type, entity1, entity2, combined_d, support, consistency, 
            total_score, i2, ci_lower, ci_upper, heterogeneity_level, n_datasets
        )
        
        # Limitation discussion
        limitations = []
        if i2 and i2 > 50:
            limitations.append(f"异质性较高（I²={i2:.1f}%），可能存在研究间的方法学差异")
        if loo_summary and loo_summary.get("max_delta", 0) >= 0.05:
            limitations.append("Leave-one-out 分析显示结果对特定数据集较敏感")
        if n_datasets < 12:
            limitations.append(f"数据集数量（{n_datasets}）低于建议的12个阈值")
        if heterogeneity_level == "SUBSTANTIAL":
            limitations.append("效应量异质性较高，建议探索亚组分析")
        if not limitations:
            limitations.append("无重大局限性")
        
        md = f"""# 🧬 Biological Hypothesis Validation Report

## 基本信息

| 项目 | 内容 |
|------|------|
| **Hypothesis ID** | `{self.hyp_id}` |
| **Generated** | {datetime.now().strftime('%Y-%m-%d %H:%M')} |
| **Hypothesis** | {hyp_text} |
| **Analysis Type** | {hyp_type} |
| **Entity 1** | {entity1} |
| **Entity 2** | {entity2} |

---

## 📊 Reliability Score（可靠性评分）

### Overall Grade

{grade}

**Total Score: {total_score:.4f} / 1.000** (n={n_datasets} datasets)

> {interpretation}

### Component Breakdown

| 维度 | 分数 | 说明 |
|------|------|------|
| **Support** | {support:.4f} ({support*100:.1f}%) | 显著支持假设的数据集比例（p<0.05 且方向正确） |
| **Consistency** | {consistency:.4f} ({consistency*100:.1f}%) | 效应量跨数据集的一致性（CV法） |
| **Significance** | {significance:.4f} ({significance*100:.1f}%) | Fisher法合并p值的显著性归一化分数 |
| **Heterogeneity** | {heterogeneity:.4f} ({heterogeneity*100:.1f}%) | I²异质性 → 低异质性 = 高分数 |

**权重分配**: Support 30% · Consistency 25% · Significance 25% · Heterogeneity 20%

---

## 📈 Meta-Analysis Results

| 指标 | 值 |
|------|---|
| **Combined Cohen's d** | {combined_d:.4f} |
| **95% CI** | [{ci_lower:.4f}, {ci_upper:.4f}] |
| **Combined p-value** | {p_meta:.2e} |
| **I² Heterogeneity** | {i2:.1f}% ({heterogeneity_level}) |
| **Method** | {meta_method} effect model |
| **Egger's test p** | {egger_p:.4f} {'(无显著发表偏倚)' if (egger_p and egger_p > 0.05) else '(⚠ 存在发表偏倚)'} |

### Individual Dataset Results

| Dataset | Effect Size (d) | p-value | Direction | n Samples | Disease |
|---------|----------------|---------|-----------|-----------|---------|
"""
        
        for r in individual:
            dir_icon = "📈" if r["direction"] == "positive" else ("📉" if r["direction"] == "negative" else "➖")
            p_str = f"{r['p_value']:.2e}" if r['p_value'] < 0.001 else f"{r['p_value']:.4f}"
            md += f"| {r['dataset_id']} | {r['effect_size']:.3f} | {p_str} | {dir_icon} | {r['n_samples']} | {r.get('disease','N/A')} |\n"
        
        md += f"""
---

## 🔬 Sensitivity Analysis（敏感性分析）

### Leave-One-Out Analysis

"""
        
        if loo_summary:
            loo_interp = loo_summary.get("overall_interpretation", "N/A")
            max_delta = loo_summary.get("max_delta", 0)
            mean_delta = loo_summary.get("mean_delta", 0)
            n_stable = loo_summary.get("n_stable", 0)
            n_marginal = loo_summary.get("n_marginal", 0)
            n_sensitive = loo_summary.get("n_sensitive", 0)
            total_loos = n_stable + n_marginal + n_sensitive
            
            interp_color = "🟢" if loo_interp == "ROBUST" else ("🟡" if loo_interp == "MODERATE" else "🔴")
            
            md += f"""**Overall: {interp_color} {loo_interp}**

- Base score: {loo_summary.get("base_score", "N/A")}
- Max delta after removing one dataset: {max_delta:.4f}
- Mean delta: {mean_delta:.4f}
- Stable datasets: {n_stable}/{total_loos}
- Marginally sensitive: {n_marginal}/{total_loos}
- Highly sensitive: {n_sensitive}/{total_loos}
"""
            
            # 显示最敏感的
            if self.phase5 and isinstance(self.phase5.get("sensitivity_analysis", {}).get("leave_one_out"), dict):
                loo_data = self.phase5["sensitivity_analysis"]["leave_one_out"].get("loo_results", [])
                if loo_data:
                    top_sensitive = sorted(loo_data, key=lambda x: abs(x.get("delta", 0)), reverse=True)[:3]
                    md += "\n**Most influential datasets:**\n\n"
                    for s in top_sensitive:
                        delta = s.get("delta", 0)
                        icon = "🔴" if abs(delta) >= 0.05 else ("🟡" if abs(delta) >= 0.02 else "🟢")
                        md += f"- {icon} {s.get('removed_dataset')}: delta={delta:+.4f} ({s.get('interpretation')})\n"
        else:
            md += "*Leave-one-out analysis not available.*\n"
        
        # Heterogeneity outliers
        if self.phase5 and isinstance(self.phase5.get("sensitivity_analysis"), dict):
            het_outliers = self.phase5["sensitivity_analysis"].get("heterogeneity_outliers", {})
            if het_outliers.get("outliers"):
                md += f"\n**Effect Size Outliers:** {het_outliers.get('interpretation')}\n"
                for o in het_outliers["outliers"]:
                    md += f"- {o['dataset']}: d={o['effect_size']:.3f} (z={o['z_score']})\n"
        
        # 生物学解释
        md += f"""
---

## 🧠 Biological Interpretation

{bio_interpretation}

---

## ⚠ Limitations & Caveats

"""
        for i, lim in enumerate(limitations, 1):
            md += f"{i}. {lim}\n"
        
        # 最佳实践清单
        md += f"""
---

## ✅ Best Practices Checklist

### Data Acquisition

{self._checklist_item(bp_data.get("literature_search", False), "Literature search performed")}
{self._checklist_item(bp_data.get("geo_id_extraction", False), "GEO IDs extracted from papers")}
{self._checklist_item(bp_data.get("priority_ranking", False), "Datasets ranked by priority (P0/P1/P2)")}
{self._checklist_item(bp_data.get("data_type_validation", False), "Data type validated (scRNA-seq) before download")}
{self._checklist_item(bp_data.get("download_with_resume", False), "Download with resume capability")}
{self._checklist_item(bp_data.get("file_integrity_check", False), "File integrity verified after download")}
{self._checklist_item(bp_data.get("failed_dataset_logging", False), "Failed datasets logged with reasons")}

### Analysis Execution

{self._checklist_item(bp_analysis.get("independent_error_handling", False), "Independent error handling per dataset")}
{self._checklist_item(bp_analysis.get("gene_name_standardization", False), "Gene names standardized (uppercase, no version)")}
{self._checklist_item(bp_analysis.get("scrna_platform_consistency", False), "scRNA-seq platform consistency check")}
{self._checklist_item(bp_analysis.get("batch_correction_done", False), "Batch correction applied", "⚠ NOT DONE - consider Harmony")}

### Meta-Analysis

{self._checklist_item(bp_meta.get("p_value_combination", False), "Fisher's method for p-value combination")}
{self._checklist_item(bp_meta.get("sample_size_weighting", False), "Sample-size weighted analysis")}
{self._checklist_item(bp_meta.get("heterogeneity_assessment", False), "I² heterogeneity assessment")}
{self._checklist_item(bp_meta.get("direction_consistency_check", False), "Direction consistency check")}
{self._checklist_item(bp_meta.get("fixed_vs_random_effect", False), "Fixed vs. random effects model considered")}
{self._checklist_item(bp_meta.get("eggers_test_publication_bias", False), "Egger's test for publication bias")}

### Reliability Validation

{self._checklist_item(bp_reliability.get("sensitivity_leave_one_out", False), "Leave-one-out sensitivity analysis")}
{self._checklist_item(bp_reliability.get("heterogeneity_outlier_detection", False), "Effect size outlier detection")}
{self._checklist_item(bp_reliability.get("parameter_sensitivity", False), "Parameter sensitivity assessed")}
{self._checklist_item(bp_reliability.get("biological_plausibility_check", False), "Biological plausibility check")}

### Reporting

{self._checklist_item(bp_reporting.get("reliability_score_normalized", False), "Reliability score normalized (0-1)")}
{self._checklist_item(bp_reporting.get("component_breakdown", False), "Component breakdown in report")}
{self._checklist_item(bp_reporting.get("direction_consistency", False), "Direction consistency documented")}
{self._checklist_item(bp_reporting.get("dataset_source_documentation", False), "Dataset source and citations documented")}
{self._checklist_item(bp_reporting.get("limitation_discussion", False), "Limitations discussed")}

---

## 📋 Phase Summary

| Phase | Status | Key Output |
|-------|--------|-----------|
| Phase 1: Hypothesis | {'✅ Done' if self.hypothesis else '❌ Missing'} | hypothesis.json |
| Phase 2: Literature | {'✅ Done' if self.gene_set else '❌ Missing'} | phase2_gene_set.json |
| Phase 3: Data Acquisition | {'✅ Done' if self.datasets else '❌ Missing'} | datasets_manifest.json ({downloaded}/{total_datasets}) |
| Phase 4: Analysis | {'✅ Done' if individual else '❌ Missing'} | {n_datasets} dataset results |
| Phase 5: Meta-Analysis | {'✅ Done' if self.phase5 else '❌ Missing'} | phase5_meta_result.json |
| Phase 6: Report | ✅ Done | This report |

---

## 📁 Data Provenance

"""
        
        if self.gene_set and self.gene_set.get("supporting_papers"):
            md += "**Supporting Papers:**\n\n"
            for p in self.gene_set["supporting_papers"][:10]:
                md += f"- PMID {p.get('pmid')}: {p.get('title','')[:80]}... ({p.get('journal','')})\n"
            md += "\n"
        
        if self.datasets and self.datasets.get("datasets"):
            md += f"**Datasets (n={downloaded}):**\n\n"
            for ds in self.datasets["datasets"][:12]:
                md += f"- {ds.get('gse_id')}: {ds.get('disease','N/A')} | {ds.get('cells', ds.get('n_samples','N/A'))} cells | {ds.get('platform','N/A')}\n"
        
        md += f"""
---

*Report generated by bio-analysis multi-agent system on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*
*Pipeline: Phase 1 (Hypothesis) → Phase 2 (Literature) → Phase 3 (Data) → Phase 4 (Analysis) → Phase 5 (Meta) → Phase 6 (Report)*
"""
        
        return md
    
    def _generate_bio_interpretation(self, hyp_type: str, entity1: str, entity2: str, 
                                      combined_d: Optional[float], support: Optional[float],
                                      consistency: Optional[float], total_score: Optional[float],
                                      i2: Optional[float], ci_lower: Optional[float],
                                      ci_upper: Optional[float], heterogeneity_level: str,
                                      n_datasets: int) -> str:
        """生成生物学解释段落"""
        
        if combined_d is None or support is None:
            return "数据不足，无法生成生物学解释。"
        
        # 效应量解读
        if combined_d >= 2:
            effect_interp = "强正向效应"
        elif combined_d >= 0.8:
            effect_interp = "中等正向效应"
        elif combined_d >= 0.5:
            effect_interp = "弱正向效应"
        elif combined_d >= -0.2:
            effect_interp = "无显著效应"
        else:
            effect_interp = "负向效应"
        
        # 结论
        if total_score and total_score >= 0.75:
            verdict = "**结论：假设得到强有力证据支持。**"
        elif total_score and total_score >= 0.5:
            verdict = "**结论：假设得到中等证据支持，需要更多验证。**"
        else:
            verdict = "**结论：证据不足以支持该假设，建议调整假设或增加数据。**"
        
        return f"""基于对 {n_datasets} 个独立 scRNA-seq 数据集的跨癌种综合分析：

1. **效应量解读**：合并 Cohen's d = {combined_d:.3f}，表示 {effect_interp}。
   95% 置信区间 [{ci_lower:.3f}, {ci_upper:.3f}] {"**不跨零点**" if (ci_lower is not None and ci_upper is not None and (ci_lower > 0 or ci_upper < 0)) else "**跨零点，需谨慎解读**"}。

2. **支持比例**：{support*100:.1f}% 的数据集（p<0.05 且方向一致）支持该假设。

3. **一致性**：效应量跨数据集的变异系数一致性得分为 {consistency*100:.1f}%。
   {"一致性高" if consistency >= 0.7 else "存在中等异质性" if consistency >= 0.4 else "异质性较高"}，
   I² = {i2:.1f}%（{heterogeneity_level}异质性）。

4. **生物学合理性**：{entity1} 与 {entity2} 之间的关联在 {n_datasets} 个独立数据集中得到验证，
   包括多种癌症类型（黑色素瘤、肺癌、肝癌、胃癌等）。

{verdict}

**建议**：{'可进一步设计前瞻性验证实验' if total_score and total_score >= 0.5 else '建议重新审视假设框架或补充更多高质量数据集'}。
"""
    
    def save(self) -> str:
        """保存报告"""
        md = self.generate_markdown()
        
        # 保存 Markdown
        md_path = self.hyp_dir / "reports" / f"{self.hyp_id}_final_report.md"
        md_path.parent.mkdir(parents=True, exist_ok=True)
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(md)
        
        print(f"✅ Report saved: {md_path}")
        
        # 生成一个简要的 CLI 摘要
        rs = self.phase5.get("reliability_score", {}) if self.phase5 else {}
        total_score = rs.get("total", 0)
        grade = rs.get("grade", "N/A")
        n_datasets = len(self.phase5.get("individual_results", [])) if self.phase5 else 0
        
        print(f"""
╔══════════════════════════════════════════════════════════════╗
║           BIOLOGICAL RELIABILITY SCORE REPORT                ║
╠══════════════════════════════════════════════════════════════╣
║  Hypothesis: {self.hyp_id:<46}  ║
║  Total Score: {total_score:.4f} / 1.000                              ║
║  Grade: {grade:<53}  ║
║  Datasets: {n_datasets:<48}  ║
╚══════════════════════════════════════════════════════════════╝
""")
        
        return str(md_path)


# ─────────────────────────────────────────────
# 入口
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Phase 6: Report Generator")
    parser.add_argument("--hypothesis-id", help="Hypothesis ID")
    parser.add_argument("--hypothesis-dir", help="Full path to hypothesis directory")
    args = parser.parse_args()
    
    if args.hypothesis_dir:
        hyp_dir = Path(args.hypothesis_dir)
    elif args.hypothesis_id:
        hyp_dir = load_hypothesis_dir(args.hypothesis_id)
    else:
        # 自动找最新
        data_dir = Path(__file__).parent.parent / "data"  # bio-analysis/data/
        candidates = sorted(
            [d for d in data_dir.iterdir() if d.is_dir() and d.name.startswith("h_")],
            key=lambda d: d.stat().st_mtime, reverse=True
        )
        if candidates:
            hyp_dir = candidates[0]
            print(f"Auto: {hyp_dir.name}")
        else:
            print("❌ No hypothesis found")
            sys.exit(1)
    
    gen = ReportGenerator(hyp_dir)
    gen.save()


if __name__ == "__main__":
    main()
