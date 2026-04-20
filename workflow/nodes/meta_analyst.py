"""MetaAnalyst Node"""
from __future__ import annotations
import sys, json
from pathlib import Path
from typing import Optional, Dict, Any

# 添加项目根路径
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from workflow.state import BioAnalysisState
from workflow.registry import StateRegistry

SCRIPT_DIR = PROJECT_ROOT / "scripts"

def meta_analyst(state: BioAnalysisState) -> BioAnalysisState:
    """
    MetaAnalyst: 调用 MetaAnalyst 脚本进行 Meta 分析。

    读取 effect_size_analyst_results/，生成 meta_analyst_meta_result.json。
    """
    hyp_id = state["hypothesis_id"]
    registry = StateRegistry(hyp_id)

    new_state = dict(state)
    new_state["meta_analyst_status"] = "running"
    registry.set(new_state, phase=None)

    try:
        script = SCRIPT_DIR / "phase5_meta_orchestrator.py"

        if script.exists():
            import subprocess
            result = subprocess.run(
                [sys.executable, str(script), "--hypothesis-dir", str(registry.hyp_dir)],
                capture_output=True, text=True, timeout=300,
            )
            print(f"[MetaAnalyst] stdout:\n{result.stdout[:800]}")
            if result.returncode != 0:
                raise RuntimeError(result.stderr[:500])

        # 读取 MetaAnalyst 输出
        meta_analyst_data = registry.read_phase5()
        if not meta_analyst_data:
            # 无有效数据 → 友好提示，不崩溃
            print(f"[MetaAnalyst] ⚠ 无可用数据完成元分析")
            print(f"  可能原因：LiteratureMiner-4 未找到足够的有效数据集")
            new_state.update({
                "meta_analyst_status": "no_data",
                "meta_analyst_reliability_score": 0.0,
                "meta_analyst_grade": "⚪ N/A",
                "meta_analyst_combined_d": 0.0,
                "meta_analyst_ci_lower": 0.0,
                "meta_analyst_ci_upper": 0.0,
                "meta_analyst_i_squared": 0.0,
                "meta_analyst_loo_robust": "N/A",
                "meta_analyst_error": "No valid results found",
                "current_phase": "report_curator",
            })
            if "meta_analyst" not in new_state["checkpoints"]:
                new_state["checkpoints"] = new_state["checkpoints"] + ["meta_analyst"]
            registry.set(new_state, phase="meta_analyst")
            print(f"[MetaAnalyst] ⚠ 跳过元分析（无数据）")
            return new_state

        rs = meta_analyst_data.get("reliability_score", {})
        ma = meta_analyst_data.get("meta_analysis", {})
        p5_sum = meta_analyst_data.get("summary", {})

        loo = meta_analyst_data.get("sensitivity_analysis", {}).get("leave_one_out", {})
        if isinstance(loo, dict) and "summary" in loo:
            loo_interp = loo["summary"].get("overall_interpretation", "N/A")
        else:
            loo_interp = "N/A"

        new_state.update({
            "meta_analyst_reliability_score": rs.get("total", 0.0),
            "meta_analyst_support": rs.get("support", 0.0),
            "meta_analyst_consistency": rs.get("consistency", 0.0),
            "meta_analyst_significance": rs.get("significance", 0.0),
            "meta_analyst_heterogeneity": rs.get("heterogeneity", 0.0),
            "meta_analyst_combined_d": ma.get("combined_effect_size", 0.0),
            "meta_analyst_ci_lower": ma.get("ci_95_lower", 0.0),
            "meta_analyst_ci_upper": ma.get("ci_95_upper", 0.0),
            "meta_analyst_combined_p": ma.get("combined_p_value", 1.0),
            "meta_analyst_i_squared": ma.get("i_squared", 0.0),
            "meta_analyst_eggertest_p": ma.get("funnel_asymmetry_p", 0.0),
            "meta_analyst_loo_robust": loo_interp,
            "meta_analyst_grade": rs.get("grade", "⚪ N/A"),
            "meta_analyst_status": "done",
            "meta_analyst_error": None,
            "current_phase": "report_curator",
        })
        if "meta_analyst" not in new_state["checkpoints"]:
            new_state["checkpoints"] = new_state["checkpoints"] + ["meta_analyst"]
        registry.set(new_state, phase="meta_analyst")

        print(f"[MetaAnalyst] ✅ Score={rs.get('total', 0):.4f} | "
              f"Combined d={ma.get('combined_effect_size', 0):.3f} | "
              f"I²={ma.get('i_squared', 0):.1f}% | "
              f"LOO={loo_interp}")

    except Exception as e:
        new_state.update({
            "meta_analyst_status": "failed",
            "meta_analyst_error": str(e),
        })
        new_state.setdefault("errors", []).append(f"Phase5: {e}")
        registry.set(new_state, phase=None)
        print(f"[MetaAnalyst] ❌ Failed: {e}")

    return new_state


# ─────────────────────────────────────────────
# ReportCurator Node - 报告生成
# ─────────────────────────────────────────────
