#!/usr/bin/env python3
"""
Bio-Analysis CLI 入口

用法:
    python3 -m main                    # 交互模式（推荐）
    python3 -m main "T cell exhaustion"
    python3 -m main list
    python3 -m main status <id>
    python3 -m main resume <id>
"""
from __future__ import annotations
import sys, json, argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from workflow.registry import StateRegistry
from workflow.state import BioAnalysisState
from workflow.graph import run_pipeline
from agents import SupervisorAgent, run_hypothesis
from skills import (
    HypothesisParsingSkill, LiteratureSearchSkill,
    GEODownloadSkill, DatasetAnalysisSkill,
    MetaAnalysisSkill, ReportGenerationSkill,
)

# ─────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────

def bold(text: str) -> str:
    return f"\033[1m{text}\033[0m"

def green(text: str) -> str:
    return f"\033[92m{text}\033[0m"

def red(text: str) -> str:
    return f"\033[91m{text}\033[0m"

def yellow(text: str) -> str:
    return f"\033[93m{text}\033[0m"

def dim(text: str) -> str:
    return f"\033[2m{text}\033[0m"

def phase_emoji(status: str) -> str:
    """状态 → emoji"""
    return {"done": "✅", "running": "⏳", "failed": "❌", "pending": "⏸"}.get(status, "⏸")

def grade_emoji(score: float) -> str:
    if score >= 0.75:
        return green("🟢 HIGH")
    elif score >= 0.50:
        return yellow("🟡 MODERATE")
    else:
        return red("🔴 LOW")

# ─────────────────────────────────────────────
# 命令实现
# ─────────────────────────────────────────────

def cmd_run(args):
    """运行完整假设验证流水线"""
    hyp_text = args.hypothesis
    hyp_id, final = run_pipeline(hyp_text, start_from=args.from_phase)
    print_summary(final)
    return 0


def cmd_list(args):
    """列出所有假设"""
    registry = StateRegistry("")
    all_hyps = registry.list_all()

    print()
    print(bold("  Bio-Analysis 假设库"))
    print("  " + "─" * 55)

    if not all_hyps:
        print(yellow("  暂无假设。输入 python3 -m main 启动新假设验证 →"))
        print()
        return 0

    for h in all_hyps:
        hid = h.get("hypothesis_id", "?")
        hyp_text = h.get("hypothesis_text", "?")
        if len(hyp_text) > 45:
            hyp_text = hyp_text[:42] + "..."

        cur = h.get("current_phase", "?")
        score = h.get("meta_analyst_reliability_score", None)
        score_str = f"{score:.2f}" if score is not None else "—"

        print(f"  {green(hid[:20])}")
        print(f"    {hyp_text}")
        print(f"    Phase: {cur:<8}  Score: {score_str}")
        print()

    print()
    return 0


def cmd_status(args):
    """显示假设详细状态"""
    registry = StateRegistry(args.hypothesis_id)
    state = registry.get()

    if not state:
        print(red(f"\n  假设 '{args.hypothesis_id}' 未找到。\n"))
        return 1

    hyp_text = state.get("hypothesis_text", "?")
    print()
    print(bold("  " + hyp_text[:60]))
    print(f"  ID: {args.hypothesis_id}")
    print("  " + "─" * 55)

    # HypothesisParser
    p1s = state.get("hypothesis_parser_status", "?")
    p1t = state.get("hypothesis_parser_type", "?")
    p1e = state.get("hypothesis_parser_entity1", "?")
    print(f"  {phase_emoji(p1s)} [HypothesisParser] 假设解析")
    print(f"     类型: {p1t}  |  实体: {p1e}")

    # LiteratureMiner
    p2s = state.get("literature_miner_status", "?")
    np2 = state.get("literature_miner_n_papers", 0)
    nd2 = state.get("literature_miner_n_datasets", 0)
    gs = state.get("literature_miner_gene_sets", {})
    n_genes = sum(len(v) for v in gs.values()) if gs else 0
    print(f"  {phase_emoji(p2s)} [LiteratureMiner] 文献挖掘")
    print(f"     文献: {np2} 篇  |  基因集: {n_genes} 个基因  |  数据集: {nd2} 个")

    # GEORetriever
    p3s = state.get("geo_retriever_status", "?")
    nd3 = state.get("geo_retriever_n_downloaded", 0)
    nf3 = state.get("geo_retriever_n_failed", 0)
    print(f"  {phase_emoji(p3s)} [GEORetriever] 数据下载")
    print(f"     已下载: {nd3}  |  失败: {nf3}")

    # EffectSizeAnalyst
    p4s = state.get("effect_size_analyst_status", "?")
    na4 = state.get("effect_size_analyst_n_analyzed", 0)
    nt4 = state.get("effect_size_analyst_n_total", 0)
    print(f"  {phase_emoji(p4s)} [EffectSizeAnalyst] 效应量分析")
    print(f"     已分析: {na4}/{nt4}")

    # MetaAnalyst
    p5s = state.get("meta_analyst_status", "?")
    rs5 = state.get("meta_analyst_reliability_score", 0)
    gd5 = state.get("meta_analyst_grade", "")
    cd5 = state.get("meta_analyst_combined_d", 0)
    i25 = state.get("meta_analyst_i_squared", 0)
    ci_l = state.get("meta_analyst_ci_lower", 0)
    ci_u = state.get("meta_analyst_ci_upper", 0)
    if rs5 > 0:
        print(f"  {phase_emoji(p5s)} [MetaAnalyst] 元分析")
        print(f"     Reliability Score: {rs5:.4f} {gd5}")
        print(f"     Cohen's d: {cd5:.3f} [{ci_l:.2f}, {ci_u:.2f}]")
        print(f"     I²: {i25:.1f}%")
    else:
        print(f"  {phase_emoji(p5s)} [MetaAnalyst] 元分析")

    # ReportCurator
    p6s = state.get("report_curator_status", "?")
    rp6 = state.get("report_curator_report_path", "—")
    if rp6 and rp6 != "—":
        print(f"  {phase_emoji(p6s)} [ReportCurator] 报告生成")
        print(f"     报告: {rp6}")
    else:
        print(f"  {phase_emoji(p6s)} [ReportCurator] 报告生成")

    print()
    return 0


def cmd_resume(args):
    """断点续跑"""
    registry = StateRegistry(args.hypothesis_id)
    state = registry.get()
    if not state:
        print(red(f"\n  假设 '{args.hypothesis_id}' 未找到。\n"))
        return 1

    from_phase = args.from_phase or state.get("current_phase", "hypothesis_parser")
    print(f"\n  从 {from_phase} 继续 {args.hypothesis_id}...\n")

    hyp_id, final = run_pipeline(state["hypothesis_text"], start_from=from_phase)
    print_summary(final)
    return 0


# ─────────────────────────────────────────────
# 结果摘要
# ─────────────────────────────────────────────

def print_summary(state: BioAnalysisState):
    """打印结果摘要"""
    grade = state.get("meta_analyst_grade", "")
    score = state.get("meta_analyst_reliability_score", 0)
    cd = state.get("meta_analyst_combined_d", 0)
    ci_l = state.get("meta_analyst_ci_lower", 0)
    ci_u = state.get("meta_analyst_ci_upper", 0)
    i2 = state.get("meta_analyst_i_squared", 0)
    loo = state.get("meta_analyst_loo_robust", "N/A")
    report = state.get("report_curator_report_path", "—")

    print()
    print(bold("  " + "═" * 55))
    print(bold("  Bio-Analysis 结果"))
    print(bold("  " + "═" * 55))

    if score > 0:
        print(f"  Reliability Score  {score:.4f} {grade_emoji(score)}")
        print(f"  Cohen's d          {cd:.3f}  [{ci_l:.2f}, {ci_u:.2f}]")
        print(f"  I² (异质性)        {i2:.1f}%")
        print(f"  LOO 鲁棒性         {loo}")
        print()
        print(f"  📄 报告            {report}")
    else:
        print(yellow(f"  当前无可用数据完成分析"))
        print(yellow(f"  可能原因：文献未找到足够 GEO 数据集"))
        print()
        print(f"  提示：可运行 python3 -m main status <id> 查看详情")

    print(bold("  " + "═" * 55))
    print()


# ─────────────────────────────────────────────
# 交互模式（无参数时进入）
# ─────────────────────────────────────────────

def interactive_run():
    """交互式引导：用户输入假设 → 自动运行"""
    print()
    print(bold("  ╔══════════════════════════════════════════════════╗"))
    print(bold("  ║         Bio-Analysis 假设验证系统               ║"))
    print(bold("  ║         输入你的生物学假设，即可开始分析         ║"))
    print(bold("  ╚══════════════════════════════════════════════════╝"))
    print()
    print(f"  {dim('示例：T cell exhaustion in tumor microenvironment')}")
    print(f"  {dim('        SPP1+ macrophage predicts poor survival in HCC')}")
    print(f"  {dim('        B cell infiltration correlates with immunotherapy response')}")
    print()
    print(f"  {dim('输入 quit 退出')}")

    while True:
        print()
        try:
            hyp_text = input(bold("  请输入你的假设 > ")).strip()
        except (EOFError, KeyboardInterrupt):
            print("\n\n  再见！")
            return 0

        if not hyp_text:
            continue
        if hyp_text.lower() in ("quit", "q", "exit", "退出"):
            print("\n  再见！")
            return 0

        print()
        print(f"  {green('▶')} 运行假设验证...")
        print(f"  {dim(hyp_text[:60])}")
        print()

        hyp_id, final = run_pipeline(hyp_text, start_from=None)
        print_summary(final)

        # 询问是否继续
        print(f"  {dim('继续输入新假设，或 quit 退出')}")
        print()


# ─────────────────────────────────────────────
# main
# ─────────────────────────────────────────────

def main():
    import argparse as _argparse

    # ─────────────────────────────────────────────────────────────
    # 自定义 Action：当第一个参数是假设文本（不以 - 开头且不是子命令名）
    # 时，自动路由到 run 子命令。
    # 这样 `python3 -m main "假设"` 等价于 `python3 -m main run "假设"`
    # ─────────────────────────────────────────────────────────────
    class _RouteRun(_argparse.Action):
        KNOWN_CMDS = {"run", "list", "status", "resume"}
        def __call__(self, parser, namespace, values, option_string=None):
            if values and not values[0].startswith("-") and values[0] not in self.KNOWN_CMDS:
                # 第一个值是假设文本，插入 run
                setattr(namespace, self.dest, "run")
                setattr(namespace, "_hypothesis", values[0])
            else:
                setattr(namespace, self.dest, values[0] if values else None)

    parser = _argparse.ArgumentParser(
        description="Bio-Analysis CLI",
        usage=_argparse.SUPPRESS,
        formatter_class=_argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python3 -m main                                    # 交互模式（推荐）
  python3 -m main "T cell exhaustion in TME"         # 直接运行（推荐）
  python3 -m main run "T cell exhaustion in TME"     # 同上，显式指定
  python3 -m main list                               # 查看所有假设
  python3 -m main status h_001                       # 查看假设状态
  python3 -m main resume h_001                      # 断点续跑
        """,
    )
    sub = parser.add_subparsers(dest="cmd")

    # run
    p = sub.add_parser("run", help="运行假设验证流水线")
    p.add_argument("hypothesis", nargs="?", default=None, help="假设文本")
    p.add_argument("--from", dest="from_phase", default=None, help="从某 Phase 开始")
    p.set_defaults(func=cmd_run)

    # list
    sub.add_parser("list", help="列出所有假设").set_defaults(func=cmd_list)

    # status
    p = sub.add_parser("status", help="查看假设状态")
    p.add_argument("hypothesis_id", nargs="?", default=None, help="假设 ID")
    p.set_defaults(func=cmd_status)

    # resume
    p = sub.add_parser("resume", help="断点续跑")
    p.add_argument("hypothesis_id", nargs="?", default=None, help="假设 ID")
    p.add_argument("--from", dest="from_phase", default=None, help="从某 Phase 开始")
    p.set_defaults(func=cmd_resume)

    # ─────────────────────────────────────────────────────────────
    # 预处理 sys.argv：检测"假设文本"路由
    # ─────────────────────────────────────────────────────────────
    import copy
    _orig_argv = copy.copy(sys.argv)
    _args = sys.argv[1:]
    if _args and not _args[0].startswith("-") and _args[0] not in {"run", "list", "status", "resume"}:
        # 假设文本直接输入：插入 run 子命令
        sys.argv = ["main", "run"] + _args
    else:
        sys.argv = ["main"] + _args

    args = parser.parse_args()

    # 如果 hypothesis 是通过 _RouteRun 传进来的，从隐藏属性恢复
    if hasattr(args, "_hypothesis") and getattr(args, "_hypothesis", None):
        args.hypothesis = args._hypothesis

    # 恢复原始 argv
    sys.argv = _orig_argv

    # 无参数 → 交互模式
    if args.cmd is None:
        interactive_run()
        return 0

    # 假设文本为空时
    if args.cmd == "run" and not args.hypothesis:
        print(yellow("  ⚠  请提供假设文本，或直接输入 python3 -m main 交互运行\n"))
        parser.parse_args(["run", "--help"])
        return 1

    if args.cmd == "status" and not args.hypothesis_id:
        print(yellow("  ⚠  请提供假设 ID\n"))
        parser.parse_args(["status", "--help"])
        return 1

    if args.cmd == "resume" and not args.hypothesis_id:
        print(yellow("  ⚠  请提供假设 ID\n"))
        parser.parse_args(["resume", "--help"])
        return 1

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
