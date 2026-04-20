"""Bio-Analysis LangGraph StateGraph

完整 6 阶段假设验证流水线。
"""
from __future__ import annotations
import sys
from pathlib import Path
from typing import Optional, Tuple
from dataclasses import dataclass

# 添加项目根路径
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from langgraph.graph import StateGraph, END
from langgraph.checkpoint.memory import MemorySaver

from workflow.state import BioAnalysisState, new_state
from workflow.registry import StateRegistry
from workflow.nodes.hypothesis_parser import hypothesis_parser
from workflow.nodes.literature_miner import literature_miner
from workflow.nodes.geo_retriever import geo_retriever
from workflow.nodes.effect_size_analyst import effect_size_analyst_supervisor
from workflow.nodes.meta_analyst import meta_analyst
from workflow.nodes.report_curator import report_curator


# ── 全局 Checkpoint ──
checkpointer = MemorySaver()


# ── 条件边 ──
def should_literature_miner(state: BioAnalysisState) -> str:
    if state.get("hypothesis_parser_status") == "failed":
        return "__end__"
    return "literature_miner"

def should_geo_retriever(state: BioAnalysisState) -> str:
    if state.get("literature_miner_status") == "failed":
        return "__end__"
    return "geo_retriever"

def should_effect_size_analyst(state: BioAnalysisState) -> str:
    if state.get("geo_retriever_status") == "failed":
        return "__end__"
    return "effect_size_analyst"

def should_meta_analyst(state: BioAnalysisState) -> str:
    if state.get("effect_size_analyst_status") == "failed":
        return "__end__"
    return "meta_analyst"

def should_report_curator(state: BioAnalysisState) -> str:
    if state.get("meta_analyst_status") == "failed":
        return "__end__"
    return "report_curator"


# ── 构建 StateGraph ──
def build_graph(checkpoint: bool = True) -> "CompiledGraph":
    graph = StateGraph(BioAnalysisState)
    
    graph.add_node("hypothesis_parser", hypothesis_parser)
    graph.add_node("literature_miner", literature_miner)
    graph.add_node("geo_retriever", geo_retriever)
    graph.add_node("effect_size_analyst", effect_size_analyst_supervisor)
    graph.add_node("meta_analyst", meta_analyst)
    graph.add_node("report_curator", report_curator)
    
    # 条件边（失败→END，成功→下一Phase）
    graph.add_conditional_edges("hypothesis_parser", should_literature_miner)
    graph.add_conditional_edges("literature_miner", should_geo_retriever)
    graph.add_conditional_edges("geo_retriever", should_effect_size_analyst)
    graph.add_conditional_edges("effect_size_analyst", should_meta_analyst)
    graph.add_conditional_edges("meta_analyst", should_report_curator)
    graph.add_edge("report_curator", END)
    
    graph.set_entry_point("hypothesis_parser")
    
    if checkpoint:
        return graph.compile(checkpointer=checkpointer)
    return graph.compile()


@dataclass
class ResumePoint:
    hypothesis_id: str
    from_phase: str


_compiled_graph = None

def get_graph():
    global _compiled_graph
    if _compiled_graph is None:
        _compiled_graph = build_graph(checkpoint=True)
    return _compiled_graph


def run_pipeline(hypothesis_text: str, start_from: Optional[str] = None) -> Tuple[str, BioAnalysisState]:
    """执行完整 6 阶段流水线"""
    registry, hyp_id = StateRegistry.init(hypothesis_text)
    initial_state = registry.get()
    
    print(f"\n{'='*60}")
    print(f"Bio-Analysis Pipeline Started")
    print(f"  Hypothesis: {hypothesis_text[:60]}...")
    print(f"  ID: {hyp_id}")
    print(f"{'='*60}\n")
    
    app = get_graph()
    config = {"configurable": {"thread_id": hyp_id}}
    
    if start_from:
        current = app.get_state(config)
        state_dict = dict(current.values)
        state_dict["current_phase"] = start_from
        state_dict["{}-status".format(start_from)] = "running"
        app.update_state(config, state_dict)
        final = app.invoke(None, config)
    else:
        final = app.invoke(initial_state, config)
    
    print(f"\n{'='*60}")
    print(f"Pipeline Completed")
    print(f"  ID: {hyp_id}")
    print(f"  Score: {final.get('meta_analyst_reliability_score', 0):.4f} {final.get('meta_analyst_grade', '')}")
    print(f"  Cohen d: {final.get('meta_analyst_combined_d', 0):.3f}")
    print(f"  I2: {final.get('meta_analyst_i_squared', 0):.1f}%")
    print(f"  Report: {final.get('report_curator_report_path', 'N/A')}")
    print(f"{'='*60}\n")
    
    return hyp_id, final
