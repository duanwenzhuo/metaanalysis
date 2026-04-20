"""HypothesisParser Node"""
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

def hypothesis_parser(state: BioAnalysisState) -> BioAnalysisState:
    """
    HypothesisParser: 解析假设文本，生成结构化 hypothesis.json。
    
    输入: hypothesis_text
    输出: hypothesis_parser_* 字段
    """
    from workflow.state import next_phase as _next_phase
    
    hyp_text = state["hypothesis_text"]
    hyp_id = state["hypothesis_id"]
    
    # ── 规则解析（不需要 LLM API）
    t = hyp_text.lower()
    
    # 类型判断
    if any(x in t for x in ["exhaust", "polariz", "differ", "state", "marker"]):
        hyp_type = "Cell State"
    elif any(x in t for x in ["interact", "signal", "ligand", "receptor", "crosstalk"]):
        hyp_type = "Cell-Cell Interaction"
    elif any(x in t for x in ["spatial", "colocal", "proximity", "niche"]):
        hyp_type = "Spatial"
    else:
        hyp_type = "Cell State"
    
    # 实体1
    entity1_map = [
        (("t cell", "cd8", "cd4"), "T cells"),
        (("macrophage", "tam", "m2"), "Macrophages"),
        (("nk cell", "nk "), "NK cells"),
        (("fibroblast", "caf"), "Cancer-Associated Fibroblasts"),
    ]
    entity1 = next((v for kws, v in entity1_map if any(k in t for k in kws)), "Immune cells")
    
    # 实体2
    entity2 = ("Tumor Microenvironment" if any(x in t for x in ["tumor", "tme", "microenviron"]) 
               else "Normal Tissue" if "normal" in t
               else "Tumor Microenvironment")
    
    # 关系
    if "exhaust" in t:
        relation, target_genes = "increased exhaustion", ["PDCD1","HAVCR2","LAG3","TIGIT","CTLA4"]
    elif "m2" in t or "polariz" in t:
        relation, target_genes = "M2 polarization", ["CD163","MRC1","ARG1","IL10"]
    elif "spp1" in t:
        relation, target_genes = "SPP1 signaling", ["SPP1","CD44","ITGB1"]
    elif "caf" in t:
        relation, target_genes = "CAF activation", ["FAP","PDGFRA","ACTA2","COL1A1"]
    else:
        relation, target_genes = "changed in tumor microenvironment", []
    
    strategy = ("GSVA / AUCell scoring" if hyp_type == "Cell State" 
                else "CellChat / NicheNet" if hyp_type == "Cell-Cell Interaction"
                else "SpaGCN / Giotto")
    
    # 更新状态
    new_state = dict(state)
    new_state.update({
        "hypothesis_parser_type": hyp_type,
        "hypothesis_parser_entity1": entity1,
        "hypothesis_parser_entity2": entity2,
        "hypothesis_parser_relation": relation,
        "hypothesis_parser_strategy": strategy,
        "hypothesis_parser_target_genes": target_genes,
        "hypothesis_parser_status": "done",
        "hypothesis_parser_error": None,
        "current_phase": "literature_miner",
    })
    if "checkpoints" not in new_state:
        new_state["checkpoints"] = []
    if "hypothesis_parser" not in new_state["checkpoints"]:
        new_state["checkpoints"] = new_state["checkpoints"] + ["hypothesis_parser"]
    
    # 持久化到 Registry
    registry = StateRegistry(hyp_id)
    registry.set(new_state, phase="hypothesis_parser")
    
    # 同时写 hypothesis.json（旧脚本兼容）
    hyp_path = registry.hyp_dir / "hypothesis.json"
    with open(hyp_path, "w", encoding="utf-8") as f:
        json.dump({
            "id": hyp_id,
            "hypothesis": hyp_text,
            "type": hyp_type,
            "entity1": entity1,
            "entity2": entity2,
            "relation": relation,
            "analysis_strategy": strategy,
            "target_gene_set": target_genes,
            "minimum_datasets": 12,
            "effect_size_metric": "Cohen's d",
        }, f, ensure_ascii=False, indent=2)
    
    print(f"[HypothesisParser] ✅ Parsed: {hyp_type} | {entity1} | {relation}")
    return new_state


# ─────────────────────────────────────────────
# LiteratureMiner Node — 文献挖掘（子 Agent）
# ─────────────────────────────────────────────
