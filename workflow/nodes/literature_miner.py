"""LiteratureMiner Node"""
from __future__ import annotations
import sys, json
from pathlib import Path
from typing import Optional, Dict, Any

# 添加项目根路径
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from workflow.state import BioAnalysisState

def _build_default_gene_sets(hyp_text: str) -> dict:
    """当文献挖掘无结果时，从假设文本构建默认基因集"""
    t = hyp_text.lower()
    if "caf" in t or "fibroblast" in t:
        return {"CAF_IDENTITY": ["FAP","PDGFRA","ACTA2","COL1A1","COL3A1"]}
    if "t cell" in t or "cd8" in t or "cd4" in t:
        return {"T_CELL_EXHAUSTION": ["PDCD1","HAVCR2","LAG3","TIGIT","CTLA4"]}
    if "macrophage" in t or "m2" in t or "tam" in t:
        return {"MACROPHAGE_POLARIZATION": ["CD163","MRC1","ARG1","IL10","SPP1"]}
    if "nk" in t:
        return {"NK_CELL_ACTIVITY": ["NKG2D","KLRK1","GNLY","PRF1","GZMB"]}
    if "b cell" in t:
        return {"B_CELL": ["CD19","MS4A1","CD27","IGHA","IGHG1"]}
    if "dendritic" in t or "dc " in t:
        return {"DC_MARKERS": ["CD11C","CD141","CD1C","IRF4","IRF8"]}
    return {"CUSTOM": []}


from workflow.registry import StateRegistry

SCRIPT_DIR = PROJECT_ROOT / "scripts"

def literature_miner(state: BioAnalysisState) -> BioAnalysisState:
    """
    LiteratureMiner: 文献挖掘，生成 gene_set.json 和推荐数据集列表。
    
    委派给 scripts/phase2_literature_miner.py 执行。
    """
    hyp_id = state["hypothesis_id"]
    hyp_text = state["hypothesis_text"]
    
    new_state = dict(state)
    new_state["literature_miner_status"] = "running"
    registry = StateRegistry(hyp_id)
    registry.set(new_state, phase=None)
    
    try:
        # 调用 LiteratureMiner 脚本
        script = SCRIPT_DIR / "literature_miner_literature_miner.py"
        if script.exists():
            import subprocess
            result = subprocess.run(
                [sys.executable, str(script), hyp_text,
                 "--output", str(registry.hyp_dir)],
                capture_output=True, text=True, timeout=300,
            )
            print(f"[LiteratureMiner] stdout:\n{result.stdout[:500]}")
            if result.returncode != 0:
                raise RuntimeError(result.stderr[:500])
        
        # 读取 LiteratureMiner 输出
        gene_set = registry.read_phase2()
        if gene_set:
            new_state.update({
                "literature_miner_supporting_papers": gene_set.get("supporting_papers", []),
                "literature_miner_gene_sets": gene_set.get("gene_sets", {}),
                "literature_miner_recommended_datasets": gene_set.get("recommended_datasets", []),
                "literature_miner_status": "done",
                "literature_miner_error": None,
                "current_phase": "geo_retriever",
            })
            if "literature_miner" not in new_state["checkpoints"]:
                new_state["checkpoints"] = new_state["checkpoints"] + ["literature_miner"]
            registry.set(new_state, phase="literature_miner")
            print(f"[LiteratureMiner] ✅ Found {len(gene_set.get('recommended_datasets', []))} recommended datasets")
        else:
            # 文献挖掘没有产出 gene_set（文献少/无数据）→ 生成空骨架，继续
            print(f"[LiteratureMiner] ⚠ 未找到 gene_set.json，使用内置基因集")
            new_state.update({
                "literature_miner_supporting_papers": [],
                "literature_miner_gene_sets": _build_default_gene_sets(hyp_text),
                "literature_miner_recommended_datasets": [],
                "literature_miner_status": "done",
                "literature_miner_error": "No gene_set produced (may have insufficient literature)",
                "current_phase": "geo_retriever",
            })
            if "literature_miner" not in new_state["checkpoints"]:
                new_state["checkpoints"] = new_state["checkpoints"] + ["literature_miner"]
            registry.set(new_state, phase="literature_miner")
            return new_state

    except Exception as e:
        new_state.update({
            "literature_miner_status": "failed",
            "literature_miner_error": str(e),
        })
        new_state.setdefault("errors", []).append(f"Phase2: {e}")
        registry.set(new_state, phase=None)
        print(f"[LiteratureMiner] ❌ Failed: {e}")
    
    return new_state


# ─────────────────────────────────────────────
# GEORetriever Node — 数据获取（子 Agent）
# ─────────────────────────────────────────────
