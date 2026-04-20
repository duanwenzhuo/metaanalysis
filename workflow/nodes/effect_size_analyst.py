"""EffectSizeAnalyst Node — 已迁移到 bioskills

改动：
1. 不再调用硬编码的 phase4_single_dataset.py
2. 改为调用 SkillPipelineEngine.run()
3. 动态基因集从 GeneSetDB 解析（替换 T cell exhaustion 硬编码）
"""
from __future__ import annotations
import sys, json
from pathlib import Path
from typing import Optional, Dict, Any
from concurrent.futures import ThreadPoolExecutor, as_completed

# 添加项目根路径
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from workflow.state import BioAnalysisState
from workflow.registry import StateRegistry


def _analyze_single_dataset_bioskills(gse_id: str, hyp_id: str, gene_sets: dict, cell_type: str, manifest_paths: dict = None) -> dict:
    """
    单个数据集分析 — 使用 bioskills 管线。
    
    这是 EffectSizeAnalyst 的 Worker 逻辑。
    返回 DatasetResult dict。
    """
    from bioskills.execution.pipeline_engine import SkillPipelineEngine
    
    try:
        # 读取已有的结果（断点续跑）
        registry = StateRegistry(hyp_id)
        existing = registry.hyp_dir / "results" / f"{gse_id}_result.json"
        if existing.exists():
            result = json.loads(open(existing).read())
            if result.get("status") != "failed":
                print(f"[Phase4 Worker {gse_id}] ⏭ Already done, skipping")
                return result
        
        # 构建 bioskills State
        # 从 manifest 读取实际路径，或回退到 data/geo/{gse_id}
        gse_path = manifest_paths.get(gse_id) or str(PROJECT_ROOT / "data" / "geo" / gse_id)
        
        # 找数据文件：支持 10X 目录、matrix 文件、H5AD
        input_path = None
        input_dir = gse_path

        # 优先找 10X filtered_feature_bc_matrix 目录
        input_dir_path = Path(input_dir)
        if input_dir_path.is_dir():
            for sub in input_dir_path.iterdir():
                if sub.is_dir() and ("filtered" in sub.name.lower() or "matrix" in sub.name.lower()):
                    # 检查里面是否有 matrix.mtx 或 features.tsv
                    if any(f.name in ("matrix.mtx.gz", "matrix.mtx", "features.tsv.gz", "features.tsv") for f in sub.iterdir()):
                        input_path = str(sub)
                        break

        # 再按扩展名找
        if not input_path and input_dir_path.is_dir():
            for ext in [".h5ad", ".h5", ".mtx", ".mtx.gz", ".csv", ".tsv"]:
                candidates = list(input_dir_path.glob(f"*{ext}"))
                if candidates:
                    input_path = str(candidates[0])
                    break

        # 再尝试 data/geo/{gse_id} 作为回退
        if not input_path:
            fallback = PROJECT_ROOT / "data" / "geo" / gse_id
            if fallback.is_dir():
                for ext in [".h5ad", ".h5", ".mtx"]:
                    cands = list(fallback.glob(f"*{ext}"))
                    if cands:
                        input_path = str(cands[0]); break
                if not input_path:
                    for sub in fallback.iterdir():
                        if sub.is_dir() and "filtered" in sub.name.lower():
                            input_path = str(sub); break
        
        if not input_path:
            return {
                "gse_id": gse_id,
                "status": "failed",
                "error": f"No data file found. Searched: {gse_path}, data/geo/{gse_id}",
            }
        
        print(f"[Phase4 Worker {gse_id}] 📂 Loading: {input_path}")
        
        # 构建 State
        state = {
            "input_path": input_path,
            "hypothesis_id": hyp_id,
            "hypothesis_parser_gene_sets": gene_sets,
            "hypothesis_parser_entity1": cell_type,
            "params": {
                "gene_sets": gene_sets,
                "cell_type": cell_type,
                "min_genes": 200,
                "min_cells": 3,
                "max_mito_pct": 0.20,
                "n_top_hvg": 2000,
                "n_pcs": 50,
                "n_neighbors": 15,
                "leiden_resolution": 0.5,
            },
            "_enable_l2": False,  # 暂时禁用 L2 自愈（需要依赖）
        }
        
        # 执行 bioskills pipeline
        engine = SkillPipelineEngine()
        result_obj = engine.run(state, pipeline="cell_state_fast")
        
        if result_obj.status == "success" or result_obj.status == "partial":
            # 提取效应量结果
            final_state = result_obj.state
            effect_result = final_state.get("effect_size_result", {})
            
            result = {
                "gse_id": gse_id,
                "status": "success",
                "effect_size_cohen_d": effect_result.get("cohens_d", 0.0),
                "p_value": effect_result.get("p_value", 1.0),
                "direction": effect_result.get("direction", "unknown"),
                "n_t_cells": effect_result.get("n_group1", 0),
                "n_normal": effect_result.get("n_group2", 0),
                "method": "bioskills_pipeline",
                "pipeline": "cell_state_fast",
                "executed_skills": result_obj.executed,
                "failed_skills": result_obj.failed,
                "duration_seconds": result_obj.duration_seconds,
            }
            
            print(f"[Phase4 Worker {gse_id}] ✅ Done: d={result['effect_size_cohen_d']:.3f}, p={result['p_value']:.2e}")
            
        else:
            result = {
                "gse_id": gse_id,
                "status": "failed",
                "error": result_obj.error or "Pipeline failed",
                "failed_skills": result_obj.failed,
                "logs": result_obj.logs[-10:],  # 最后10条日志
            }
            print(f"[Phase4 Worker {gse_id}] ❌ Failed: {result['error']}")
        
        return result
        
    except Exception as e:
        import traceback
        tb = traceback.format_exc()[:500]
        print(f"[Phase4 Worker {gse_id}] ❌ Exception: {e}")
        return {
            "gse_id": gse_id,
            "status": "failed",
            "error": str(e),
            "traceback": tb,
        }


def effect_size_analyst_supervisor(state: BioAnalysisState) -> BioAnalysisState:
    """
    EffectSizeAnalyst Supervisor Node：
    - 读取待分析数据集列表
    - 从 hypothesis_parser 获取基因集和细胞类型
    - 并行启动 Worker（每个数据集一个线程）
    - 合并结果到 state
    """
    hyp_id = state["hypothesis_id"]
    registry = StateRegistry(hyp_id)
    
    new_state = dict(state)
    new_state["effect_size_analyst_status"] = "running"
    registry.set(new_state, phase=None)
    
    try:
        # ── Step 1: 获取基因集和细胞类型 ────────────
        gene_sets = state.get("hypothesis_parser_gene_sets") or state.get("literature_miner_gene_sets", {})
        cell_type = state.get("hypothesis_parser_entity1", "T cell")
        
        if not gene_sets:
            # 动态解析
            from bioskills.knowledge.gene_set_db import resolve_gene_set
            relation = state.get("hypothesis_parser_relation", "")
            hyp_text = state.get("hypothesis_text", "")
            gene_sets_list = resolve_gene_set(cell_type)  # returns List[str]
            gene_sets = {"CUSTOM": gene_sets_list}
            print(f"[EffectSizeAnalyst] GeneSetDB → {list(gene_sets.keys())} ({len(gene_sets_list)} genes)")
        
        # ── Step 2: 获取数据集列表 ──────────────────
        datasets = state.get("geo_retriever_datasets", [])
        if not datasets:
            manifest = registry.read_manifest()
            if manifest:
                datasets = manifest.get("datasets", [])
        
        if not datasets:
            print("[EffectSizeAnalyst] ⚠ No datasets to analyze, skipping to MetaAnalyst")
            new_state.update({
                "effect_size_analyst_status": "done",
                "effect_size_analyst_n_analyzed": 0,
                "effect_size_analyst_n_total": 0,
                "current_phase": "meta_analyst",
            })
            if "effect_size_analyst" not in new_state["checkpoints"]:
                new_state["checkpoints"] = new_state["checkpoints"] + ["effect_size_analyst"]
            registry.set(new_state, phase="effect_size_analyst")
            return new_state
        
        # ── Step 3: 过滤已完成的 ────────────────────
        # 构建 manifest_paths: gse_id -> path
        manifest_paths = {}
        for d in datasets:
            if d.get("path") and d.get("status") in ("downloaded", "raw_tar"):
                manifest_paths[d["gse_id"]] = d["path"]
        print(f"[EffectSizeAnalyst] manifest_paths: {list(manifest_paths.keys())[:5]}")

        existing_results = registry.list_phase4_results()
        existing_ids = {r.get("gse_id") for r in existing_results if r.get("status") != "failed"}
        
        to_analyze = [d["gse_id"] for d in datasets 
                      if d.get("gse_id") not in existing_ids]
        
        already_done = len(datasets) - len(to_analyze)
        print(f"[EffectSizeAnalyst] 📊 {len(datasets)} total, {already_done} done, {len(to_analyze)} to analyze")
        
        results_map = {r.get("gse_id"): r for r in existing_results}
        errors = []
        
        # ── Step 4: 并行执行 ────────────────────────
        MAX_WORKERS = min(4, len(to_analyze))  # 降低并发数避免内存问题
        
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {
                executor.submit(_analyze_single_dataset_bioskills, gse_id, hyp_id, gene_sets, cell_type, manifest_paths): gse_id
                for gse_id in to_analyze
            }
            
            for future in as_completed(futures):
                gse_id = futures[future]
                try:
                    result = future.result(timeout=600)
                    
                    registry.write_dataset_result(gse_id, result)
                    
                    if result.get("status") == "failed":
                        errors.append({"gse_id": gse_id, "error": result.get("error", "unknown")})
                    else:
                        results_map[gse_id] = result
                    
                    done = len(results_map)
                    total = len(datasets)
                    print(f"[EffectSizeAnalyst] 📈 Progress: {done}/{total} ({done/total*100:.0f}%)")
                    
                except Exception as e:
                    err = {"gse_id": gse_id, "error": str(e)}
                    errors.append(err)
                    print(f"[Phase4 Worker {gse_id}] ❌ Exception: {e}")
        
        # ── Step 5: 更新状态 ────────────────────────
        new_state.update({
            "effect_size_analyst_results": results_map,
            "effect_size_analyst_n_analyzed": len(results_map),
            "effect_size_analyst_n_total": len(datasets),
            "effect_size_analyst_errors": errors,
            "effect_size_analyst_status": "done",
            "effect_size_analyst_error": None,
            "current_phase": "meta_analyst",
            "hypothesis_parser_gene_sets": gene_sets,  # 确保基因集写入 state
        })
        if "effect_size_analyst" not in new_state["checkpoints"]:
            new_state["checkpoints"] = new_state["checkpoints"] + ["effect_size_analyst"]
        registry.set(new_state, phase="effect_size_analyst")
        
        print(f"[EffectSizeAnalyst] ✅ Done: {len(results_map)} analyzed, {len(errors)} failed")
        
    except Exception as e:
        new_state.update({
            "effect_size_analyst_status": "failed",
            "effect_size_analyst_error": str(e),
        })
        new_state.setdefault("errors", []).append(f"Phase4: {e}")
        registry.set(new_state, phase=None)
        print(f"[EffectSizeAnalyst] ❌ Failed: {e}")
    
    return new_state
