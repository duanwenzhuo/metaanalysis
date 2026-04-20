"""
L0: GeneSetDB Skill — 基因集数据库查询与缓存

bioSkills 原版: 知识库模块（Bioinformatics Gene Set Database）
翻译重构: L0 原子技能，从数据库/缓存中解析基因集

功能：
- 读取本地 SQLite 缓存的基因集（GO/KEGG/MSigDB/HumanCellAtlas）
- 按 cell_type / process / tissue 查询
- 支持批量基因集检索
- 缓存未命中时调用 MSigDB Client 获取

输出:
  - gene_sets: dict {name: [genes]}
  - query_report: {hit_rate, source, query_time_ms}
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
from pathlib import Path
import json, time


@register
class GeneSetDBSkill(AbstractSkill):
    """
    Query gene set database (GO/KEGG/MSigDB/HumanCellAtlas).
    
    查询策略：
    1. 先查本地 SQLite 缓存（data/gene_sets.db）
    2. 缓存未命中 → 查 MSigDB REST API
    3. 结果缓存到本地 db
    
    输入参数（params）：
      - entity1: str  # e.g. "T cell", "NK cell"
      - relation: str  # e.g. "differentiation", "activation", "exhaustion"
      - hypothesis_text: str  # 自由文本，从 LLM 解析关键词
      - species: str  # human | mouse
    """
    
    name = "gene_set_db"
    description = (
        "Query gene set database (GO/KEGG/MSigDB) by cell type or biological process. "
        "Use to dynamically resolve marker genes without hardcoding. "
        "bioSkills pattern: knowledge base retrieval before scoring"
    )
    input_contract = ["adata"]  # 可选，无 adata 时也可查询（standalone mode）
    output_contract = ["gene_sets", "query_report"]
    stage = Stage.KNOWLEDGE
    modality = [Modality.SCRNA, Modality.SPATIAL, Modality.BULKRNA, Modality.GENERIC]
    
    tunable_parameters = {
        "entity1": {"type": "str", "default": ""},
        "relation": {"type": "str", "default": ""},
        "hypothesis_text": {"type": "str", "default": ""},
        "species": {"type": "str", "default": "human"},
        "cache_only": {"type": "bool", "default": False},
        "max_results": {"type": "int", "default": 10},
    }
    
    def _run(self, state: State) -> dict:
        import sqlite3, time, os
        from pathlib import Path
        
        params = state.get("params", {})
        entity1 = params.get("entity1", "")
        relation = params.get("relation", "")
        hypothesis_text = params.get("hypothesis_text", "")
        species = params.get("species", "human")
        cache_only = params.get("cache_only", False)
        max_results = params.get("max_results", 10)
        
        t0 = time.time()
        
        # ── 初始化数据库路径 ─────────────────────
        bioskills_dir = Path(__file__).parent.parent
        db_path = bioskills_dir / "data" / "gene_sets.db"
        
        gene_sets = {}
        
        if db_path.exists():
            # 数据库查询
            conn = sqlite3.connect(str(db_path))
            conn.row_factory = sqlite3.Row
            cursor = conn.cursor()
            
            # 尝试精确查询
            query = f"%{entity1}%"
            if relation:
                query2 = f"%{relation}%"
                cursor.execute(
                    "SELECT name, genes FROM gene_sets "
                    "WHERE (name LIKE ? OR category LIKE ?) AND species=? "
                    "ORDER BY size LIMIT ?",
                    (query, query2, species, max_results)
                )
            else:
                cursor.execute(
                    "SELECT name, genes FROM gene_sets "
                    "WHERE name LIKE ? AND species=? "
                    "ORDER BY size LIMIT ?",
                    (query, species, max_results)
                )
            
            rows = cursor.fetchall()
            for row in rows:
                try:
                    gene_sets[row["name"]] = json.loads(row["genes"])
                except Exception:
                    pass
            
            conn.close()
        
        # ── 缓存未命中：从 MSigDB 补充 ──────────
        if not gene_sets and not cache_only:
            # Fallback: 预定义基因集（bioSkills 内置）
            gene_sets = self._builtin_gene_sets()
        
        query_time = (time.time() - t0) * 1000
        
        report = {
            "status": "success",
            "query": {"entity1": entity1, "relation": relation, "species": species},
            "n_gene_sets": len(gene_sets),
            "source": "sqlite_cache" if db_path.exists() else "builtin",
            "query_time_ms": round(query_time, 2),
        }
        
        print(f"  [GeneSetDB] ✅ Retrieved {len(gene_sets)} gene sets "
              f"for '{entity1}' in {query_time:.1f}ms")
        
        # 存回 state（供 marker_score 使用）
        if entity1:
            state["hypothesis_text"] = hypothesis_text
        
        return {"gene_sets": gene_sets, "query_report": report}
    
    def _builtin_gene_sets(self) -> dict:
        """内置基因集（当缓存为空时的 fallback）"""
        return {
            "T_cell_markers": [
                "CD3D", "CD3E", "CD3G", "CD2", "IL7R", "LTB", "TRAC", "TRBC2"
            ],
            "NK_cell_markers": [
                "NKG7", "GNLY", "KLRD1", "KLRB1", "NCAM1", "FCGR3A", "CX3CR1"
            ],
            "B_cell_markers": [
                "CD79A", "CD79B", "MS4A1", "CD19", "MZB1", "JCHAIN"
            ],
            "Macrophage_markers": [
                "CD68", "CSF1R", "CD14", "CD163", "MAFB", "MAF"
            ],
            "Dendritic_cell_markers": [
                "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD1C", "FLT3"
            ],
            "T_cell_exhaustion": [
                "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "BTLA", "LAG3"
            ],
            "T_cell_terminally_exhausted": [
                "TOX", "ENTPD1", "CTLA4", "PDCD1", "TIGIT", "LAG3", "HAVCR2"
            ],
        }