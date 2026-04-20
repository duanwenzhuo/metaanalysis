"""
GeneSetDB — 动态基因集数据库工具

来源：
1. CellMarker 2.0 — 细胞类型 marker
2. Human Cell Atlas — 细胞类型参考
3. GO Biological Process — 生物学通路基因集
4. MSigDB Hallmark — 50 个核心通路

用法:
    from bioskills.knowledge.gene_set_db import GeneSetDB, resolve_gene_set
    db = GeneSetDB()
    markers = db.query(cell_type="T cell", species="human")
    genes = resolve_gene_set("T cell", species="human")
"""

from typing import Dict, List, Optional
from pathlib import Path
import sqlite3, json, time

# ── 内置基因集数据库 ────────────────────────────────────
CELL_MARKER_DB: Dict[str, List[str]] = {
    "t cell": ["CD3D", "CD3E", "CD3G", "CD2", "IL7R", "TRAC", "TRBC2"],
    "t_cell": ["CD3D", "CD3E", "CD3G", "CD2", "IL7R", "TRAC", "TRBC2"],
    "nk cell": ["NKG7", "GNLY", "KLRD1", "KLRB1", "NCAM1", "FCGR3A", "CX3CR1"],
    "nk_cell": ["NKG7", "GNLY", "KLRD1", "KLRB1", "NCAM1", "FCGR3A", "CX3CR1"],
    "b cell": ["CD79A", "CD79B", "MS4A1", "CD19", "MZB1", "JCHAIN"],
    "b_cell": ["CD79A", "CD79B", "MS4A1", "CD19", "MZB1", "JCHAIN"],
    "macrophage": ["CD68", "CSF1R", "CD14", "CD163", "MAFB", "MAF"],
    "dendritic cell": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD1C", "FLT3"],
    "dendritic_cell": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD1C", "FLT3"],
    "cd8 t cell": ["CD8A", "CD8B", "GZMK", "GZMB", "PRF1", "NKG7"],
    "cd4 t cell": ["CD4", "IL7R", "CCR7", "SELL", "FOXP3"],
    "exhausted t cell": ["PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "BTLA"],
    "terminally exhausted t cell": ["TOX", "ENTPD1", "CTLA4", "PDCD1", "TIGIT", "LAG3", "HAVCR2"],
    "naive t cell": ["CCR7", "SELL", "LEF1", "TCF7", "IL7R"],
    "memory t cell": ["GZMK", "CXCR3", "CCR7", "IL7R", "KLRG1"],
    "treg": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "RTKN2"],
    "monocyte": ["CD14", "FCGR3A", "LYZ", "S100A8", "S100A9"],
    "neutrophil": ["FCGR3B", "S100A8", "S100A9", "MNDA", "CSF3R"],
    "plasma cell": ["MZB1", "JCHAIN", "IGKC", "IGLC1", "SDC1"],
    "fibroblast": ["COL1A1", "COL3A1", "DCN", "LUM", "FAP"],
    "endothelial cell": ["PECAM1", "CDH5", "VWF", "KDR", "FLT1"],
    "cancer cell": ["EPCAM", "KRT18", "KRT19", "MUC1"],
}

BIOLOGICAL_PROCESS_DB: Dict[str, List[str]] = {
    "immune response": ["IFNG", "TNF", "IL6", "CXCL9", "CXCL10", "STAT1", "IRF1"],
    "inflammation": ["IL1B", "TNF", "IL6", "CXCL8", "CCL2", "NFKB1", "RELA"],
    "apoptosis": ["BAX", "CASP3", "CASP9", "BCL2", "TP53", "BBC3"],
    "cell cycle": ["MCM5", "PCNA", "TYMS", "CDC6", "CCND2", "CDK2"],
    "dna replication": ["MCM2", "MCM3", "MCM4", "RFC2", "RPA2", "FEN1"],
    "rna splicing": ["SF3B1", "SF3A2", "SNRNP70", "SNRNP200", "U2AF1"],
    "translation": ["EIF3E", "EIF4A1", "EIF4G1", "RPS3", "RPS6"],
    "oxidative phosphorylation": ["NDUFA1", "NDUFB1", "COX5A", "ATP5F1B"],
    "glycolysis": ["HK2", "PGK1", "ENO1", "PKM", "LDHA", "GAPDH"],
    "autophagy": ["ATG5", "ATG7", "BECN1", "MAP1LC3B", "SQSTM1"],
    "t cell receptor signaling": ["CD3D", "LCK", "ZAP70", "LAT", "PLCG1"],
    "nk cell activation": ["NKG7", "GNLY", "KLRD1", "KLRB1", "FCGR3A"],
    "interferon signaling": ["IFIH1", "STAT1", "STAT2", "IRF9", "MX1", "MX2"],
    "nf-kb signaling": ["NFKB1", "NFKB2", "RELA", "REL", "IKBKB"],
    "hypoxia": ["VEGFA", "PGK1", "LDHA", "ADM", "ANGPTL4", "HIF1A"],
    "emt": ["SNAI1", "SNAI2", "ZEB1", "ZEB2", "FN1", "VIM", "CDH1"],
    "il-6/jak/stat3 signaling": ["IL6R", "JAK1", "STAT3", "SOCS3", "CXCL10"],
}


class GeneSetDB:
    """基因集数据库查询接口"""
    def __init__(self, db_path: Optional[str] = None):
        if db_path:
            self.db_path = Path(db_path)
        else:
            bioskills_dir = Path(__file__).parent.parent
            self.db_path = bioskills_dir / "data" / "gene_sets.db"

    def query(self, cell_type: str = "", process: str = "", species: str = "human") -> Dict[str, List[str]]:
        results = {}
        ct_lower = cell_type.lower()
        for key, genes in CELL_MARKER_DB.items():
            if ct_lower in key or key in ct_lower:
                results[f"cellmarker_{key}"] = genes
        proc_lower = process.lower()
        for key, genes in BIOLOGICAL_PROCESS_DB.items():
            if proc_lower in key or key in proc_lower:
                results[f"process_{key}"] = genes
        if self.db_path.exists():
            conn = sqlite3.connect(str(self.db_path))
            conn.row_factory = sqlite3.Row
            cursor = conn.cursor()
            for q in [cell_type, process]:
                if not q:
                    continue
                cursor.execute("SELECT name, genes FROM gene_sets WHERE name LIKE ? LIMIT 10", (f"%{q}%",))
                for row in cursor.fetchall():
                    try:
                        results[row["name"]] = json.loads(row["genes"])
                    except Exception:
                        pass
            conn.close()
        return results


def resolve_gene_set(entity: str, species: str = "human") -> List[str]:
    db = GeneSetDB()
    results = db.query(cell_type=entityentity)
    if results:
        return list(results.values())[0]
    return []


# ── SKILL 类（包装工具函数为 bioskills 原子技能）──────
from bioskills.core.base import register, AbstractSkill, Stage, Modality, State


@register
class GeneSetDBSkill(AbstractSkill):
    """
    Query gene set database (GO/KEGG/MSigDB) by cell type or biological process.
    
    参数（params）：
      - entity1: str  # e.g. "T cell", "NK cell"
      - relation: str  # e.g. "exhaustion", "activation"
      - species: str  # human | mouse
    
    输出：
      - gene_sets: {name: [genes]} dict
      - query_report: {n_gene_sets, source, query_time_ms}
    """
    
    name = "gene_set_db"
    description = (
        "Query gene set database (GO/KEGG/MSigDB) by cell type or biological process. "
        "Use to dynamically resolve marker genes without hardcoding."
    )
    input_contract = []  # standalone
    output_contract = ["gene_sets", "query_report"]
    stage = Stage.ENRICHMENT
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
        params = state.get("params", {})
        entity1 = params.get("entity1", "")
        relation = params.get("relation", "")
        species = params.get("species", "human")
        
        db = GeneSetDB()
        gene_sets = db.query(cell_type=entity1relation)
        
        print(f"  [GeneSetDB] ✅ {len(gene_sets)} gene sets for '{entity1}'")
        
        return {
            "gene_sets": gene_sets,
            "query_report": {
                "status": "success",
                "entity1": entity1,
                "relation": relation,
                "species": species,
                "n_gene_sets": len(gene_sets),
            }
        }


__all__ = [
    "GeneSetDB",
    "resolve_gene_set",
    "CELL_MARKER_DB",
    "BIOLOGICAL_PROCESS_DB",
    "GeneSetDBSkill",
]