"""
L1: MSigDBClient Skill — MSigDB 基因集下载与缓存

包装 bioskills.knowledge.msigdb_client.MSigDBClient 为原子技能。
下载 Broad Institute MSigDB 官方基因集（.gmt 格式）。

输出:
  - gene_sets: {name: [genes]} dict
  - collection_info: 集合元数据
  - download_report: {collection, n_gene_sets, cache_hit}
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State


@register
class MSigDBClientSkill(AbstractSkill):
    """
    Download and query MSigDB gene set collections.
    
    常用 collection：
      - h: Hallmark (50 most relevant pathways)
      - c2: Curated (KEGG, Reactome, BioCarta)
      - c5: GO (Gene Ontology)
      - c7: Immunologic signatures (immune-related)
      - c8: Cell type signatures (cell type markers)
    
    参数（params）：
      - collection: str  # "h" | "c2" | "c5" | "c7" | "c8"
      - action: str  # "download" | "search" | "get_gene_set"
      - query: str  # 搜索关键词（action=search 时）
      - gene_set_name: str  # 特定基因集名（action=get_gene_set 时）
      - species: str  # "human" | "mouse"
    """
    
    name = "msigdb_client"
    description = (
        "Download MSigDB gene sets (Hallmark, KEGG, GO, Immune signatures) "
        "from Broad Institute. Auto-caches to ~/.cache/bioskills/msigdb/. "
        "bioSkills: knowledge base download before enrichment analysis"
    )
    input_contract = []  # standalone（无 adata 也可运行）
    output_contract = ["gene_sets", "download_report"]
    stage = Stage.ENRICHMENT
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.SPATIAL, Modality.GENERIC]
    
    tunable_parameters = {
        "collection": {"type": "str", "default": "h"},
        "action": {"type": "str", "default": "download"},
        "query": {"type": "str", "default": ""},
        "gene_set_name": {"type": "str", "default": ""},
        "species": {"type": "str", "default": "human"},
        "cache_dir": {"type": "str", "default": ""},
    }
    
    def _run(self, state: State) -> dict:
        from bioskills.knowledge.msigdb_client import MSigDBClient
        params = state.get("params", {})
        
        collection = params.get("collection", "h")
        action = params.get("action", "download")
        query = params.get("query", "")
        gene_set_name = params.get("gene_set_name", "")
        species = params.get("species", "human")
        cache_dir = params.get("cache_dir", "")
        
        client_kwargs = {"timeout": 30}
        if cache_dir:
            client_kwargs["cache_dir"] = cache_dir
        
        client = MSigDBClient(**client_kwargs)
        gene_sets = {}
        cache_hit = False
        
        if action == "download":
            # 检查缓存
            import os
            from pathlib import Path
            cache_path = Path(
                os.path.expanduser(cache_dir) if cache_dir 
                else f"~/.cache/bioskills/msigdb/{collection}_{species}.gmt"
            )
            cache_hit = cache_path.exists()
            
            gene_sets = client.download(collection=collection, species=species)
            
        elif action == "search":
            gene_sets = client.search(query=query, collection=collection, species=species)
            
        elif action == "get_gene_set":
            genes = client.get_gene_set(gene_set_name, collection=collection, species=species)
            gene_sets = {gene_set_name: genes}
        
        report = {
            "status": "success",
            "action": action,
            "collection": collection,
            "species": species,
            "n_gene_sets": len(gene_sets),
            "cache_hit": cache_hit,
            "gene_set_names": list(gene_sets.keys())[:20],  # 前20个
        }
        
        print(f"  [MSigDB] ✅ {action}: {len(gene_sets)} gene sets "
              f"from collection '{collection}' (cache_hit={cache_hit})")
        
        return {"gene_sets": gene_sets, "download_report": report}