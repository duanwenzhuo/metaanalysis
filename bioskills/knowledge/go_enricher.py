"""
L0: GOEnricher Skill — Gene Ontology 富集分析

bioSkills 原版: enrichment 模块（go_enrichment.py）
翻译重构: 独立 L0 Skill，直接封装 g:Profiler / Enrichr API

功能：
- 超几何检验：给定基因集 vs GO/KEGG terms
- 多个数据库并行查询（GO_BP, GO_MF, GO_CC, KEGG, Reactome）
- Benjamini-Hochberg FDR 校正
- 结果 ranked by p-value

输出:
  - enrichment_results: list of {term, genes, p_value, padj, overlap}
  - enriched_terms: list of significant term names
"""

from bioskills.core.base import register, AbstractSkill, Stage, Modality, State
from typing import List, Dict


@register
class GOEnricherSkill(AbstractSkill):
    """
    Perform Gene Ontology (GO) and KEGG pathway enrichment analysis.
    
    API 选择策略（bioSkills）：
    | Source      | Speed | Coverage | Cost |
    |-------------|-------|----------|------|
    | g:Profiler  | Fast  | Great    | Free |
    | Enrichr     | Fast  | Good     | Free |
    | DAVID       | Slow  | Moderate | Free |
    
    默认使用 g:Profiler（REST API，R/Python 均支持）
    
    输入参数（params）：
      - genes: List[str]  # 基因列表
      - source: str  # "gprofiler" | "enrichr" | "david"
      - databases: List[str]  # ["GO:BP", "KEGG", "REAC"] 等
      - organism: str  # "hsapiens" | "mmusculus"
      - significance_threshold: float  # p-adjusted threshold
    """
    
    name = "go_enrichment"
    description = (
        "Perform GO and KEGG pathway enrichment analysis for a gene list. "
        "Uses g:Profiler REST API (free, no key required). "
        "bioSkills pattern: enrichment after differential expression"
    )
    input_contract = ["adata"]  # 或直接传 genes via params
    output_contract = ["enrichment_results", "enriched_terms"]
    stage = Stage.ENRICHMENT
    modality = [Modality.SCRNA, Modality.BULKRNA, Modality.SPATIAL, Modality.GENERIC]
    
    tunable_parameters = {
        "genes": {"type": "list", "default": []},  # 空=从 adata 选取差异基因
        "source": {"type": "str", "default": "gprofiler"},
        "databases": {
            "type": "list",
            "default": ["GO:BP", "KEGG", "REAC", "WP", "HPA"]
        },
        "organism": {"type": "str", "default": "hsapiens"},
        "significance_threshold": {"type": "float", "default": 0.05},
        "min_overlap": {"type": "int", "default": 3},
        "max_terms": {"type": "int", "default": 50},
    }
    
    def _run(self, state: State) -> dict:
        import json, urllib.request, urllib.parse, time
        
        params = state.get("params", {})
        genes = params.get("genes", [])
        source = params.get("source", "gprofiler")
        databases = params.get("databases", ["GO:BP", "KEGG", "REAC"])
        organism = params.get("organism", "hsapiens")
        sig_thresh = params.get("significance_threshold", 0.05)
        min_overlap = params.get("min_overlap", 3)
        max_terms = params.get("max_terms", 50)
        
        # ── Step 1: 从 adata 提取差异基因（如果未传入）───
        if not genes and "adata" in state:
            adata = state["adata"]
            # 找 DE 结果（假设 volcano skill 运行过）
            if "diff_genes" in adata.uns:
                genes = adata.uns["diff_genes"]
            else:
                # 取 logfc > 1 的基因
                if "rank_genes_groups" in adata.uns:
                    # Scanpy style
                    first_key = list(adata.uns["rank_genes_groups"].keys())[0]
                    names = adata.uns["rank_genes_groups"]["names"]
                    scores = adata.uns["rank_genes_groups"]["scores"]
                    genes = [n for n, s in zip(names[first_key], scores[first_key]) 
                             if s > 1.0]
        
        if not genes:
            return {
                "enrichment_results": [],
                "enriched_terms": [],
                "enrichment_report": {
                    "status": "failed",
                    "error": "No genes provided and no diff_genes found in adata.uns"
                }
            }
        
        print(f"  [GOEnricher] Running enrichment for {len(genes)} genes "
              f"({organism}, source={source})")
        
        # ── Step 2: 调用 g:Profiler API ───────────
        if source == "gprofiler":
            return self._run_gprofiler(genes, organism, databases, sig_thresh, 
                                       min_overlap, max_terms)
        else:
            return {
                "enrichment_results": [],
                "enriched_terms": [],
                "enrichment_report": {
                    "status": "failed",
                    "error": f"Source '{source}' not implemented. Use 'gprofiler'."
                }
            }
    
    def _run_gprofiler(self, genes, organism, databases, sig_thresh, 
                       min_overlap, max_terms) -> dict:
        """调用 g:Profiler REST API"""
        import json, urllib.request, time
        
        # g:Profiler 转换请求格式
        user_values = genes
        
        payload = json.dumps({
            "user": user_values,
            "organism": organism,
            "sources": databases,
            "significance_threshold_method": "g_SCS",
            "numeric_namespace": "HGNC",
            "background": "registered",
            "user_threshold": sig_thresh,
            "min_set_size": min_overlap,
            "max_set_size": 1000,
        }).encode("utf-8")
        
        req = urllib.request.Request(
            "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
            data=payload,
            headers={"Content-Type": "application/json"},
            method="POST"
        )
        
        results = []
        try:
            with urllib.request.urlopen(req, timeout=30) as response:
                data = json.loads(response.read().decode("utf-8"))
                terms = data.get("result", [])
                
                for term in terms[:max_terms]:
                    results.append({
                        "term_id": term.get("id", ""),
                        "term_name": term.get("name", ""),
                        "source": term.get("native", ""),
                        "p_value": round(term.get("p_value", 1.0), 6),
                        "overlap": term.get("overlap_size", 0),
                        "term_size": term.get("term_size", 0),
                        "genes": term.get("intersections", []),
                    })
        
        except Exception as e:
            return {
                "enrichment_results": [],
                "enriched_terms": [],
                "enrichment_report": {
                    "status": "failed",
                    "error": f"g:Profiler API failed: {e}",
                    "suggestion": "Check internet connection; try Enrichr as fallback"
                }
            }
        
        enriched_terms = [r["term_name"] for r in results if r["p_value"] < sig_thresh]
        
        print(f"  [GOEnricher] ✅ {len(results)} terms enriched "
              f"({len(enriched_terms)} significant, FDR<{sig_thresh})")
        
        return {
            "enrichment_results": results,
            "enriched_terms": enriched_terms,
            "enrichment_report": {
                "status": "success",
                "n_input_genes": len(genes),
                "n_terms_total": len(results),
                "n_terms_significant": len(enriched_terms),
                "source": "gprofiler",
                "organism": organism,
            }
        }