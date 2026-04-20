#!/usr/bin/env python3
"""
Phase 2: Literature Research & Gene Set Builder

从文献中提取基因集 + 优先级数据集列表，输出 phase2_gene_set.json。

Usage:
    python scripts/phase2_literature_miner.py "T cell exhaustion in tumor microenvironment"
    python scripts/phase2_literature_miner.py --hypothesis "Macrophage M2 polarization in TME" --output data/h_xxx/
"""

import os
import sys
import re
import json
import time
import argparse
import gzip
import urllib.request
import urllib.parse
import urllib.error
from datetime import datetime
from pathlib import Path
from collections import Counter, defaultdict
from typing import List, Dict, Any, Optional, Tuple

# ─────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────

def setup_proxy():
    """设置代理"""
    for key in ["HTTP_PROXY", "HTTPS_PROXY", "http_proxy", "https_proxy"]:
        os.environ[key] = os.environ.get(key) or "http://127.0.0.1:7890"


def fetch_url(url: str, timeout: int = 20) -> Optional[str]:
    """带代理的URL抓取"""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", errors="replace")
    except Exception as e:
        print(f"  ⚠ fetch failed: {url[:60]}... → {e}")
        return None


def extract_geo_ids_from_text(text: str) -> List[str]:
    """从文本中提取GEO ID（支持GSE和GSM）"""
    gse_ids = re.findall(r'GSE\d{5,8}', text.upper())
    return list(dict.fromkeys(gse_ids))


def extract_pmids_from_text(text: str) -> List[str]:
    """从文本中提取PMID"""
    return re.findall(r'\b(?:PMID|DOI)[:\s]*(\d{7,10})\b', text, re.IGNORECASE)


def fetch_pubmed_abstracts(pmids: List[str]) -> Dict[str, Dict]:
    """批量获取PubMed摘要"""
    results = {}
    batch_size = 100
    
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i+batch_size]
        ids_str = ",".join(batch)
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={ids_str}&rettype=abstract&retmode=xml"
        
        xml_text = fetch_url(url, timeout=30)
        if not xml_text:
            continue
        
        try:
            import xml.etree.ElementTree as ET
            root = ET.fromstring(xml_text)
            for article in root.findall(".//PubmedArticle"):
                pmid_el = article.find(".//PMID")
                if pmid_el is None:
                    continue
                pmid = pmid_el.text
                
                title_el = article.find(".//ArticleTitle")
                title = title_el.text if title_el is not None else ""
                
                abstract_el = article.find(".//AbstractText")
                abstract = abstract_el.text if abstract_el is not None else ""
                
                # 提取期刊
                journal_el = article.find(".//Journal/Title")
                journal = journal_el.text if journal_el is not None else ""
                
                results[pmid] = {
                    "title": title,
                    "abstract": abstract,
                    "journal": journal,
                    "pmid": pmid,
                    "text": f"{title} {abstract}"
                }
        except Exception as e:
            print(f"  ⚠ XML解析失败 (batch {i//batch_size}): {e}")
    
    return results


def fetch_geo_metadata(gse_id: str) -> Dict[str, Any]:
    """获取GEO元数据"""
    metadata = {
        "accession": gse_id,
        "exists": False,
        "title": None,
        "summary": None,
        "n_samples": 0,
        "organism": None,
        "platform": None,
        "technology": None,
        "disease": None,
    }
    
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={gse_id}&retmode=json"
    data = fetch_url(url)
    if not data:
        return metadata
    
    try:
        import json as json_lib
        j = json_lib.loads(data)
        result = j.get("result", {}).get(gse_id, {})
        if result:
            metadata["exists"] = True
            metadata["title"] = result.get("title", "")
            metadata["summary"] = result.get("summary", "")
            metadata["n_samples"] = result.get("n_samples", 0)
            metadata["organism"] = result.get("taxon", "")
            
            # 提取平台
            gpl = result.get("gpl", "")
            metadata["platform"] = gpl
            
            # 判断技术
            text = f"{metadata['title']} {metadata['summary']}".lower()
            if any(x in text for x in ["10x", "chromium", "single cell", "scrna", "single-cell", "droplet"]):
                metadata["technology"] = "10x/scRNA-seq"
            elif any(x in text for x in ["smart-seq", "smartseq", "facs", "sort"]):
                metadata["technology"] = "Smart-seq2/FACS"
            elif any(x in text for x in ["spatial", "visium", "slide-seq", "st"]):
                metadata["technology"] = "Spatial"
            elif any(x in text for x in ["microarray", "affymetrix", "agilent"]):
                metadata["technology"] = "Microarray"
            else:
                metadata["technology"] = "Unknown"
            
            # 提取疾病
            disease_kw = {
                "melanoma": "Melanoma", "breast cancer": "Breast Cancer",
                "lung cancer": "Lung Cancer", "nsclc": "NSCLC",
                "colorectal": "Colorectal Cancer", "colon cancer": "Colorectal Cancer",
                "pancreatic": "Pancreatic Cancer", "pdac": "Pancreatic Cancer",
                "gastric": "Gastric Cancer", "liver": "Liver Cancer", "hcc": "HCC",
                "ovarian": "Ovarian Cancer", "kidney": "Renal Cell Carcinoma",
                "rcc": "Renal Cell Carcinoma", "prostate": "Prostate Cancer",
                "glioma": "Glioma", "glioblastoma": "Glioblastoma",
                "head and neck": "Head and Neck SCC", "hnscc": "HNSCC",
                "endometrial": "Endometrial Cancer", "cervical": "Cervical Cancer",
            }
            for kw, disease in disease_kw.items():
                if kw in text:
                    metadata["disease"] = disease
                    break
    except Exception as e:
        print(f"  ⚠ GEO metadata parse failed {gse_id}: {e}")
    
    return metadata


def search_pubmed(query: str, max_results: int = 20) -> List[str]:
    """搜索PubMed，返回PMID列表"""
    encoded = urllib.parse.quote(query)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={encoded}&retmax={max_results}&retmode=json&sort=relevance"
    
    data = fetch_url(url, timeout=20)
    if not data:
        return []
    
    try:
        import json as json_lib
        j = json_lib.loads(data)
        return j.get("esearchresult", {}).get("idlist", [])
    except:
        return []


def search_geo(keywords: str, max_results: int = 15) -> List[str]:
    """搜索GEO数据库"""
    encoded = urllib.parse.quote(keywords)
    url = f"https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&search={encoded}&zsort=date&display={max_results}"
    
    html = fetch_url(url)
    if not html:
        return []
    
    gse_ids = re.findall(r'GSE\d{5,}', html.upper())
    return list(dict.fromkeys(gse_ids))[:max_results]


# ─────────────────────────────────────────────
# 主逻辑
# ─────────────────────────────────────────────

class LiteratureMiner:
    """
    Phase 2: 文献挖掘，输出 gene_set.json + recommended_datasets (P0/P1/P2)
    """
    
    def __init__(self, hypothesis_text: str, output_dir: str):
        self.hypothesis_text = hypothesis_text
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logs_dir = Path(os.path.dirname(__file__)).parent / "logs"
        self.logs_dir.mkdir(exist_ok=True)
        
        self.geo_ids_counter = Counter()  # 每个GSE被多少篇文献引用
        self.geo_to_papers = defaultdict(list)  # GSE → 论文列表
        self.all_papers = {}  # PMID → 论文详情
        self.gene_sets = {}  # 基因集
        self.recommended_datasets = []  # 优先级数据集
        self.failed_pmids = []
        self.failed_geo = []
        self.search_queries = []
        
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    def log(self, msg: str):
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] {msg}")
        log_path = self.logs_dir / f"phase2_{self.timestamp}.log"
        with open(log_path, "a") as f:
            f.write(f"[{ts}] {msg}\n")
    
    # ── Step 1: 生成搜索词 ──────────────────────
    
    def generate_search_queries(self) -> List[str]:
        """根据假设生成多个搜索词"""
        h = self.hypothesis_text.lower()
        
        queries = []
        
        if "exhaust" in h and ("t cell" in h or "cd8" in h or "tumor" in h):
            queries = [
                "T cell exhaustion tumor microenvironment scRNA-seq",
                "PD-1 CD8 T cell exhaustion single cell cancer",
                "HAVCR2 LAG3 TIGIT exhausted T cells tumor",
                "tumor infiltrating lymphocytes exhaustion scRNA",
            ]
            self.gene_sets = {
                "T_EXHAUSTION": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "BTLA", "CD160", "KLRG1"],
                "T_CELL_IDENTITY": ["CD3E", "CD3D", "CD8A", "CD4"],
                "TUMOR_MICROENVIRONMENT": ["EPCAM", "KRT8", "KRT18"],
            }
        elif "macrophage" in h or ("m2" in h and "polariz" in h):
            queries = [
                "tumor associated macrophage M2 polarization scRNA-seq",
                "TAM M2 CD163 MRC1 ARG1 tumor microenvironment",
                "macrophage polarization cancer single cell RNA",
                "SPP1 CD44 macrophage tumor immune microenvironment",
            ]
            self.gene_sets = {
                "M2_POLARIZATION": ["CD163", "MRC1", "ARG1", "IL10", "TGFB1", "CD200R1", "FOLR2", "SEPP1", "CCL18", "CCL22"],
                "M1_POLARIZATION": ["NOS2", "IL12B", "CXCL9", "CXCL10", "TNF", "IL1B"],
                "MACROPHAGE_IDENTITY": ["CD68", "CSF1R", "LYZ", "AIF1"],
                "SPP1_AXIS": ["SPP1", "CD44", "ITGB1", "ITGA4", "ITGAV"],
            }
        elif "fibroblast" in h or "caf" in h:
            queries = [
                "cancer associated fibroblast CAF scRNA-seq tumor",
                "myCAF iCAF apCAF fibroblast heterogeneity cancer",
                "tumor fibroblast CAF subtypes single cell",
                "CAF FAP COL1A1 PDGFRA cancer microenvironment",
            ]
            self.gene_sets = {
                "myCAF": ["ACTA2", "TAGLN", "MYH11", "TPM1", "COL1A1", "COL3A1", "FAP", "POSTN", "MMP11"],
                "iCAF": ["IL6", "CXCL12", "CXCL14", "CCL2", "PDGFRA", "HAS1", "LMNA"],
                "apCAF": ["CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1"],
                "CAF_PAN": ["FAP", "PDGFRA", "PDGFRB", "COL1A2", "DCN", "LUM", "VIM"],
            }
        elif "nk cell" in h:
            queries = [
                "NK cell tumor microenvironment scRNA-seq",
                "natural killer cell exhaustion cancer single cell",
                "NK cell dysfunction tumor infiltrating scRNA",
            ]
            self.gene_sets = {
                "NK_CELL_IDENTITY": ["NCR1", "NCR2", "NCR3", "KLRD1", "KLRF1", "GNLY", "NKG7"],
                "NK_ACTIVATION": ["IFNG", "TNF", "GZMB", "PRF1", "FCGR3A"],
                "NK_EXHAUSTION": ["PDCD1", "HAVCR2", "LAG3", "TIGIT"],
            }
        elif "dendritic" in h or "dc" in h:
            queries = [
                "dendritic cell tumor microenvironment scRNA-seq",
                "cDC1 cDC2 plasmacytoid DC cancer single cell",
                "DC antigen presentation tumor immune scRNA",
            ]
            self.gene_sets = {
                "cDC1": ["XCR1", "CLEC9A", "IRF8", "BATF3", "CD141"],
                "cDC2": ["CD1C", "FCER1A", "IRF4", "CD11B"],
                "pDC": ["IRF7", "IRF8", "LILRA4", "CLEC4C"],
                "DC_ACTIVATION": ["CD80", "CD83", "CD86", "HLA-DRA"],
            }
        else:
            # 通用癌症免疫关键词
            queries = [
                f"{h} scRNA-seq tumor microenvironment",
                f"{h} single cell RNA-seq cancer",
                f"tumor immune microenvironment {h} scRNA",
            ]
            self.gene_sets = {"CUSTOM": []}
        
        self.search_queries = queries
        self.log(f"Generated {len(queries)} search queries")
        for q in queries:
            self.log(f"  Q: {q}")
        
        return queries
    
    # ── Step 2: 搜索PubMed ─────────────────────
    
    def search_literature(self, queries: List[str]) -> List[str]:
        """搜索文献，返回所有PMID"""
        all_pmids = set()
        
        for query in queries:
            self.log(f"Searching PubMed: {query[:60]}")
            pmids = search_pubmed(query, max_results=20)
            self.log(f"  Found {len(pmids)} PMIDs")
            all_pmids.update(pmids)
            time.sleep(0.3)  # NCBI API限速
        
        self.log(f"Total unique PMIDs: {len(all_pmids)}")
        return list(all_pmids)
    
    # ── Step 3: 提取GEO ID ──────────────────────
    
    def extract_geo_ids_from_papers(self, pmids: List[str]) -> Dict[str, List[str]]:
        """从论文中提取GEO ID"""
        self.log(f"Fetching abstracts for {len(pmids)} PMIDs...")
        
        papers = fetch_pubmed_abstracts(pmids)
        self.all_papers = papers
        self.log(f"Successfully fetched {len(papers)} papers")
        
        # 从每篇论文中提取GEO ID
        for pmid, paper in papers.items():
            text = paper.get("text", "")
            gse_ids = extract_geo_ids_from_text(text)
            
            for gse_id in gse_ids:
                self.geo_ids_counter[gse_id] += 1
                self.geo_to_papers[gse_id].append({
                    "pmid": pmid,
                    "title": paper.get("title", ""),
                    "journal": paper.get("journal", ""),
                })
        
        self.log(f"Extracted GEO IDs from papers: {len(self.geo_ids_counter)}")
        return dict(self.geo_ids_counter)
    
    # ── Step 4: GEO搜索兜底 ─────────────────────
    
    def geo_search_fallback(self, queries: List[str]) -> List[str]:
        """GEO搜索兜底找更多数据集"""
        geo_ids = []
        
        for query in queries:
            self.log(f"Searching GEO: {query[:60]}")
            gse_ids = search_geo(query, max_results=10)
            geo_ids.extend(gse_ids)
            self.log(f"  Found {len(gse_ids)} GEO series")
            time.sleep(0.3)
        
        return list(dict.fromkeys(geo_ids))
    
    # ── Step 5: 验证数据集类型 ──────────────────
    
    def validate_and_enrich_datasets(self, gse_ids: List[str]) -> List[Dict]:
        """验证数据集类型，过滤非scRNA-seq"""
        validated = []
        
        for gse_id in gse_ids:
            self.log(f"Validating {gse_id}...")
            meta = fetch_geo_metadata(gse_id)
            
            if not meta["exists"]:
                self.log(f"  ✗ {gse_id} does not exist")
                self.failed_geo.append({"gse_id": gse_id, "reason": "Not found in GEO"})
                continue
            
            tech = meta.get("technology", "Unknown")
            
            # 优先选择 scRNA-seq
            if tech in ["10x/scRNA-seq", "Smart-seq2/FACS", "Spatial"]:
                priority = 0 if tech == "10x/scRNA-seq" else 1
                status = "✓ scRNA-seq"
            elif tech == "Unknown":
                # 不知道类型的，标记为待确认
                priority = 2
                status = "? Unknown"
            else:
                # 芯片/Bulk → 降权但保留
                priority = 3
                status = f"~ {tech} (low priority)"
            
            # 构建数据源
            if gse_id in self.geo_ids_counter:
                paper_sources = self.geo_to_papers.get(gse_id, [])
                citation_count = self.geo_ids_counter[gse_id]
                priority = min(priority, 0 if citation_count >= 2 else 1)
            else:
                paper_sources = []
                citation_count = 0
            
            item = {
                "gse_id": gse_id,
                "priority": priority,
                "source_pmids": list(set(p["pmid"] for p in paper_sources)),
                "source_papers": [{"pmid": p["pmid"], "title": p["title"][:100]} for p in paper_sources],
                "disease": meta.get("disease"),
                "technology": tech,
                "n_samples": meta.get("n_samples", 0),
                "organism": meta.get("organism"),
                "title": meta.get("title", ""),
                "citation_count": citation_count,
                "platform": meta.get("platform", ""),
            }
            
            validated.append(item)
            self.log(f"  {status} {gse_id} ({meta.get('disease', 'unknown')}) [{citation_count} citations]")
            time.sleep(0.2)
        
        return validated
    
    # ── Step 6: 生成优先级列表 ──────────────────
    
    def build_priority_list(self, validated_datasets: List[Dict]) -> List[Dict]:
        """生成 P0/P1/P2 优先级数据集列表"""
        
        # 按 priority + citation_count + n_samples 排序
        def sort_key(d):
            return (d["priority"], -d["citation_count"], -d.get("n_samples", 0))
        
        validated_datasets.sort(key=sort_key)
        
        # 分配 P0/P1/P2 标签
        p0_list = []
        p1_list = []
        p2_list = []
        
        for d in validated_datasets:
            if d["citation_count"] >= 3:
                d["tier"] = "P0"
                p0_list.append(d)
            elif d["citation_count"] >= 1:
                d["tier"] = "P1"
                p1_list.append(d)
            elif d["technology"] in ["10x/scRNA-seq", "Smart-seq2/FACS"]:
                d["tier"] = "P2"
                p2_list.append(d)
            else:
                d["tier"] = "P2-FALLBACK"
                p2_list.append(d)
        
        # 优先保留P0/P1（文献引用过的）
        priority_order = p0_list + p1_list + p2_list
        
        self.recommended_datasets = priority_order
        
        self.log(f"\nPriority Dataset List:")
        self.log(f"  P0 (highly cited): {len(p0_list)}")
        self.log(f"  P1 (cited):        {len(p1_list)}")
        self.log(f"  P2 (relevant):     {len(p2_list)}")
        self.log(f"  Total:             {len(priority_order)}")
        
        return priority_order
    
    # ── Step 7: 收集支持文献 ─────────────────────
    
    def build_supporting_papers(self) -> List[Dict]:
        """整理支持文献列表"""
        papers = []
        seen_pmids = set()
        
        for gse_id, paper_list in self.geo_to_papers.items():
            for p in paper_list:
                if p["pmid"] in seen_pmids:
                    continue
                seen_pmids.add(p["pmid"])
                
                paper_detail = self.all_papers.get(p["pmid"], {})
                papers.append({
                    "pmid": p["pmid"],
                    "title": paper_detail.get("title", p["title"]),
                    "journal": paper_detail.get("journal", ""),
                    "gse_ids": list(set(
                        g for g, pl in self.geo_to_papers.items()
                        if p["pmid"] in [x["pmid"] for x in pl]
                    )),
                })
        
        papers.sort(key=lambda x: -len(x.get("gse_ids", [])))
        return papers[:30]  # 最多30篇
    
    # ── Step 8: 保存输出 ─────────────────────────
    
    def save_output(self) -> str:
        """保存 phase2_gene_set.json"""
        
        output = {
            "schema_version": "2.0",
            "generated_at": datetime.now().isoformat(),
            "hypothesis": self.hypothesis_text,
            "gene_sets": self.gene_sets,
            "supporting_papers": self.build_supporting_papers(),
            "recommended_datasets": self.recommended_datasets,
            "statistics": {
                "total_papers_found": len(self.all_papers),
                "total_geo_ids_extracted": len(self.geo_ids_counter),
                "p0_datasets": len([d for d in self.recommended_datasets if d.get("tier") == "P0"]),
                "p1_datasets": len([d for d in self.recommended_datasets if d.get("tier") == "P1"]),
                "p2_datasets": len([d for d in self.recommended_datasets if d.get("tier") in ("P2", "P2-FALLBACK")]),
                "scRNA_datasets": len([d for d in self.recommended_datasets if d.get("technology") in ("10x/scRNA-seq", "Smart-seq2/FACS")]),
                "failed_pmids": len(self.failed_pmids),
                "failed_geo": len(self.failed_geo),
            },
            "search_queries": self.search_queries,
            "metadata": {
                "best_practices_checklist": {
                    "literature_search": True,
                    "geo_id_extraction": True,
                    "priority_ranking": True,
                    "data_type_validation": True,
                    "multi_source_correlation": True,
                }
            }
        }
        
        output_path = self.output_dir / "phase2_gene_set.json"
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(output, f, ensure_ascii=False, indent=2)
        
        self.log(f"\n✅ Phase 2 output saved: {output_path}")
        return str(output_path)
    
    # ── 主流程 ──────────────────────────────────
    
    def run(self) -> str:
        """执行完整的 Phase 2 流程"""
        self.log("=" * 60)
        self.log("PHASE 2: Literature Research & Gene Set Builder")
        self.log(f"Hypothesis: {self.hypothesis_text}")
        self.log("=" * 60)
        
        # Step 1: 生成搜索词
        queries = self.generate_search_queries()
        
        # Step 2: 搜索PubMed
        pmids = self.search_literature(queries)
        
        # Step 3: 从论文中提取GEO ID
        if pmids:
            self.extract_geo_ids_from_papers(pmids)
        
        # Step 4: GEO搜索兜底
        geo_ids = self.geo_search_fallback(queries)
        
        # 合并：文献来源优先，GEO搜索补充
        all_gse_ids = list(set(
            list(self.geo_ids_counter.keys()) + geo_ids
        ))
        
        # Step 5: 验证数据集
        validated = self.validate_and_enrich_datasets(all_gse_ids)
        
        # Step 6: 优先级排序
        priority_list = self.build_priority_list(validated)
        
        # Step 7: 保存输出
        output_path = self.save_output()
        
        # 打印摘要
        self.log("\n" + "=" * 60)
        self.log("PHASE 2 SUMMARY")
        self.log(f"  Papers analyzed:     {len(self.all_papers)}")
        self.log(f"  GEO IDs extracted:  {len(self.geo_ids_counter)}")
        self.log(f"  Recommended datasets: {len(self.recommended_datasets)}")
        self.log(f"  P0/P1/P2:            "
                 f"{len([d for d in self.recommended_datasets if d.get('tier')=='P0'])}/"
                 f"{len([d for d in self.recommended_datasets if d.get('tier')=='P1'])}/"
                 f"{len([d for d in self.recommended_datasets if d.get('tier') in ('P2','P2-FALLBACK')])}")
        self.log("=" * 60)
        
        return output_path


# ─────────────────────────────────────────────
# 入口
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Phase 2: Literature Mining & Gene Set Builder")
    parser.add_argument("hypothesis", nargs="?", help="Hypothesis text")
    parser.add_argument("--hypothesis", dest="h2", help="Hypothesis (alternative flag)")
    parser.add_argument("--output", "-o", default=None, help="Output directory")
    parser.add_argument("--hypothesis-id", help="Custom hypothesis ID (e.g. h_003)")
    args = parser.parse_args()
    
    hypothesis = args.hypothesis or args.h2
    if not hypothesis:
        print("Usage: python phase2_literature_miner.py \"your hypothesis here\" [--output data/]")
        sys.exit(1)
    
    setup_proxy()
    
    script_base = os.path.dirname(os.path.abspath(__file__))
    output_dir = args.output or os.path.join(script_base, "..", "data")
    if args.hypothesis_id:
        output_dir = os.path.join(script_base, "..", "data", args.hypothesis_id)
    
    miner = LiteratureMiner(hypothesis_text=hypothesis, output_dir=output_dir)
    output_path = miner.run()
    print(f"\n📦 Output: {output_path}")


if __name__ == "__main__":
    main()
