"""
L0: GOEnrichmentSkill — GO/KEGG 通路富集分析

bioSkills: bio-pathway-analysis-go-enrichment + kegg-pathways → bioskills L0
支持: GO-BP / GO-MF / GO-CC / KEGG / Reactome
工具: gseapy / goatools / enrichr

bioSkills Method Comparison:
| 工具        | 语言 | 数据源          | 适用场景                |
|------------|------|----------------|------------------------|
| gseapy     | Py   | Enrichr API    | 快速在线查询            |
| goatools   | Py   | 本地 OBO        | 离线/批量               |
| clusterProfiler | R  | 本地            | R标准                   |
"""

import os
from typing import Optional, List, Dict
from bioskills.core.base import AbstractSkill


class GOEnrichmentSkill(AbstractSkill):
    name = "go_enrichment"
    stage = "enrichment"
    input_contract = ["gene_list"]
    output_contract = ["go_results", "kegg_results", "enrichment_plot"]

    tunables = {
        "databases": ["GO_Biological_Process_2023", "KEGG_2021_Human"],
        "organism": "human",               # human | mouse | yeast | fly | fish | worm
        "threshold": 0.05,
        "cutoff": 0.05,                     # adjusted p-value cutoff
        "top_n": 20,                        # show top N pathways
        "background": None,                 # background gene count
        "backend": "gseapy",                # gseapy | goatools
        "output_path": None,
        "format": "png",
        "show_plots": False,
    }

    # Enrichr organism mapping
    ENRICHR_LIBS = {
        "human": ["GO_Biological_Process_2023", "GO_Molecular_Function_2023",
                   "GO_Cellular_Component_2023", "KEGG_2021_Human",
                   "Reactome_2022", "MSigDB_Hallmark_2020"],
        "mouse": ["GO_Biological_Process_2023", "GO_Molecular_Function_2023",
                   "GO_Cellular_Component_2023", "KEGG_2019_Mouse",
                   "Reactome_2022"],
    }

    def _run(self, state: dict) -> dict:
        params = {**self.tunables, **state.get("params", {})}

        gene_list = state.get("sig_genes") or state.get("gene_list") or state.get("up_genes")
        if gene_list is None:
            return {"error": "gene_list / sig_genes required in state"}

        if isinstance(gene_list, str):
            gene_list = [g.strip() for g in gene_list.split(",")]

        backend = params.get("backend", "gseapy")

        if backend == "gseapy":
            return self._run_gseapy(gene_list, params)
        elif backend == "goatools":
            return self._run_goatools(gene_list, params)
        else:
            return {"error": f"Unknown backend: {backend}"}

    def _run_gseapy(self, gene_list, params):
        try:
            import gseapy as gp
        except ImportError:
            return {"error": "gseapy required: pip install gseapy"}

        organism = params.get("organism", "human")
        databases = params.get("databases")
        if not databases:
            databases = self.ENRICHR_LIBS.get(organism, self.ENRICHR_LIBS["human"])[:3]

        outdir = params.get("output_path") or "enrichr_output"
        cutoff = params.get("cutoff", 0.05)

        all_results = {}
        all_enriched = {}

        for db in databases:
            try:
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=db,
                    organism=organism,
                    outdir=os.path.join(outdir, db.replace("/", "_")),
                    cutoff=cutoff,
                    no_plot=not params.get("show_plots", False),
                )

                res_df = enr.results
                if res_df is not None and len(res_df) > 0:
                    # Filter significant
                    if "Adjusted P-value" in res_df.columns:
                        sig = res_df[res_df["Adjusted P-value"].astype(float) < cutoff]
                    else:
                        sig = res_df

                    all_results[db] = res_df.to_dict(orient="records")
                    enriched = []
                    for _, row in sig.head(params.get("top_n", 20)).iterrows():
                        enriched.append({
                            "pathway": row.get("Term", ""),
                            "pval": float(row.get("P-value", 1)),
                            "fdr": float(row.get("Adjusted P-value", 1)),
                            "overlap": row.get("Overlap", ""),
                            "genes": row.get("Genes", ""),
                            "db": db,
                        })
                    all_enriched[db] = enriched
            except Exception as e:
                all_enriched[db] = [{"error": str(e)}]

        # Aggregate
        total_enriched = []
        for db, items in all_enriched.items():
            total_enriched.extend(items)
        total_enriched.sort(key=lambda x: x.get("fdr", 1))

        return {
            "state_updates": {
                "go_results": all_results,
                "enriched_pathways": total_enriched,
                "go_report": {
                    "status": "success",
                    "databases": databases,
                    "n_total_enriched": len([x for x in total_enriched if "error" not in x]),
                    "organism": organism,
                    "backend": "gseapy",
                    "cutoff": cutoff,
                }
            }
        }

    def _run_goatools(self, gene_list, params):
        try:
            from goatools import GOEnrichmentStudy
            from goatools.obo_parser import GODag
        except ImportError:
            return {"error": "goatools required: pip install goatools"}

        organism = params.get("organism", "human")
        cutoff = params.get("cutoff", 0.05)

        # Load GO DAG
        try:
            obodag = GODag("go-basic.obo")
        except Exception:
            return {"error": "go-basic.obo not found. Download: wget http://geneontology.org/ontology/go-basic.obo"}

        # Background genes
        background = params.get("background")
        if background is None:
            background = gene_list  # minimal, should be all tested genes

        goeaobj = GOEnrichmentStudy(
            background,
            obodag,
            methods=["fdr_bh"],
        )

        results = goeaobj.run_study(gene_list)

        # Filter significant
        sig_results = [r for r in results if r.p_fdr_bh < cutoff]

        enriched = []
        for r in sig_results[:params.get("top_n", 20)]:
            enriched.append({
                "go_id": r.GO,
                "name": r.name,
                "namespace": r.NS,
                "pval": r.p_uncorrected,
                "fdr": r.p_fdr_bh,
                "ratio_in_study": f"{r.ratio_in_study[0]}/{r.ratio_in_study[1]}",
                "ratio_in_pop": f"{r.ratio_in_pop[0]}/{r.ratio_in_pop[1]}",
                "study_items": list(r.study_items),
            })

        return {
            "state_updates": {
                "go_results": enriched,
                "enriched_pathways": enriched,
                "go_report": {
                    "status": "success",
                    "n_enriched": len(enriched),
                    "organism": organism,
                    "backend": "goatools",
                    "cutoff": cutoff,
                }
            }
        }
