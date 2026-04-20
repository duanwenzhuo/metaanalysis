"""
L0: DESeq2Skill — Bulk RNA-seq 差异表达分析

bioSkills: bio-de-deseq2-basics → bioskills L0
支持: DESeq2 完整流程（创建 → 过滤 → 标准化 → DE检验）
设计公式: 单因素/多因素/交互作用/LRT

bioSkills Method Comparison:
| 工具       | 语言 | 速度 | 适用场景                    |
|-----------|------|------|----------------------------|
| DESeq2    | R    | 慢   | 标准 DE，count 数据         |
| apeglm    | R    | 慢   | LFC shrinkage（推荐）      |
| ashr      | R    | 慢   | LFC shrinkage（备选）      |
| normal    | R    | 快   | LFC shrinkage（旧方法）    |
"""

import subprocess
import sys
from typing import Optional
from bioskills.core.base import register, AbstractSkill

# ── rpy2 bridge ──────────────────────────────────────────────────────────────

def _r_install(pkg: str) -> bool:
    """Install R package via BiocManager."""
    result = subprocess.run(
        ["Rscript", "-e",
         f"if(!require('{pkg}', quietly=TRUE)) {{ BiocManager::install('{pkg}') }}"],
        capture_output=True, text=True
    )
    return result.returncode == 0


def _ensure_rpy2():
    try:
        import rpy2.robjects as ro
        return ro
    except ImportError:
        try:
            subprocess.run([sys.executable, "-m", "pip3", "install", "rpy2", "-q"],
                           check=True)
            import rpy2.robjects as ro
            return ro
        except Exception:
            return None


def _try_r(cmd: str) -> tuple:
    """Execute R command via rpy2, returns (success, result_dict)."""
    ro = _ensure_rpy2()
    if ro is None:
        return False, {"error": "rpy2 not available. Install: pip install rpy2"}

    try:
        ro.globalenv["__output__"] = ro.r("list()")
        full_script = f"""
library(DESeq2)
library(apeglm)
{cmd}
"""
        ro.r(full_script)
        output = ro.globalenv["__output__"]
        result = {}
        if hasattr(output, "names"):
            for i, name in enumerate(output.names):
                val = list(output)[i]
                if hasattr(val, "__iter__") and not isinstance(val, (str, bytes)):
                    result[name] = list(val)
                else:
                    result[name] = val
        return True, result
    except Exception as e:
        return False, {"error": str(e)}


# ── Skill ────────────────────────────────────────────────────────────────────

@register
class DESeq2Skill(AbstractSkill):
    name = "deseq2"
    stage = "differential_expression"
    input_contract = ["counts", "metadata", "design_formula"]
    output_contract = ["deseq2_results", "normalized_counts", "vsd", "rld"]

    tunables = {
        "alpha": 0.05,
        "lfc_threshold": 0.0,
        "filter_min_count": 10,
        "filter_min_samples": 3,
        "shrinkage_method": "apeglm",   # apeglm | ashr | normal
        "test_type": "Wald",            # Wald | LRT
        "reduced_formula": None,        # for LRT
        "reference_level": None,
        "parallel": False,
        "blind": False,
    }

    def _build_dds_script(self, params: dict) -> str:
        alpha = params.get("alpha", 0.05)
        lfc = params.get("lfc_threshold", 0.0)
        filt_count = params.get("filter_min_count", 10)
        filt_samp = params.get("filter_min_samples", 3)
        test = params.get("test_type", "Wald")
        ref = params.get("reference_level")
        shrink = params.get("shrinkage_method", "apeglm")
        reduced = params.get("reduced_formula")
        parallel = params.get("parallel", False)

        script = f"""
counts <- __counts__
metadata <- __metadata__
design <- __design__

dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = metadata,
                               design = design)

# Pre-filter
keep <- rowSums(counts(dds)) >= {filt_count}
dds <- dds[keep, ]

# Reference level
{ f"dds$condition <- relevel(dds$condition, ref = '{ref}')" if ref else "# no ref set" }

# Run DESeq2
{ f"dds <- DESeq(dds, test = 'LRT', reduced = {reduced}, parallel = {parallel})"
  if test == "LRT" else f"dds <- DESeq(dds, parallel = {parallel})" }

# Results
res <- results(dds, alpha = {alpha}, lfcThreshold = {lfc})

# LFC shrinkage
coef_name <- resultsNames(dds)[2]
resLFC <- lfcShrink(dds, coef = coef_name, type = '{shrink}')

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# VST
vsd <- vst(dds, blind = {params.get('blind', False)})

# Rlog
rld <- rlog(dds, blind = {params.get('blind', False)})

# Summary
sig_up <- sum(resLFC$padj < {alpha} & resLFC$log2FoldChange > {lfc}, na.rm = TRUE)
sig_down <- sum(resLFC$padj < {alpha} & resLFC$log2FoldChange < -{lfc}, na.rm = TRUE)
total_tested <- sum(!is.na(resLFC$padj))

__output__ <- list(
    results_df = as.data.frame(resLFC),
    norm_counts = norm_counts,
    vsd_mat = assay(vsd),
    rld_mat = assay(rld),
    design_formula = as.character(design)[2],
    n_tested = total_tested,
    n_sig_up = sig_up,
    n_sig_down = sig_down,
    n_genes_remain = nrows(dds),
    shrinkage_method = '{shrink}',
    size_factors = sizeFactors(dds)
)
"""
        return script

    def _run(self, state: dict) -> dict:
        params = {**self.tunables, **state.get("params", {})}

        counts = state.get("counts")
        metadata = state.get("metadata")
        design_formula = state.get("design_formula", "~ condition")

        if counts is None:
            return {"error": "counts matrix required in state"}

        # Convert pandas/numpy to R matrix
        ro = _ensure_rpy2()
        if ro is None:
            return {"error": "rpy2 not available: pip install rpy2", "state_updates": {}}

        try:
            import pandas as pd
            import numpy as np

            if hasattr(counts, "to_numpy"):
                count_np = counts.to_numpy()
            elif hasattr(counts, "values"):
                count_np = counts.values
            else:
                count_np = np.array(counts)

            counts_df = pd.DataFrame(count_np)

            if hasattr(metadata, "to_dict"):
                meta_df = pd.DataFrame(metadata)
            elif hasattr(metadata, "values"):
                meta_df = pd.DataFrame(metadata)
            else:
                meta_df = pd.DataFrame(metadata)

            script = self._build_dds_script(params)
            script = script.replace("__counts__", "as.matrix(read.csv(text='"
                                     + counts_df.to_csv(index=False, header=False)
                                     + "', header=FALSE))")
            script = f"read.csv_all <- function(file) {{ df <- read.csv(file, header=FALSE); as.matrix(df) }}\n" + script.replace(
                "__counts__", "read.csv_all(text='"
                + counts_df.to_csv(index=False, header=False).replace("\n", "', header=FALSE)\nread.csv_all(text='")
                + "', header=FALSE))")

            # Build metadata inline
            meta_r = "data.frame("
            for col in meta_df.columns:
                vals = meta_df[col].tolist()
                if isinstance(vals[0], str):
                    meta_r += f'{col}=factor(c({",".join(repr(str(v)) for v in vals)})), '
                elif isinstance(vals[0], bool):
                    meta_r += f'{col}=c({",".join(str(v) for v in vals)}), '
                else:
                    meta_r += f'{col}=c({",".join(str(v) for v in vals)}), '
            meta_r = meta_r.rstrip(", ") + ")"
            script = script.replace("__metadata__", meta_r)
            script = script.replace("__design__", design_formula.replace("~", "~"))

            ro.globalenv["__output__"] = ro.r("list()")
            ro.r(script)
            output = ro.globalenv["__output__"]

            results_dict = {}
            if hasattr(output, "names"):
                for i, name in enumerate(output.names):
                    val = list(output)[i]
                    if hasattr(val, "__iter__") and hasattr(val, "__len__"):
                        if len(val) > 0:
                            try:
                                results_dict[name] = list(val)
                            except:
                                results_dict[name] = str(val)
                    else:
                        results_dict[name] = val

            n_tested = results_dict.get("n_tested", 0)
            n_sig = results_dict.get("n_sig_up", 0) + results_dict.get("n_sig_down", 0)

            # Update state
            if "adata" in state:
                adata = state["adata"]
                if "X" in adata.layers:
                    adata.layers["deseq2_normalized"] = results_dict.get("norm_counts", adata.X)

            return {
                "state_updates": {
                    "deseq2_results": results_dict.get("results_df", []),
                    "normalized_counts": results_dict.get("norm_counts", []),
                    "vsd": results_dict.get("vsd_mat", []),
                    "rld": results_dict.get("rld_mat", []),
                    "deseq2_report": {
                        "status": "success",
                        "n_tested": n_tested,
                        "n_sig_up": results_dict.get("n_sig_up", 0),
                        "n_sig_down": results_dict.get("n_sig_down", 0),
                        "design": results_dict.get("design_formula", design_formula),
                        "shrinkage": results_dict.get("shrinkage_method", "apeglm"),
                        "method": "DESeq2",
                        "tool_type": "r",
                    }
                }
            }

        except Exception as e:
            return {
                "error": f"DESeq2 failed: {str(e)}",
                "state_updates": {
                    "deseq2_report": {
                        "status": "error",
                        "error": str(e),
                        "method": "DESeq2",
                        "tool_type": "r",
                    }
                }
            }
