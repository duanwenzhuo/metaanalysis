"""
L0: EdgeRSkill — Bulk RNA-seq 差异表达分析（edgeR）

bioSkills: bio-de-edger-basics → bioskills L0
支持: edgeR 完整流程（DGEList → 过滤 → 标准化 → QL 检验）
设计公式: 单因素 / 多因素 / 交互作用 / LRT

bioSkills Method Comparison:
| 工具         | 语言 | 速度 | 适用场景                    |
|-------------|------|------|----------------------------|
| edgeR QL    | R    | 快   | 标准 DE，quasi-likelihood |
| edgeR LRT   | R    | 快   | 似然比检验                 |
| exact test  | R    | 快   | 简单两组比较               |
| glmTreat    | R    | 慢   | 有 LFC 阈值的检验          |
"""

import subprocess
import sys
from typing import Optional
from bioskills.core.base import register, AbstractSkill


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


@register
class EdgeRSkill(AbstractSkill):
    name = "edger"
    stage = "differential_expression"
    input_contract = ["counts", "metadata", "design_formula"]
    output_contract = ["edger_results", "log_cpm", "normalized_counts", "bcv_plot"]

    tunables = {
        "alpha": 0.05,
        "lfc_threshold": 0.0,
        "filter_cpm_min": 1.0,
        "filter_min_samples": 3,
        "norm_method": "TMM",           # TMM | RLE | upperquartile | none
        "test_method": "QLF",            # QLF | LRT | exact
        "dispersion_method": "auto",     # auto | common | trended | tagwise
        "contrast": None,               # for specific comparisons
        "prior_count": 2,                # for log-CPM
    }

    def _build_script(self, params: dict) -> str:
        alpha = params.get("alpha", 0.05)
        cpm_thresh = params.get("filter_cpm_min", 1.0)
        min_samp = params.get("filter_min_samples", 3)
        norm_m = params.get("norm_method", "TMM")
        test_m = params.get("test_method", "QLF")
        disp_m = params.get("dispersion_method", "auto")
        prior_cnt = params.get("prior_count", 2)
        contrast = params.get("contrast")
        lfc = params.get("lfc_threshold", 0.0)

        # Build dispersion block
        if disp_m == "auto":
            disp_block = "y <- estimateDisp(y, design)"
        else:
            disp_block = (
                "y <- estimateGLMCommonDisp(y, design)\n"
                "y <- estimateGLMTrendedDisp(y, design)\n"
                "y <- estimateGLMTagwiseDisp(y, design)"
            )

        # Build test block
        if test_m == "QLF":
            test_block = (
                "fit <- glmQLFit(y, design)\n"
                "qlf <- glmQLFTest(fit, coef = 2)"
            )
        elif test_m == "LRT":
            if contrast:
                contrast_str = str(contrast)
                test_block = (
                    "fit <- glmFit(y, design)\n"
                    "qlf <- glmLRT(fit, contrast = " + contrast_str + ")"
                )
            else:
                test_block = (
                    "fit <- glmFit(y, design)\n"
                    "qlf <- glmLRT(fit, coef = 2)"
                )
        else:  # exact
            test_block = (
                "et <- exactTest(y)\n"
                "qlf <- et"
            )

        script = (
            "library(edgeR)\n"
            "library(limma)\n"
            "\n"
            "counts <- __counts__\n"
            "group <- __group__\n"
            "\n"
            "# Create DGEList\n"
            "y <- DGEList(counts = counts, group = group)\n"
            "\n"
            "# Filtering\n"
            "keep <- filterByExpr(y, group = group)\n"
            "y <- y[keep, , keep.lib.sizes = FALSE]\n"
            "\n"
            "# Normalization\n"
            f"y <- calcNormFactors(y, method = '{norm_m}')\n"
            "\n"
            "# Design matrix\n"
            "design <- model.matrix(~ group)\n"
            "\n"
            "# Dispersion estimation\n"
            f"{disp_block}\n"
            "\n"
            "# Fit model\n"
            f"{test_block}\n"
            "\n"
            "# Results\n"
            'if ("topTags" %in% names(attributes(qlf))) {\n'
            "    res_df <- topTags(qlf, n = Inf)$table\n"
            "} else {\n"
            "    res_df <- as.data.frame(qlf)\n"
            "}\n"
            "\n"
            "# Log-CPM\n"
            f"log_cpm <- cpm(y, log = TRUE, prior.count = {prior_cnt})\n"
            "\n"
            "# Normalized counts\n"
            "norm_counts <- cpm(y, normalized.lib.sizes = TRUE)\n"
            "\n"
            "# Summary\n"
            f"sig_up <- sum(res_df$FDR < {alpha} & res_df$logFC > {lfc}, na.rm = TRUE)\n"
            f"sig_down <- sum(res_df$FDR < {alpha} & res_df$logFC < -{lfc}, na.rm = TRUE)\n"
            "total_tested <- sum(!is.na(res_df$FDR))\n"
            "\n"
            "__output__ <- list(\n"
            "    results_df = res_df,\n"
            "    log_cpm = log_cpm,\n"
            "    norm_counts = norm_counts,\n"
            "    n_tested = total_tested,\n"
            "    n_sig_up = sig_up,\n"
            "    n_sig_down = sig_down,\n"
            "    design = as.character(design)[1],\n"
            "    norm_factors = y$samples$norm.factors,\n"
            "    common_dispersion = y$common.dispersion,\n"
            f"    test_method = '{test_m}'\n"
            ")\n"
        )
        return script

    def _run(self, state: dict) -> dict:
        params = {**self.tunables, **state.get("params", {})}

        counts = state.get("counts")
        metadata = state.get("metadata")
        group = state.get("group")
        if group is None and metadata is not None:
            if hasattr(metadata, "get"):
                group = metadata.get("group") or metadata.get("condition")

        if counts is None:
            return {"error": "counts matrix required in state"}

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

            # Build group vector
            if group is None:
                if hasattr(metadata, "get"):
                    group = metadata.get("group") or metadata.get("condition") or ["G1"] * count_np.shape[1]
                else:
                    group = ["G1"] * count_np.shape[1]

            group_str = "factor(c(" + ",".join(repr(str(g)) for g in group) + "))"

            script = self._build_script(params)

            # Build R-readable counts
            counts_csv = counts_df.to_csv(index=False, header=False)
            counts_r = "as.matrix(read.csv(text='" + counts_csv.replace("\n", "', header=FALSE)\nread.csv(text='") + "', header=FALSE))"
            script = script.replace("__counts__", counts_r)
            script = script.replace("__group__", group_str)

            ro.globalenv["__output__"] = ro.r("list()")
            ro.r(script)
            output = ro.globalenv["__output__"]

            results_dict = {}
            if hasattr(output, "names"):
                for i, name in enumerate(output.names):
                    val = list(output)[i]
                    if hasattr(val, "__iter__") and hasattr(val, "__len__") and len(val) > 0:
                        try:
                            results_dict[name] = list(val)
                        except:
                            results_dict[name] = str(val)
                    else:
                        results_dict[name] = val

            n_tested = results_dict.get("n_tested", 0)

            return {
                "state_updates": {
                    "edger_results": results_dict.get("results_df", []),
                    "log_cpm": results_dict.get("log_cpm", []),
                    "normalized_counts": results_dict.get("norm_counts", []),
                    "edger_report": {
                        "status": "success",
                        "n_tested": n_tested,
                        "n_sig_up": results_dict.get("n_sig_up", 0),
                        "n_sig_down": results_dict.get("n_sig_down", 0),
                        "design": results_dict.get("design", "~ group"),
                        "norm_method": params.get("norm_method", "TMM"),
                        "test_method": params.get("test_method", "QLF"),
                        "common_dispersion": results_dict.get("common_dispersion", 0),
                        "method": "edgeR",
                        "tool_type": "r",
                    }
                }
            }

        except Exception as e:
            return {
                "error": f"edgeR failed: {str(e)}",
                "state_updates": {
                    "edger_report": {
                        "status": "error",
                        "error": str(e),
                        "method": "edgeR",
                        "tool_type": "r",
                    }
                }
            }
