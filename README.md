# MetaAnalysisi — 文献驱动的单细胞元分析系统

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

从 PubMed 文献提取研究假设 → 自动检索 GEO 单细胞数据集 → 执行标准化分析流水线 → 元分析合并效应量 → 生成报告。

**核心特点：**
- 📖 **假设驱动**：输入自然语言假设，自动拆解为可检验命题
- 🔬 **无监督分析**：无需标注，GSVA + Marker Score 自动细胞类型注释
- 🛠️ **42+ 生物信息学技能**：模块化、可组合、契约验证
- 📊 **Meta 分析**：Cohen's d + 随机效应模型，跨数据集效应量合并


---

## 🏗️ 系统架构

```
metaanalysis/
├── main/                        # CLI 入口
│   └── cli.py                   # run / list / resume / status
├── workflow/                    # LangGraph 状态机（6阶段）
│   ├── state.py                 # BioAnalysisState (TypedDict, 50字段)
│   ├── graph.py                 # StateGraph + 条件边 + MemorySaver
│   ├── registry.py              # 断点续跑 StateRegistry
│   └── nodes/
│       ├── hypothesis_parser.py      # 假设解析
│       ├── literature_miner.py      # PubMed 文献挖掘
│       ├── geo_retriever.py          # GEO 数据集检索
│       ├── effect_size_analyst.py    # 效应量分析
│       ├── meta_analyst.py           # 元分析
│       └── report_curator.py         # 报告生成
├── bioskills/                   # 三层自进化技能体系
│   ├── core/
│   │   ├── base.py             # AbstractSkill + @register + SkillRegistry
│   │   └── contract.py          # ContractValidator 契约验证
│   ├── execution/
│   │   └── pipeline_engine.py  # SkillPipelineEngine（顺序/并行/断点续跑）
│   ├── l0/                     # 42+ 原子技能（按 Stage 分类）
│   ├── l1/                     # Pipeline 模板
│   │   └── cell_state_pipeline.py
│   ├── l2/                     # 自愈层
│   │   └── self_heal.py        # GapAnalyzer + SkillSynthesizer
│   └── knowledge/
│       ├── gene_set_db.py     # GeneSetDB (70 cell types × 55 gene sets)
│       ├── go_enricher.py     # GO 富集
│       └── msigdb_client.py   # MSigDB
├── scripts/
│   ├── phase2_literature_miner.py
│   ├── phase5_meta_orchestrator.py
│   ├── phase6_report_generator.py
│   └── geo_scseq_downloader.py
└── core/                       # 工具库
    ├── config.py, constants.py
    ├── types.py, utils.py, exceptions.py
```

### LangGraph 工作流

```
hypothesis
    ↓
hypothesis_parser → literature_miner → geo_retriever
                                            ↓
                              effect_size_analyst (SkillPipelineEngine)
                                   ↓
                              meta_analyst → report_curator → END
```

---

## 🛠️ bioskills L0 技能（42+ 个）

所有技能继承 `AbstractSkill`，带 **input_contract / output_contract 契约验证**。

| Stage | 技能 | 说明 |
|-------|------|------|
| **preprocessing** | `data_io` | 自动识别 10X/H5AD/Loom/CSV 加载 |
| | `qc` | QC 过滤（min_genes/min_cells/MT%过滤） |
| | `normalize` | Total Count / CPM / DESeq2 / TMM 标准化 |
| | `hvg` | Seurat v3 / M3Drop / Brennecke 高变异基因 |
| | `scale` | Z-score / Pearson / None 缩放 |
| **dimension** | `pca` | PCA 降维（50 PCs 默认） |
| | `neighbors` | KNN 邻居图（SNN-louvain） |
| | `umap` | UMAP 2D 可视化 |
| **cluster** | `leiden` | Leiden 聚类（多分辨率） |
| | `louvain` | Louvain 聚类 |
| **annotation** | `marker_score` | Marker Gene Score 细胞类型注释 |
| | `cell_annotation_ref` | 参考数据集注释 |
| **enrichment** | `gsva` | GSVA 基因集变异分析 |
| | `gsea` | GSEA 基因集富集分析 |
| | `go_enrichment` | GO/KEGG 富集（gseapy） |
| **statistics** | `cohens_d` | Cohen's d 效应量计算 |
| | `mann_whitney` | Mann-Whitney U 检验 |
| **differential_expression** | `deseq2` | DESeq2 差异表达（R） |
| | `edger` | edgeR 差异表达（R） |
| | `de_results` | 差异表达结果汇总 |
| **visualization** | `volcano` | 火山图（matplotlib + plotly） |
| | `de_visualization` | 热图/MA图/VennDiagram |
| **integration** | `batch_integration` | Harmony/scVI/Seurat 批次整合 |
| | `batch_correction` | ComBat/limma removeBatchEffect |
| **communication** | `cell_chat` | CellChat/NicheNet/LIANA |
| **trajectory** | `trajectory_inference` | Monocle3/scVelo/PAGA/Slingshot |
| **qc_advanced** | `doublet_detection` | Scrublet/DoubletFinder/scDblFinder |
| | `cell_cycle` | CellCycleScoring 细胞周期 |
| **splicing** | *(更多技能...) | |



---

## 📊 数据结构

```
data/                          # 用户本地数据（gitignore）
├── hypotheses/                # 每个假设的完整输出
│   ├── h_001_t_cell_exhaustion/
│   │   ├── manifest.json      # 阶段状态记录
│   │   ├── phase2/            # 文献挖掘结果
│   │   ├── phase3/            # GEO 下载记录
│   │   ├── phase4/            # 效应量分析结果
│   │   ├── phase5/            # Meta 分析结果
│   │   └── phase6/            # 最终报告
│   └── h_002_xxx/
└── references/               # 参考数据集
```

**manifest.json 格式：**
```json
{
  "hypothesis_id": "h_001_t_cell_exhaustion",
  "created_at": "2026-04-17T...",
  "status": "phase5_complete",
  "datasets": [...],
  "phases": {
    "phase2": {"status": "done", "papers": 47},
    "phase4": {"status": "done", "n_datasets": 5},
    "phase5": {"status": "done", "combined_cohens_d": 1.78}
  }
}
```

---

## 🔧 配置

编辑 `core/config.py` 或设置环境变量：

```bash
export BIOANALYSIS_DATA_DIR="./data"
export BIOANALYSIS_CACHE_DIR="./cache"
export BIOANALYSIS_LOG_LEVEL="INFO"
```

---

## 📝 许可

MIT License — 可自由使用、修改、分发。

