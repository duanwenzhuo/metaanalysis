# 2026-04-21 — MetaAnalysisi 整理为独立项目

## 任务
将 `/Users/tang/.qclaw/workspace/bio-analysis/` 整理为独立项目，迁移到 `/Users/tang/.qclaw/workspace-taizi/metaanalysis/`。

## 整理原则
1. **去掉冗余旧脚本**（analyze_gse120575_v*.py, batch_analyze_v3.py 等已废弃）
2. **保留三层架构核心**：workflow/LangGraph → bioskills L0/L1/L2 → core工具库
3. **去掉旧硬编码路径**，改为相对路径/动态路径
4. **数据外置**：data/ 目录留在原处，不复制

## 最终目录结构

```
/Users/tang/.qclaw/workspace-taizi/metaanalysis/
├── README.md                    # 项目说明
├── main/                        # CLI 入口
│   ├── __init__.py
│   ├── __main__.py
│   └── cli.py
├── workflow/                    # LangGraph 状态机
│   ├── state.py                 # TypedDict BioAnalysisState（50字段）
│   ├── graph.py                 # StateGraph + 条件边 + END跳转
│   ├── registry.py              # 断点续跑 StateRegistry（原子写入）
│   └── nodes/
│       ├── hypothesis_parser.py
│       ├── literature_miner.py
│       ├── geo_retriever.py
│       ├── effect_size_analyst.py   # bioskills SkillPipelineEngine
│       ├── meta_analyst.py
│       └── report_curator.py
├── bioskills/                   # 三层自进化技能体系
│   ├── __init__.py
│   ├── core/                    # AbstractSkill + SkillRegistry + Contract
│   │   ├── __init__.py
│   │   ├── base.py              # 508行，@register装饰器，自动发现
│   │   ├── contract.py          # ContractValidator 契约验证
│   │   └── pipeline.py          # PipelineResult
│   ├── execution/               # 统一执行引擎
│   │   ├── __init__.py
│   │   └── pipeline_engine.py   # SkillPipelineEngine（含data_io自动加载）
│   ├── l0/                      # 原子技能（~42个，在原 bio-analysis bioskills/）
│   ├── l1/                      # Pipeline模板
│   │   └── cell_state_pipeline.py
│   ├── l2/                      # 自愈层
│   │   └── self_heal.py        # GapAnalyzer + SkillSynthesizer
│   └── knowledge/               # GeneSetDB（未复制，依赖原目录）
├── scripts/                     # Phase 2/5/6 独立脚本
│   ├── phase2_literature_miner.py
│   ├── phase5_meta_orchestrator.py
│   ├── phase6_report_generator.py
│   └── geo_scseq_downloader.py
├── core/                        # 工具库（~550行）
│   ├── config.py, constants.py, types.py
│   ├── utils.py, exceptions.py
│   └── __init__.py
└── tests/                       # 测试（从原目录复制）

## 关键修复
1. `bioskills/core/base.py` 硬编码路径 → `Path(__file__).parent.parent` 动态路径
2. `bioskills/l2/self_heal.py` 硬编码路径 → 同上
3. `scripts/phase5_meta_orchestrator.py` 硬编码前缀 → 移除
4. `scripts/phase6_report_generator.py` 同上
5. `workflow/registry.py` 同上

## 当前阻塞问题（未解决）
- **data_io ContractError**：pipeline_engine 执行 data_io 时，execute() 的 _validate_output 检查 output_contract，但 data_io 的 _run() 成功返回 adata，execute() 层面却报 "Output contract violated"
  - 原因：data_io 加载了真实 adata (737280×33694)，但 _validate_output 可能没找到返回的 adata
  - 待查：execute() 调用 _run() 后的 result dict 里是否包含 adata
- **Cohen's d = 0**：skill chain 没拿到 adata，导致所有技能都报 ContractError

## 验证结果
```bash
cd /Users/tang/.qclaw/workspace-taizi/metaanalysis
python3 -c "from bioskills import SkillRegistry, SkillPipelineEngine; print('OK')"
# ✅ bioskills core import OK
```
