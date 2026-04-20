# Execution Layer — bioskills 执行引擎

## 概述

`execution/` 是 bioskills 的管道执行基础设施，提供统一的技能调度引擎和预定义管道模板。

## 核心组件

### SkillPipelineEngine

管道引擎类，统一调用所有已注册的 bioskills 技能。

**方法：**
| 方法 | 说明 |
|------|------|
| `execute_pipeline(state, pipeline)` | 执行完整管道 |
| `get_pipeline_info()` | 获取可用管道列表 |

### execute_pipeline(state, pipeline)

快捷函数，执行指定管道。

```python
from bioskills.execution import execute_pipeline

state = {
    "adata": adata,
    "params": {"resolution": 0.5, "species": "human"}
}
result = execute_pipeline(state, pipeline="cell_state_full")
```

### pipeline_info()

返回所有可用管道定义。

### PipelineResult

管道执行结果 dataclass，包含：
- `status`: `"success"` / `"failed"` / `"partial"`
- `logs`: 执行日志列表
- `artifacts`: 输出产物字典
- `adata`: 处理后的 AnnData（如果适用）
- `duration_seconds`: 执行耗时

## 预定义管道

| 管道名 | 描述 |
|--------|------|
| `cell_state_full` | 完整细胞状态分析（QC→PCA→Neighbors→Leiden→UMAP→MarkerScore） |
| `differential_expression` | 差异表达分析 |
| `trajectory_analysis` | 轨迹分析 |

## 目录结构

```
execution/
├── SKILL.md              ← 本文件
├── __init__.py
└── pipeline_engine.py   ← 核心引擎实现
```

## 使用场景

- **批量执行**：不想手动串联技能时，直接调用 `execute_pipeline`
- **管道组合**：查看 `pipeline_info()` 了解可用管道，修改参数
- **结果包装**：所有执行结果统一为 `PipelineResult`，便于日志追踪