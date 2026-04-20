# L2 Synthesized Skills — _synthesized/

## 概述

`synthesized` 层是 bioskills 的高层技能合成区——由多个原子技能组合而成的复杂技能。

## 技能清单

| 文件 | 技能名 | 类名 | 描述 |
|------|--------|------|------|
| `custom_gsva_scores.py` | `custom_gsva_scores` | `GsvaScoresSkill` | 自定义基因集 GSVA 评分 |

### GsvaScoresSkill

**用途：** 灵活的基因集变异分析（GSVA），支持任意自定义基因集。

**输入合同：** `["adata", "params"]`
**输出合同：** `["gsva_scores", "adata"]`
**阶段：** `ENRICHMENT`

**参数：**
| 参数 | 类型 | 默认 | 说明 |
|------|------|------|------|
| `gene_sets` | dict | `{}` | 字典，key=基因集名，val=基因列表 |
| `sample_norm` | str | `"rank"` | 样本归一化方式 |
| `min_size` | int | `5` | 最小基因数 |
| `max_size` | int | `500` | 最大基因数 |

**示例：**
```python
state = {
    "adata": adata,
    "params": {
        "gene_sets": {
            "T_cell_signature": ["CD3D", "CD3E", "CD8A", "CD8B", "PDCD1"],
            "Exhaustion": ["LAG3", "TIGIT", "HAVCR2", "CTLA4"]
        }
    }
}
skill = GsvaScoresSkill()
result = skill.execute(state)
```

**技术栈：** scanpy + gseapy (ssGSEA)

## 目录结构

```
_synthesized/
├── SKILL.md            ← 本文件
├── custom_gsva_scores.py  ← 唯一真实实现
├── templates/           ← 模板目录（未实现）
│   └── (future skills)
└── contracts/           ← 合约目录（未实现）
    └── (future contracts)
```

## 开发指南

新增合成技能时：
1. 在 `_synthesized/` 下创建 `{skill_name}.py`
2. 继承 `AbstractSkill`，使用 `@register`
3. 在 `SKILL.md` 中添加条目
4. 测试与真实 scanpy AnnData 联合验证