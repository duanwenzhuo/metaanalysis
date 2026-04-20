# bioskills / cluster

Cell clustering and community detection for single-cell data.

## Skills

### leiden.py
Leiden algorithm clustering via scanpy. Supports multi-resolution and directed Leiden.
- **bioSkills**: `bio-single-cell-leiden-clustering`

### louvain.py
Louvain community detection via scanpy. Well-suited for large datasets.
- **bioSkills**: `bio-single-cell-louvain-clustering`

## Trigger Phrases

- "Cluster the cells"
- "Community detection"
- "Identify cell populations"
- "Leiden / Louvain clustering"
- "Cell grouping"

## Method Comparison

| Algorithm | Speed | Resolution | Best For |
|-----------|-------|------------|----------|
| Leiden    | Medium | High flexibility | Most use cases |
| Louvain   | Fast  | Moderate | Large datasets |

## bioSkills Reference

- `bio-single-cell-leiden-clustering`
- `bio-single-cell-louvain-clustering`