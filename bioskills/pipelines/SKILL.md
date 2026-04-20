# bioskills / pipelines

End-to-end single-cell analysis pipelines orchestrated by bioskills core.

## Skills

### cell_state_pipeline.py
Complete cell state identification pipeline: QC → Normalize → HVG → PCA → Neighbors → Leiden → Annotation → Visualization.
- **bioSkills**: `bio-cell-state-pipeline`

## Trigger Phrases

- "Run the complete pipeline"
- "Full analysis from raw data"
- "Cell state identification"
- "End-to-end single-cell analysis"

## Pipeline Stages

1. **QC**: Filter low-quality cells, remove doublets
2. **Normalize**: Library size normalization + log transform
3. **HVG**: Select top 2000 highly variable genes
4. **PCA**: Linear dimensionality reduction (50 PCs)
5. **Neighbors**: Build SNN graph (k=20)
6. **Cluster**: Leiden clustering (res=0.5)
7. **Annotation**: Marker gene / reference-based cell typing
8. **Visualize**: UMAP + marker expression heatmap

## bioSkills Reference

- `bio-cell-state-pipeline`