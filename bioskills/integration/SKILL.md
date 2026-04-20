# bioskills / integration

Batch integration and data fusion skills.

## Skills

### batch_correction.py
Remove batch effects using Harmony, ComBat-Seq, or Seurat CCA.
- **Classes**: `BatchCorrectionSkill`
- **bioSkills**: `bio-single-cell-batch-integration`

### scvi.py
Deep learning-based batch correction using scVI (single-cell Variational Inference).
- **Classes**: `ScVISkill`, `ScanVISkill`
- **scVI**: Unsupervised integration via variational autoencoder
- **scANVI**: Semi-supervised integration with label transfer
- **Dependencies**: `scvi-tools`
- **bioSkills**: `bio-scvi-integration`, `bio-scanvi-label-transfer`

## Trigger Phrases

- "Integrate my batches"
- "Remove batch effects"
- "scVI integration"
- "Label transfer between batches"
- "Correct for sample heterogeneity"

## Method Comparison

| Tool       | Speed | Scalability | Best For |
|------------|-------|-------------|----------|
| Harmony    | Fast  | Good        | Quick integration |
| scVI       | Slow  | Excellent   | Large datasets (>500k) |
| scANVI     | Slow  | Excellent   | + Label transfer |
| ComBat-Seq | Fast | Good        | Count data (RNA-seq) |
| Seurat CCA | Slow  | Good        | Conserved biology |

## bioSkills Reference

- `bio-single-cell-batch-integration`
- `bio-differential-expression-batch-correction`
- `bio-scvi-integration`
- `bio-scanvi-label-transfer`
