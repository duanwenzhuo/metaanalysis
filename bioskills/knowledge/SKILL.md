# bioskills / knowledge

Gene set databases, enrichment tooling, and regulatory network inference.

## Skills

### gene_set_db.py
Local gene set database supporting Cell Marker, biological process, and custom collections.
- **Classes**: `GeneSetDBSkill`
- **bioSkills**: `bio-gene-set-db`

### go_enricher.py
Gene Ontology overrepresentation analysis using biological process annotations.
- **Classes**: `GOEnricherSkill`
- **bioSkills**: `bio-go-enrichment`

### msigdb_client.py
MSigDB API client for fetching hallmark, C2, C5 gene sets.
- **Classes**: `MSigDBClientSkill`
- **bioSkills**: `bio-msigdb`

### scenic.py
SCENIC (Single-Cell Regulatory Network Inference and Clustering) for gene regulatory network analysis.
- **Classes**: `SCENICSkill`, `RegulonScoreSkill`
- **SCENIC**: Infers transcription factor regulons and gene regulatory networks
- **RegulonScore**: Scores cells by regulon activity (AUCell-like)
- **Dependencies**: `pyscenic` (optional), `scanpy`
- **bioSkills**: `bio-scenic-grn`, `bio-regulon-activity`

## Trigger Phrases

- "Query gene set database"
- "Find marker genes"
- "Cell marker database"
- "MSigDB gene sets"
- "GO enrichment"
- "Gene regulatory network"
- "SCENIC analysis"
- "Regulon activity"
- "Transcription factor targets"

## bioSkills Reference

- `bio-gene-set-db`
- `bio-go-enrichment`
- `bio-msigdb`
- `bio-scenic-grn`
- `bio-regulon-activity`
