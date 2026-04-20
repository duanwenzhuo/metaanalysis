# bioskills / core

Base infrastructure for the bioskills framework. Provides contract-based state dict system.

## Components

### base.py
Core classes and interfaces for bioskills operations.

### contract.py
State Dict contract definitions. Defines the interface between pipeline stages.
- **bioSkills**: `bio-state-dict-contract`

### pipeline.py
Pipeline orchestration engine. Connects bioskills modules into end-to-end workflows.
- **bioSkills**: `bio-single-cell-pipeline`

## Architecture

The bioskills framework uses a **State Dict** pattern:
1. Each skill receives a state dict
2. Transforms or annotates it
3. Returns modified state dict

This enables modular, composable single-cell analysis pipelines.

## bioSkills Reference

- `bio-state-dict-contract`
- `bio-single-cell-pipeline`