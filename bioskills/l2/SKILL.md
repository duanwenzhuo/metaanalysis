# bioskills / l2

Level-2 meta-learning and self-healing pipeline for single-cell analysis.

## Skills

### self_heal.py
Adaptive pipeline that uses LLM reasoning to heal sub-optimal analysis results.
- **bioSkills**: `bio-self-heal`

## Trigger Phrases

- "Self-heal the results"
- "Fix failed clustering"
- "Adaptive pipeline"
- "Re-run with corrections"

## Description

The `self_heal` module implements a closed-loop pipeline:
1. Runs initial analysis
2. Evaluates quality metrics
3. Uses LLM reasoning to identify issues
4. Adjusts parameters and re-runs
5. Repeats until quality threshold met

## bioSkills Reference

- `bio-self-heal`