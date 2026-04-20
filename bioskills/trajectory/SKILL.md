# bioskills / trajectory

Trajectory inference skills for single-cell data.

## Skills

### trajectory_inference.py
Trajectory and pseudotime inference using Monocle3 or scVelo.
- **Classes**: `TrajectoryInferenceSkill`
- **bioSkills**: `bio-single-cell-trajectory-inference`

### paga.py
PAGA (Partition-based Graph Abstraction) for trajectory topology analysis.
- **Classes**: `PAGASkill`, `PAGAForkSkill`
- **PAGA**: Fast graph-based trajectory structure inference
- **PAGAFork**: Branch point detection and fate decision analysis
- **Dependencies**: `scanpy`, `networkx`
- **bioSkills**: `bio-paga-trajectory`, `bio-branch-detection`

## Trigger Phrases

- "Find the developmental trajectory"
- "Infer pseudotime ordering"
- "Identify cell trajectories"
- "Which cells are most differentiated"
- "PAGA analysis"
- "Detect branch points"
- "Fate decision analysis"

## Method Comparison

| Tool      | Speed | Scalability | Best For |
|-----------|-------|-------------|----------|
| Monocle3  | Slow  | Good        | Trajectory + pseudotime |
| scVelo    | Slow  | Moderate    | RNA velocity (dynamics) |
| PAGA      | Fast  | Excellent   | Large-scale trajectories |
| PAGAFork  | Fast  | Excellent   | Branch detection |

## bioSkills Reference

- `bio-single-cell-trajectory-inference`
- `bio-paga-trajectory`
- `bio-branch-detection`
