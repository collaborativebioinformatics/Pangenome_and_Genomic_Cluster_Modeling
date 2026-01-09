# Pangenome Graph Visualization Tool

Compares two pangenome graphs (GFA format) and generates 10 publication-ready visualizations.

## Quick Start
```bash
pip install matplotlib numpy --break-system-packages

python3 analyze_gfa.py <graph1.gfa> <graph2.gfa> [output_dir]
```

## Labels
- **Graph Iteration 1** (Gold) - Local/chunk graph from federated construction
- **Fully Converged Graph (HPRC)** (Green) - Reference HPRC pangenome

## Output Files
1. `01_node_length_distribution.png` - Side-by-side node length histograms
2. `02_node_length_overlay.png` - Overlayed distributions
3. `03_degree_distribution.png` - Node connectivity
4. `04_structure_comparison.png` - Nodes/Edges/Paths comparison
5. `05_size_metrics.png` - N50/N90/Total BP
6. `06_radar_comparison.png` - Normalized metrics radar
7. `07_node_types_pie.png` - Linear vs branching nodes
8. `08_cumulative_length.png` - N50/N90 curves
9. `09_dashboard.png` - Complete summary dashboard
10. `10_ratio_comparison.png` - Metric ratios

## Key Metrics
- **N50**: Length where 50% of sequence is in nodes >= this size (higher = better contiguity)
- **Branch Ratio**: Proportion of branching nodes (higher = more complex)
- **Edge/Node Ratio**: Graph density measure
