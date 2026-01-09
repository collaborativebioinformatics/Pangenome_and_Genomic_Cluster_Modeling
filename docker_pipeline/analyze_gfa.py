#!/usr/bin/env python3
import sys
import gzip
import json
import os
from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

class GFAParser:
    def __init__(self, filepath):
        self.filepath = filepath
        self.nodes = {}
        self.edges = []
        self.paths = {}
        self.node_lengths = []
        self.samples = set()
        print(f"Parsing {filepath}...")
        self.parse()
        print(f"  Done: {len(self.nodes):,} nodes, {len(self.edges):,} edges, {len(self.paths):,} paths")
    
    def parse(self):
        opener = gzip.open if str(self.filepath).endswith('.gz') else open
        mode = 'rt' if str(self.filepath).endswith('.gz') else 'r'
        with opener(self.filepath, mode) as f:
            for line in f:
                line = line.strip()
                if not line: continue
                parts = line.split('\t')
                if parts[0] == 'S':
                    self.nodes[parts[1]] = len(parts[2]) if len(parts) > 2 else 0
                    self.node_lengths.append(len(parts[2]) if len(parts) > 2 else 0)
                elif parts[0] == 'L':
                    self.edges.append({'from': parts[1], 'to': parts[3]})
                elif parts[0] == 'P':
                    self.paths[parts[1]] = parts[2].split(',')
                    if '#' in parts[1]: self.samples.add(parts[1].split('#')[0])

class GraphAnalyzer:
    def __init__(self, gfa, name):
        self.gfa = gfa
        self.name = name
        self.stats = {}
        self.degree_dist = []
        self.analyze()
    
    def analyze(self):
        s = self.stats
        s['name'] = self.name
        s['num_nodes'] = len(self.gfa.nodes)
        s['num_edges'] = len(self.gfa.edges)
        s['num_paths'] = len(self.gfa.paths)
        s['num_samples'] = len(self.gfa.samples)
        
        lengths = self.gfa.node_lengths
        if lengths:
            s['total_bp'] = sum(lengths)
            s['mean_node_len'] = np.mean(lengths)
            s['median_node_len'] = np.median(lengths)
            s['min_node_len'] = min(lengths)
            s['max_node_len'] = max(lengths)
            s['std_node_len'] = np.std(lengths)
            s['n50'] = self.calc_nx(lengths, 50)
            s['n90'] = self.calc_nx(lengths, 90)
        
        in_deg, out_deg = defaultdict(int), defaultdict(int)
        for e in self.gfa.edges:
            out_deg[e['from']] += 1
            in_deg[e['to']] += 1
        
        self.degree_dist = [in_deg[n] + out_deg[n] for n in self.gfa.nodes]
        if self.degree_dist:
            s['mean_degree'] = np.mean(self.degree_dist)
            s['max_degree'] = max(self.degree_dist)
            s['isolated_nodes'] = sum(1 for d in self.degree_dist if d == 0)
            s['linear_nodes'] = sum(1 for d in self.degree_dist if d == 2)
            s['branch_nodes'] = sum(1 for d in self.degree_dist if d > 2)
            s['branch_ratio'] = s['branch_nodes'] / s['num_nodes'] if s['num_nodes'] else 0
        
        s['edge_node_ratio'] = s['num_edges'] / s['num_nodes'] if s['num_nodes'] else 0
    
    def calc_nx(self, lengths, x):
        sorted_l = sorted(lengths, reverse=True)
        total = sum(sorted_l)
        target = total * x / 100
        cumsum = 0
        for l in sorted_l:
            cumsum += l
            if cumsum >= target: return l
        return 0

def get_graph_name(filepath):
    name = Path(filepath).name
    for ext in ['.gz', '.gfa', '.fa']:
        name = name.replace(ext, '')
    for suffix in ['.smooth.final', '.seqwish', '.14b29fc.11fba48', '.14b29fc']:
        name = name.replace(suffix, '')
    return name

def create_visualizations(gfa1, gfa2, a1, a2, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    c1, c2 = '#D4A017', '#2E8B57'
    
    print("Creating visualizations...")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    axes[0].hist(gfa1.node_lengths, bins=100, color=c1, alpha=0.8, edgecolor='black', linewidth=0.3)
    axes[0].set_title(f'{a1.name}\nNode Length Distribution', fontsize=14, fontweight='bold')
    axes[0].set_xlabel('Node Length (bp)'); axes[0].set_ylabel('Frequency'); axes[0].set_yscale('log'); axes[0].grid(True, alpha=0.3)
    axes[1].hist(gfa2.node_lengths, bins=100, color=c2, alpha=0.8, edgecolor='black', linewidth=0.3)
    axes[1].set_title(f'{a2.name}\nNode Length Distribution', fontsize=14, fontweight='bold')
    axes[1].set_xlabel('Node Length (bp)'); axes[1].set_ylabel('Frequency'); axes[1].set_yscale('log'); axes[1].grid(True, alpha=0.3)
    plt.tight_layout(); plt.savefig(f'{output_dir}/01_node_length_distribution.png', dpi=200, bbox_inches='tight'); plt.close()
    print("  01_node_length_distribution.png")
    
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.hist(gfa1.node_lengths, bins=100, color=c1, alpha=0.6, label=a1.name, edgecolor='black', linewidth=0.2)
    ax.hist(gfa2.node_lengths, bins=100, color=c2, alpha=0.6, label=a2.name, edgecolor='black', linewidth=0.2)
    ax.set_title('Node Length Distribution Comparison', fontsize=16, fontweight='bold')
    ax.set_xlabel('Node Length (bp)', fontsize=12); ax.set_ylabel('Frequency', fontsize=12); ax.set_yscale('log')
    ax.legend(fontsize=12); ax.grid(True, alpha=0.3)
    plt.tight_layout(); plt.savefig(f'{output_dir}/02_node_length_overlay.png', dpi=200, bbox_inches='tight'); plt.close()
    print("  02_node_length_overlay.png")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    deg1, deg2 = Counter(a1.degree_dist), Counter(a2.degree_dist)
    axes[0].bar(list(deg1.keys())[:20], list(deg1.values())[:20], color=c1, alpha=0.8, edgecolor='black')
    axes[0].set_title(f'{a1.name}\nNode Degree Distribution', fontsize=14, fontweight='bold')
    axes[0].set_xlabel('Degree'); axes[0].set_ylabel('Count'); axes[0].set_yscale('log'); axes[0].grid(True, alpha=0.3)
    axes[1].bar(list(deg2.keys())[:20], list(deg2.values())[:20], color=c2, alpha=0.8, edgecolor='black')
    axes[1].set_title(f'{a2.name}\nNode Degree Distribution', fontsize=14, fontweight='bold')
    axes[1].set_xlabel('Degree'); axes[1].set_ylabel('Count'); axes[1].set_yscale('log'); axes[1].grid(True, alpha=0.3)
    plt.tight_layout(); plt.savefig(f'{output_dir}/03_degree_distribution.png', dpi=200, bbox_inches='tight'); plt.close()
    print("  03_degree_distribution.png")
    
    fig, ax = plt.subplots(figsize=(14, 8))
    metrics = ['Nodes', 'Edges', 'Paths', 'Samples', 'Branch\nNodes']
    vals1 = [a1.stats['num_nodes'], a1.stats['num_edges'], a1.stats['num_paths'], a1.stats['num_samples'], a1.stats.get('branch_nodes', 0)]
    vals2 = [a2.stats['num_nodes'], a2.stats['num_edges'], a2.stats['num_paths'], a2.stats['num_samples'], a2.stats.get('branch_nodes', 0)]
    x = np.arange(len(metrics)); width = 0.35
    bars1 = ax.bar(x - width/2, vals1, width, label=a1.name, color=c1, alpha=0.8, edgecolor='black')
    bars2 = ax.bar(x + width/2, vals2, width, label=a2.name, color=c2, alpha=0.8, edgecolor='black')
    ax.set_ylabel('Count (log scale)', fontsize=12); ax.set_title('Graph Structure Comparison', fontsize=16, fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels(metrics, fontsize=11); ax.legend(fontsize=12); ax.set_yscale('log'); ax.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            h = bar.get_height()
            if h > 0: ax.annotate(f'{int(h):,}', xy=(bar.get_x() + bar.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=8, rotation=45)
    plt.tight_layout(); plt.savefig(f'{output_dir}/04_structure_comparison.png', dpi=200, bbox_inches='tight'); plt.close()
    print("  04_structure_comparison.png")
    
    fig, ax = plt.subplots(figsize=(14, 8))
    metrics = ['Total BP', 'Mean\nNode Len', 'Median\nNode Len', 'N50', 'N90', 'Max\nNode Len']
    vals1 = [a1.stats['total_bp'], a1.stats['mean_node_len'], a1.stats['median_node_len'], a1.stats['n50'], a1.stats['n90'], a1.stats['max_node_len']]
    vals2 = [a2.stats['total_bp'], a2.stats['mean_node_len'], a2.stats['median_node_len'], a2.stats['n50'], a2.stats['n90'], a2.stats['max_node_len']]
    x = np.arange(len(metrics))
    ax.bar(x - width/2, vals1, width, label=a1.name, color=c1, alpha=0.8, edgecolor='black')
    ax.bar(x + width/2, vals2, width, label=a2.name, color=c2, alpha=0.8, edgecolor='black')
    ax.set_ylabel('Base Pairs (log scale)', fontsize=12); ax.set_title('Sequence Size Metrics', fontsize=16, fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels(metrics, fontsize=10); ax.legend(fontsize=12); ax.set_yscale('log'); ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout(); plt.savefig(f'{output_dir}/05_size_metrics.png', dpi=200, bbox_inches='tight'); plt.close()
    print("  05_size_metrics.png")
    
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
    categories = ['Nodes\n(norm)', 'Edges\n(norm)', 'Paths\n(norm)', 'Branch\nRatio*10', 'Edge/Node\nRatio', 'Mean\nDegree/10']
    max_n = max(a1.stats['num_nodes'], a2.stats['num_nodes'])
    max_e = max(a1.stats['num_edges'], a2.stats['num_edges'])
    max_p = max(a1.stats['num_paths'], a2.stats['num_paths']) or 1
    vals1 = [a1.stats['num_nodes']/max_n, a1.stats['num_edges']/max_e, a1.stats['num_paths']/max_p, min(a1.stats['branch_ratio']*10, 1), min(a1.stats['edge_node_ratio'], 1), min(a1.stats.get('mean_degree',0)/10, 1)]
    vals2 = [a2.stats['num_nodes']/max_n, a2.stats['num_edges']/max_e, a2.stats['num_paths']/max_p, min(a2.stats['branch_ratio']*10, 1), min(a2.stats['edge_node_ratio'], 1), min(a2.stats.get('mean_degree',0)/10, 1)]
    angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False).tolist()
    vals1 += vals1[:1]; vals2 += vals2[:1]; angles += angles[:1]
    ax.plot(angles, vals1, 'o-', linewidth=2, label=a1.name, color=c1); ax.fill(angles, vals1, alpha=0.25, color=c1)
    ax.plot(angles, vals2, 'o-', linewidth=2, label=a2.name, color=c2); ax.fill(angles, vals2, alpha=0.25, color=c2)
    ax.set_xticks(angles[:-1]); ax.set_xticklabels(categories, fontsize=10)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=11)
    ax.set_title('Normalized Graph Metrics Comparison', fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout(); plt.savefig(f'{output_dir}/06_radar_comparison.png', dpi=200, bbox_inches='tight'); plt.close()
    print("  06_radar_comparison.png")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for ax, a, c, colors in [(axes[0], a1, c1, ['#FFD700', '#DAA520', '#B8860B']), (axes[1], a2, c2, ['#90EE90', '#3CB371', '#2E8B57'])]:
        sizes = [a.stats.get('isolated_nodes', 0), a.stats.get('linear_nodes', 0), a.stats.get('branch_nodes', 0)]
        labels = [f'Isolated\n(deg=0)\n{sizes[0]:,}', f'Linear\n(deg=2)\n{sizes[1]:,}', f'Branching\n(deg>2)\n{sizes[2]:,}']
        explode = (0.05, 0, 0.05)
        ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=colors, explode=explode, textprops={'fontsize': 10})
        ax.set_title(f'{a.name}\nNode Types by Connectivity', fontsize=14, fontweight='bold')
    plt.tight_layout(); plt.savefig(f'{output_dir}/07_node_types_pie.png', dpi=200, bbox_inches='tight'); plt.close()
    print("  07_node_types_pie.png")
    
    fig, ax = plt.subplots(figsize=(14, 7))
    sorted1, sorted2 = np.sort(gfa1.node_lengths)[::-1], np.sort(gfa2.node_lengths)[::-1]
    cum1, cum2 = np.cumsum(sorted1) / sum(sorted1) * 100, np.cumsum(sorted2) / sum(sorted2) * 100
    ax.plot(range(len(cum1)), cum1, color=c1, linewidth=2, label=a1.name)
    ax.plot(range(len(cum2)), cum2, color=c2, linewidth=2, label=a2.name)
    ax.axhline(y=50, color='red', linestyle='--', alpha=0.7, linewidth=2, label='N50 threshold')
    ax.axhline(y=90, color='orange', linestyle='--', alpha=0.7, linewidth=2, label='N90 threshold')
    ax.set_xlabel('Node Rank (sorted by length, descending)', fontsize=12)
    ax.set_ylabel('Cumulative Sequence (%)', fontsize=12)
    ax.set_title('Cumulative Node Length Distribution (N50/N90 Analysis)', fontsize=16, fontweight='bold')
    ax.legend(fontsize=11); ax.grid(True, alpha=0.3); ax.set_xscale('log')
    plt.tight_layout(); plt.savefig(f'{output_dir}/08_cumulative_length.png', dpi=200, bbox_inches='tight'); plt.close()
    print("  08_cumulative_length.png")
    
    fig = plt.figure(figsize=(24, 18))
    fig.suptitle('PANGENOME GRAPH COMPARISON DASHBOARD', fontsize=24, fontweight='bold', y=0.98)
    gs = fig.add_gridspec(4, 4, hspace=0.35, wspace=0.3)
    
    ax1 = fig.add_subplot(gs[0, 0:2]); ax1.axis('off')
    txt1 = f"{a1.name}\n{'='*40}\nNodes:        {a1.stats['num_nodes']:>15,}\nEdges:        {a1.stats['num_edges']:>15,}\nPaths:        {a1.stats['num_paths']:>15,}\nSamples:      {a1.stats['num_samples']:>15,}\nTotal BP:     {a1.stats['total_bp']:>15,}\nN50:          {a1.stats['n50']:>15,}\nN90:          {a1.stats['n90']:>15,}\nBranch Nodes: {a1.stats.get('branch_nodes',0):>15,}\nBranch Ratio: {a1.stats['branch_ratio']:>15.4f}\nMean Degree:  {a1.stats.get('mean_degree',0):>15.2f}"
    ax1.text(0.05, 0.5, txt1, transform=ax1.transAxes, fontsize=11, verticalalignment='center', fontfamily='monospace', bbox=dict(boxstyle='round,pad=0.5', facecolor=c1, alpha=0.3))
    
    ax2 = fig.add_subplot(gs[0, 2:4]); ax2.axis('off')
    txt2 = f"{a2.name}\n{'='*40}\nNodes:        {a2.stats['num_nodes']:>15,}\nEdges:        {a2.stats['num_edges']:>15,}\nPaths:        {a2.stats['num_paths']:>15,}\nSamples:      {a2.stats['num_samples']:>15,}\nTotal BP:     {a2.stats['total_bp']:>15,}\nN50:          {a2.stats['n50']:>15,}\nN90:          {a2.stats['n90']:>15,}\nBranch Nodes: {a2.stats.get('branch_nodes',0):>15,}\nBranch Ratio: {a2.stats['branch_ratio']:>15.4f}\nMean Degree:  {a2.stats.get('mean_degree',0):>15.2f}"
    ax2.text(0.05, 0.5, txt2, transform=ax2.transAxes, fontsize=11, verticalalignment='center', fontfamily='monospace', bbox=dict(boxstyle='round,pad=0.5', facecolor=c2, alpha=0.3))
    
    ax3 = fig.add_subplot(gs[1, 0:2])
    m = ['Nodes', 'Edges', 'Paths', 'Branch']; v1 = [a1.stats['num_nodes'], a1.stats['num_edges'], a1.stats['num_paths'], a1.stats.get('branch_nodes',0)]
    v2 = [a2.stats['num_nodes'], a2.stats['num_edges'], a2.stats['num_paths'], a2.stats.get('branch_nodes',0)]; x = np.arange(len(m))
    ax3.bar(x-0.2, v1, 0.4, label=a1.name, color=c1, alpha=0.8); ax3.bar(x+0.2, v2, 0.4, label=a2.name, color=c2, alpha=0.8)
    ax3.set_xticks(x); ax3.set_xticklabels(m); ax3.set_yscale('log'); ax3.legend(); ax3.grid(True, alpha=0.3, axis='y'); ax3.set_title('Structure Comparison', fontweight='bold')
    
    ax4 = fig.add_subplot(gs[1, 2:4])
    ax4.hist(gfa1.node_lengths, bins=50, alpha=0.6, color=c1, label=a1.name)
    ax4.hist(gfa2.node_lengths, bins=50, alpha=0.6, color=c2, label=a2.name)
    ax4.set_yscale('log'); ax4.legend(); ax4.grid(True, alpha=0.3); ax4.set_title('Node Length Distribution', fontweight='bold'); ax4.set_xlabel('Length (bp)')
    
    ax5 = fig.add_subplot(gs[2, 0:2])
    deg1_k = sorted(deg1.keys())[:15]; deg2_k = sorted(deg2.keys())[:15]
    ax5.bar([k-0.2 for k in deg1_k], [deg1[k] for k in deg1_k], 0.4, color=c1, alpha=0.8, label=a1.name)
    ax5.bar([k+0.2 for k in deg2_k], [deg2[k] for k in deg2_k], 0.4, color=c2, alpha=0.8, label=a2.name)
    ax5.set_yscale('log'); ax5.grid(True, alpha=0.3); ax5.set_title('Degree Distribution', fontweight='bold'); ax5.set_xlabel('Degree'); ax5.legend()
    
    ax6 = fig.add_subplot(gs[2, 2:4])
    m = ['Total BP', 'N50', 'N90', 'Mean Len']
    v1 = [a1.stats['total_bp'], a1.stats['n50'], a1.stats['n90'], a1.stats['mean_node_len']]
    v2 = [a2.stats['total_bp'], a2.stats['n50'], a2.stats['n90'], a2.stats['mean_node_len']]
    x = np.arange(len(m))
    ax6.bar(x-0.2, v1, 0.4, label=a1.name, color=c1, alpha=0.8); ax6.bar(x+0.2, v2, 0.4, label=a2.name, color=c2, alpha=0.8)
    ax6.set_xticks(x); ax6.set_xticklabels(m); ax6.set_yscale('log'); ax6.legend(); ax6.grid(True, alpha=0.3, axis='y'); ax6.set_title('Size Metrics', fontweight='bold')
    
    ax7 = fig.add_subplot(gs[3, 0:2])
    ax7.plot(range(len(cum1)), cum1, color=c1, linewidth=2, label=a1.name)
    ax7.plot(range(len(cum2)), cum2, color=c2, linewidth=2, label=a2.name)
    ax7.axhline(y=50, color='red', linestyle='--', alpha=0.5); ax7.axhline(y=90, color='orange', linestyle='--', alpha=0.5)
    ax7.set_xscale('log'); ax7.legend(); ax7.grid(True, alpha=0.3); ax7.set_title('Cumulative Length (N50/N90)', fontweight='bold'); ax7.set_xlabel('Node Rank')
    
    ax8 = fig.add_subplot(gs[3, 2:4]); ax8.axis('off')
    n1, n2 = a1.stats['num_nodes'], a2.stats['num_nodes']
    e1, e2 = a1.stats['num_edges'], a2.stats['num_edges']
    larger_n = a2.name if n2 > n1 else a1.name
    larger_e = a2.name if e2 > e1 else a1.name
    ratio_n = max(n1,n2)/min(n1,n2) if min(n1,n2) > 0 else 0
    ratio_e = max(e1,e2)/min(e1,e2) if min(e1,e2) > 0 else 0
    more_complex = a2.name if a2.stats['branch_ratio'] > a1.stats['branch_ratio'] else a1.name
    insights = f"KEY INSIGHTS\n{'='*50}\n\nScale Difference:\n  - {larger_n} has {ratio_n:.1f}x MORE nodes\n  - {larger_e} has {ratio_e:.1f}x MORE edges\n\nComplexity Analysis:\n  - {a1.name}: Branch ratio = {a1.stats['branch_ratio']:.4f}\n  - {a2.name}: Branch ratio = {a2.stats['branch_ratio']:.4f}\n  - {more_complex} is MORE COMPLEX\n\nContiguity (N50):\n  - {a1.name}: N50 = {a1.stats['n50']:,} bp\n  - {a2.name}: N50 = {a2.stats['n50']:,} bp\n\nSample Coverage:\n  - {a1.name}: {a1.stats['num_samples']} samples, {a1.stats['num_paths']} paths\n  - {a2.name}: {a2.stats['num_samples']} samples, {a2.stats['num_paths']} paths"
    ax8.text(0.02, 0.98, insights, transform=ax8.transAxes, fontsize=10, verticalalignment='top', fontfamily='monospace', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.savefig(f'{output_dir}/09_dashboard.png', dpi=200, bbox_inches='tight'); plt.close()
    print("  09_dashboard.png")
    
    fig, ax = plt.subplots(figsize=(10, 12))
    metrics_h = ['num_nodes', 'num_edges', 'num_paths', 'num_samples', 'total_bp', 'mean_node_len', 'median_node_len', 'n50', 'n90', 'branch_nodes', 'branch_ratio', 'edge_node_ratio', 'mean_degree']
    labels_h = ['Nodes', 'Edges', 'Paths', 'Samples', 'Total BP', 'Mean Node Len', 'Median Node Len', 'N50', 'N90', 'Branch Nodes', 'Branch Ratio', 'Edge/Node', 'Mean Degree']
    ratios = []
    for m in metrics_h:
        v1, v2 = a1.stats.get(m, 0), a2.stats.get(m, 0)
        ratios.append(v1/v2 if v2 != 0 else 1)
    
    colors = ['#ff4444' if r < 0.5 else '#ffaa44' if r < 0.9 else '#88cc88' if r < 1.1 else '#44aa44' if r < 2 else '#228822' for r in ratios]
    y_pos = np.arange(len(labels_h))
    bars = ax.barh(y_pos, ratios, color=colors, edgecolor='black', alpha=0.8)
    ax.axvline(x=1, color='black', linestyle='--', linewidth=2, label='Equal (1.0)')
    ax.set_yticks(y_pos); ax.set_yticklabels(labels_h, fontsize=11)
    ax.set_xlabel('Ratio (Graph1 / Graph2)', fontsize=12)
    ax.set_title(f'Metric Ratios: {a1.name} / {a2.name}', fontsize=14, fontweight='bold')
    for i, (bar, r) in enumerate(zip(bars, ratios)):
        ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2, f'{r:.3f}x', va='center', fontsize=10, fontweight='bold')
    ax.set_xlim(0, max(ratios) * 1.3)
    ax.grid(True, alpha=0.3, axis='x')
    plt.tight_layout(); plt.savefig(f'{output_dir}/10_ratio_comparison.png', dpi=200, bbox_inches='tight'); plt.close()
    print("  10_ratio_comparison.png")

def generate_report(a1, a2, output_dir):
    r = []
    r.append("=" * 80)
    r.append("              PANGENOME GRAPH COMPARISON REPORT")
    r.append(f"              Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    r.append("=" * 80)
    
    r.append(f"\n{'-'*80}")
    r.append(f"GRAPH 1: {a1.name}")
    r.append("-" * 80)
    for k, v in sorted(a1.stats.items()):
        if k == 'name': continue
        r.append(f"  {k:30s}: {v:>20,.4f}" if isinstance(v, float) else f"  {k:30s}: {v:>20,}")
    
    r.append(f"\n{'-'*80}")
    r.append(f"GRAPH 2: {a2.name}")
    r.append("-" * 80)
    for k, v in sorted(a2.stats.items()):
        if k == 'name': continue
        r.append(f"  {k:30s}: {v:>20,.4f}" if isinstance(v, float) else f"  {k:30s}: {v:>20,}")
    
    r.append(f"\n{'='*80}")
    r.append("COMPARISON TABLE")
    r.append("=" * 80)
    r.append(f"{'Metric':<25s} | {a1.name[:18]:>18s} | {a2.name[:18]:>18s} | {'Ratio':>10s}")
    r.append("-" * 80)
    
    for k in sorted(a1.stats.keys()):
        if k == 'name': continue
        v1, v2 = a1.stats.get(k, 0), a2.stats.get(k, 0)
        ratio = v1/v2 if v2 != 0 else 0
        if isinstance(v1, float):
            r.append(f"{k:<25s} | {v1:>18,.4f} | {v2:>18,.4f} | {ratio:>9.3f}x")
        else:
            r.append(f"{k:<25s} | {v1:>18,} | {v2:>18,} | {ratio:>9.3f}x")
    
    r.append(f"\n{'='*80}")
    r.append("KEY INSIGHTS")
    r.append("=" * 80)
    
    n1, n2 = a1.stats['num_nodes'], a2.stats['num_nodes']
    e1, e2 = a1.stats['num_edges'], a2.stats['num_edges']
    
    r.append(f"\nSCALE COMPARISON:")
    if n2 > n1:
        r.append(f"   - {a2.name} has {n2/n1:.1f}x MORE nodes ({n2:,} vs {n1:,})")
    else:
        r.append(f"   - {a1.name} has {n1/n2:.1f}x MORE nodes ({n1:,} vs {n2:,})")
    if e2 > e1:
        r.append(f"   - {a2.name} has {e2/e1:.1f}x MORE edges ({e2:,} vs {e1:,})")
    else:
        r.append(f"   - {a1.name} has {e1/e2:.1f}x MORE edges ({e1:,} vs {e2:,})")
    
    r.append(f"\nCOMPLEXITY ANALYSIS:")
    br1, br2 = a1.stats['branch_ratio'], a2.stats['branch_ratio']
    r.append(f"   - {a1.name} branch ratio: {br1:.4f} ({a1.stats.get('branch_nodes',0):,} branching nodes)")
    r.append(f"   - {a2.name} branch ratio: {br2:.4f} ({a2.stats.get('branch_nodes',0):,} branching nodes)")
    r.append(f"   - {a2.name if br2 > br1 else a1.name} shows HIGHER complexity")
    
    r.append(f"\nCONTIGUITY (N50/N90):")
    r.append(f"   - {a1.name}: N50={a1.stats['n50']:,} bp, N90={a1.stats['n90']:,} bp")
    r.append(f"   - {a2.name}: N50={a2.stats['n50']:,} bp, N90={a2.stats['n90']:,} bp")
    
    r.append(f"\nSAMPLE COVERAGE:")
    r.append(f"   - {a1.name}: {a1.stats['num_samples']} unique samples, {a1.stats['num_paths']} paths")
    r.append(f"   - {a2.name}: {a2.stats['num_samples']} unique samples, {a2.stats['num_paths']} paths")
    
    r.append(f"\nSEQUENCE STATISTICS:")
    r.append(f"   - {a1.name}: Total={a1.stats['total_bp']/1e6:.2f}Mb, Mean node={a1.stats['mean_node_len']:.0f}bp")
    r.append(f"   - {a2.name}: Total={a2.stats['total_bp']/1e6:.2f}Mb, Mean node={a2.stats['mean_node_len']:.0f}bp")
    
    r.append(f"\nCONNECTIVITY:")
    r.append(f"   - {a1.name}: Mean degree={a1.stats.get('mean_degree',0):.2f}, Edge/node ratio={a1.stats['edge_node_ratio']:.2f}")
    r.append(f"   - {a2.name}: Mean degree={a2.stats.get('mean_degree',0):.2f}, Edge/node ratio={a2.stats['edge_node_ratio']:.2f}")
    
    r.append("\n" + "=" * 80)
    
    txt = "\n".join(r)
    print(txt)
    
    with open(f'{output_dir}/comparison_report.txt', 'w') as f:
        f.write(txt)
    
    json_data = {
        'graph1': {'name': a1.name, 'stats': {k: float(v) if isinstance(v, (int, float, np.integer, np.floating)) else v for k, v in a1.stats.items()}},
        'graph2': {'name': a2.name, 'stats': {k: float(v) if isinstance(v, (int, float, np.integer, np.floating)) else v for k, v in a2.stats.items()}},
        'generated': datetime.now().isoformat()
    }
    with open(f'{output_dir}/comparison_data.json', 'w') as f:
        json.dump(json_data, f, indent=2)
    
    print(f"\nReport saved: {output_dir}/comparison_report.txt")
    print(f"JSON saved: {output_dir}/comparison_data.json")

def main():
    if len(sys.argv) < 3:
        print("Usage: python analyze_gfa.py <gfa1> <gfa2> [output_dir]")
        print("Example: python analyze_gfa.py chr19_chunk1.gfa chr19.hprc-v1.0-pggb.gfa.gz ./comparison")
        sys.exit(1)
    
    gfa1_path, gfa2_path = sys.argv[1], sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "./comparison"
    
    name1 = get_graph_name(gfa1_path)
    name2 = get_graph_name(gfa2_path)
    
    print("\n" + "="*60)
    print("     PANGENOME GRAPH COMPARISON ANALYSIS")
    print("="*60 + "\n")
    
    gfa1 = GFAParser(gfa1_path)
    gfa2 = GFAParser(gfa2_path)
    
    print("\nAnalyzing graphs...")
    a1 = GraphAnalyzer(gfa1, name1)
    a2 = GraphAnalyzer(gfa2, name2)
    
    print(f"\nGenerating visualizations to {output_dir}/...")
    create_visualizations(gfa1, gfa2, a1, a2, output_dir)
    
    print("\nGenerating report...")
    generate_report(a1, a2, output_dir)
    
    print(f"\n{'='*60}")
    print(f"ANALYSIS COMPLETE!")
    print(f"Output directory: {output_dir}/")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
