#!/usr/bin/env python3
"""Extract APOE locus coordinates from Alzheimer's GWAS data"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import argparse

DEFAULT_LOCUS_PADDING = 50000
DEFAULT_P_VALUE_THRESHOLD = 1e-100

def main():
    parser = argparse.ArgumentParser(description='Extract APOE locus from GWAS data')
    parser.add_argument('--gwas', default='../gwas/data/GCST013197.tsv',
                       help='Path to GWAS data file')
    parser.add_argument('--output', default='../gwas/data',
                       help='Output directory for BED and plots')
    parser.add_argument('--padding', type=int, default=DEFAULT_LOCUS_PADDING,
                       help=f'Base pairs padding around SNPs (default: {DEFAULT_LOCUS_PADDING:,})')
    parser.add_argument('--pval', type=float, default=DEFAULT_P_VALUE_THRESHOLD,
                       help=f'P-value threshold (default: {DEFAULT_P_VALUE_THRESHOLD})')
    args = parser.parse_args()
    
    GWAS_FILE = args.gwas
    OUTPUT_DIR = args.output
    LOCUS_PADDING = args.padding
    P_VALUE_THRESHOLD = args.pval
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("Loading GWAS data...")
    try:
        gwas_df = pd.read_csv(GWAS_FILE, delimiter='\t')
    except FileNotFoundError:
        print(f"Error: GWAS file not found: {GWAS_FILE}")
        print("Please download the file first or check the path.")
        exit(1)
    except Exception as e:
        print(f"Error reading GWAS file: {e}")
        exit(1)

    print(f"Total SNPs: {len(gwas_df)}")

    chr19_df = gwas_df[gwas_df['chromosome'] == 19].copy()
    sig_snps = chr19_df[chr19_df['p_value'] < P_VALUE_THRESHOLD].sort_values(by='p_value').reset_index(drop=True)
    print(f"Significant SNPs (chr19, p < {P_VALUE_THRESHOLD}): {len(sig_snps)}")

    if len(sig_snps) > 0:
        bp_min = sig_snps['base_pair_location'].min() - LOCUS_PADDING
        bp_max = sig_snps['base_pair_location'].max() + LOCUS_PADDING
        
        print(f"\nAPOE Locus: chr19:{bp_min:,}-{bp_max:,} ({(bp_max - bp_min) / 1e6:.2f} Mb)")
        
        # Save BED file for odgi extract (zero-based, half-open intervals)
        bed_file = os.path.join(OUTPUT_DIR, 'apoe_locus.bed')
        with open(bed_file, 'w') as f:
            f.write(f"chr19\t{int(bp_min)}\t{int(bp_max)}\tAPOE_locus\n")
        
        # Save interval data
        interval_df = pd.DataFrame({
            'chromosome': ['19'],
            'start': [int(bp_min)],
            'end': [int(bp_max)],
            'width': [int(bp_max - bp_min)],
            'padding': [LOCUS_PADDING],
            'num_sig_snps': [len(sig_snps)]
        })
        interval_df.to_csv(os.path.join(OUTPUT_DIR, 'apoe_locus_interval.tsv'), sep='\t', index=False)
        
        # Plot locus
        chr19_df['p_value_log'] = -np.log10(np.clip(chr19_df['p_value'], np.finfo(float).tiny, None))
        locus_region = chr19_df[(chr19_df['base_pair_location'] >= bp_min) & 
                                (chr19_df['base_pair_location'] <= bp_max)]
        
        fig, ax = plt.subplots(figsize=(12, 6), dpi=150)
        sns.scatterplot(data=locus_region, x='base_pair_location', y='p_value_log', ax=ax, s=15, alpha=0.6)
        ax.axhline(-np.log10(5e-8), color="red", linestyle="--", label="Genome-wide significance")
        ax.set_xlabel('Position on Chromosome 19 (bp)')
        ax.set_ylabel('-log10(p-value)')
        ax.set_title("Alzheimer's Disease GWAS - APOE Locus")
        ax.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, 'apoe_locus_plot.png'), bbox_inches='tight')
        
        print(f"Saved: {bed_file}")
        print("\nExtraction complete!")

if __name__ == "__main__":
    main()