#!/usr/bin/env python3
"""
OmniGenome: Genomic Dot Plot Generator
--------------------------------------
Optimized for APOE haploblock clusters (~50kb).
Uses k-mer hashing 

"""

import sys
import gzip
import argparse
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser(description="Generate a genomic dot plot.")
    parser.add_argument("fasta", help="Input FASTA or FASTA.GZ file")
    parser.add_argument("output", help="Output PNG filename")
    parser.add_argument("-k", "--kmer", type=int, default=15, 
                        help="K-mer size (15 is recommended for APOE)")
    return parser.parse_args()

def find_kmer_matches(seq1, seq2, k):
    """
    Builds a hash map of seq1 and scans seq2 for matches.
    
    """
    # 1. Index the first sequence
    lookup = defaultdict(list)
    for i in range(len(seq1) - k + 1):
        kmer = seq1[i : i + k]
        lookup[kmer].append(i)
    
    # 2. Search using the second sequence
    x_points, y_points = [], []
    for j in range(len(seq2) - k + 1):
        kmer = seq2[j : j + k]
        if kmer in lookup:
            for i in lookup[kmer]:
                x_points.append(j)
                y_points.append(i)
                
    return x_points, y_points

def main():
    args = get_args()

    # Load FASTA (Handles both .gz or .fa)
    opener = gzip.open if args.fasta.endswith(".gz") else open
    try:
        with opener(args.fasta, "rt") as f:
            records = list(SeqIO.parse(f, "fasta"))
    except Exception as e:
        sys.exit(f"Error reading file: {e}")

    if len(records) < 2:
        sys.exit("Error: Need at least 2 sequences to compare.")

    # Select the cluster 
    rec1, rec2 = records[0], records[1]
    s1_seq, s2_seq = str(rec1.seq).upper(), str(rec2.seq).upper()

    print(f"Comparing {rec1.id} vs {rec2.id}...")

    # Run hashing logic
    x, y = find_kmer_matches(s1_seq, s2_seq, args.kmer)

    # Visualization
    plt.figure(figsize=(10, 10))
    plt.scatter(x, y, s=0.5, c='black', marker='.', alpha=0.6)
    
    # Labels & Styling
    plt.title(f"OmniGenome Dot Plot\n{rec1.id} vs {rec2.id} (k={args.kmer})")
    plt.xlabel(f"{rec2.id} (bp)")
    plt.ylabel(f"{rec1.id} (bp)")
    
    # Invert Y-axis to match genome browsers (0 at top)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    
    plt.savefig(args.output, dpi=300)
    print(f"Plot saved to {args.output}")

if __name__ == "__main__":

    main()

