#!/usr/bin/env python3
"""Build APOE locus pangenome graph from 1000 Genomes data"""

import subprocess
import os
import sys
from pathlib import Path

APOE_START = 45274138
APOE_END = 45479708
CHR = "19"
OUTPUT_DIR = "../data"
GRAPH_OUTPUT = "../results/graphs/APOE_test"

# 1000 Genomes samples to use (diverse ancestry)
# Selected to represent broad population diversity:
# - HG00096, HG00097: GBR (British, European ancestry)
# - NA19625, NA19648: ASW (African American, admixed African/European)
# - HG01879: ACB (African Caribbean, African ancestry)
# - NA12878: CEU (Utah European, commonly used reference sample)
SAMPLES = [
    "HG00096",
    "HG00097",
    "NA19625",
    "NA19648",
    "HG01879",
    "NA12878",
]

def run_cmd(cmd, desc):
    """Execute shell command and check for errors"""
    print(f"\n{'='*60}")
    print(desc)
    print(f"{'='*60}")
    print(f"$ {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: Command failed with exit code {result.returncode}")
        print(f"stderr: {result.stderr}")
        return False
    if result.stdout:
        print(result.stdout)
    return True

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(GRAPH_OUTPUT, exist_ok=True)
    
    print("APOE Pangenome Graph Construction")
    print("="*60)
    print(f"Region: chr{CHR}:{APOE_START}-{APOE_END}")
    print(f"Samples: {len(SAMPLES)}")
    
    # Download reference genome for APOE region
    ref_url = f"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr{CHR}.fa.gz"
    ref_file = f"{OUTPUT_DIR}/chr{CHR}_ref.fa.gz"
    
    if not Path(ref_file).exists():
        print(f"\nDownloading reference chr{CHR}...")
        run_cmd(f"wget -q -O {ref_file} {ref_url}", "Download reference")
    
    # Extract APOE region from reference
    apoe_ref = f"{OUTPUT_DIR}/apoe_ref.fa"
    if not Path(apoe_ref).exists():
        print(f"\nExtracting APOE region from reference...")
        run_cmd(f"gunzip -c {ref_file} | samtools faidx - chr{CHR}:{APOE_START}-{APOE_END} > {apoe_ref}", 
                "Extract APOE region")
    
    print("\n" + "="*60)
    print("NOTE: Full 1000 Genomes download requires samtools/bcftools")
    print("For demonstration, we'll use the reference sequence only")
    print("="*60)
    
    # Build graph from reference
    print("\nBuilding pangenome graph with PGGB...")
    graph_cmd = f"""docker compose -f ../docker/docker-compose.yml run pggb pggb \\
        sys.exit(1)ta/apoe_ref.fa \\
        -o /output/APOE_test \\
        -n 2 -t 4 -p 95 -s 5000"""
    
    if run_cmd(graph_cmd, "Build APOE graph"):
        print("\n✓ Graph construction complete!")
        print(f"Output: {GRAPH_OUTPUT}/")
    else:
        print("\n✗ Graph construction failed")
        return
    
    # Find the generated graph file
    og_files = list(Path(GRAPH_OUTPUT).glob("*.og"))
    if og_files:
        print(f"\nGenerated graph: {og_files[0].name}")
        print("\nNext: Run extract_graph.sh to extract and index")
    
if __name__ == "__main__":
    main()
