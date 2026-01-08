# Docker Setup for Pangenome Graph Construction

This directory contains Docker configurations for running PGGB and vg tools.

## Prerequisites

- Docker CE with Compose plugin (v5.0+)
- Minimum 8GB RAM recommended

## Available Services

### PGGB (Pangenome Graph Builder)
Builds pangenome graphs from multi-sample FASTA files.
```bash
docker compose run pggb pggb \
  -i /data/<input.fa.gz> \
  -o /output/<run_name> \
  -n <num_haplotypes> \
  -t 8 -p 90 -s 10000
```

### VG (Variation Graph Toolkit)
Converts GFA graphs to Giraffe-compatible formats (GBZ, dist, min).
```bash
docker compose run vg autoindex --workflow giraffe \
  -g /data/<graph.gfa> \
  -p /data/<output_prefix>
```

## Quick Start
```bash
# 1. Start Docker daemon
systemctl start docker

# 2. Build containers
docker compose build

# 3. Run PGGB to build graph
docker compose run pggb pggb \
  -i /data/input.fa.gz \
  -o /output/my_graph \
  -n 12 -t 8 -p 90 -s 10000

# 4. Convert to Giraffe format
docker compose run vg autoindex --workflow giraffe \
  -g /data/my_graph/*.smooth.final.gfa \
  -p /data/my_graph/giraffe_index
```
## Input Data Format

Input FASTA files (.fa.gz) should contain multiple sequences from the same genomic region:
```
>Individual1#Chr1
ACTGACTGACTG...
>Individual2#Chr1
ACTGACTGACTG...
```

### Example: HLA Test Data
The HLA-DRB1 test dataset contains 12 sequences from chromosome 6:
- 1 primary reference (GRCh38)
- 5 alternate locus haplotypes (COX, DBB, MANN, QBL, SSTO)
- 4 MHC class II haplotypes (DR51, DR52, DR53, HLA-DRB1)
- 2 alternate assemblies (CHM1, HuRef)

Total: 163,416 bp across 12 sequences

## Output Files

| File | Description | Used By |
|------|-------------|---------|
| .gfa | Graph in GFA format | Visualization, analysis |
| .og | ODGI binary format | odgi tools |
| .gbz | Giraffe graph format | vg giraffe |
| .dist | Distance index | vg giraffe |
| .min | Minimizer index | vg giraffe |
| .zipcodes | Zipcodes | vg giraffe |

## Volume Mounts

- ../data -> /data (input files, read-only)
- ../results/graphs -> /output (output files)
