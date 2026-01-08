> CMU x NVIDIA Hackathon, January 7-9, 2026

# Pangenome and Genomic Cluster Modeling
# Problem 
Pangenomes have the potential to more acurately catpure genomic diversity across human, however current state of the art graphs are only on 10's of individuals. Large biobanks (UK Biobank, AllOfUS) have long read sequencing data on thousands of people that could be used in the service constructing a more diverse pangenome graph, however methods do not currently exists for pangome graph construction that enables cross-talk between biobanks. 

# Intro and aim

For this project we are aiming to complete two aims
1) Perform federated pangenome graph construction (using HPRC data as a proof of principle)
2) Perform federated genomic background hashing for phenotype association of APOE locus (using 1k genomes data)
   
# Contributors
- 

# How to use this repo

# Methods

# 1) Federated pangenome graph construction
## Flowchart

<img width="942" height="490" alt="image" src="https://github.com/user-attachments/assets/1e736b17-d779-4cd9-aa68-78d06b2d1740" />

## a) Dowload Data from HPRC

Dowload data from HPRC for generating graphs from:
```
download_path=/path/to/download/destination/
python ./HPRC_download_prep/download_hprc.py \
./HPRC_download_prep/assemblies.tsv \
$download_path
```
### Extract Chromsome 19 and 22 from HPRC samples

As pangenome graph construction is very comptuationally intensive we will be running the process on chromosome19 and chromsome22 as a more tractable dataset.

### Install Entrez Direct using conda
The contig names in the fasta files for HPRC are NCBI identifiers and need to be queried using edirect tools to convert these to convetional chromosme numbers (e.g. CM102454.1 -> chr22). This can be installed with conda using:
```
conda install bioconda::entrez-direct
```

From this then the following python script can be used to extract chr19 and chr22 from the assembly FASTA files for HPRC:
```
# This should be the path to the assemblies downloaded above
download_path='/path/to/download/destination/'
python ./HPRC_download_prep/make_hrcp_chr22_fasta.py \
/space/ceph/1/ABCD/users/rloughna/pangenome_construction/hprc_chr22_pansn_full.fa.gz \
--output-chr19 /space/ceph/1/ABCD/users/rloughna/pangenome_construction/hprc_chr19_pansn_full.fa.gz \
--n-individuals 20 \
--bgzip
```

## b) Create Partitions of Data to Simulate Biobank Cohorts
**WIP**

## c) Graph Construction 
### Prerequisites

- Docker CE with Compose plugin (v5.0+)
- Minimum 8GB RAM recommended

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

### Quick Start
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

### d) Aggregate Graphs from Individual Cohorts
**WIP**
`vg combine` from vg toolkit looks promising.

# 2) Genomic background hashing for phenotype association of APOE locus
<img width="980" height="490" alt="image" src="https://github.com/user-attachments/assets/9986a83d-8513-4914-8dbc-5a3087eb969c" />

We extracted this locus around the APOE gene which will be the region we use for localized pangenome graph mapping. Pangome graph mapping may then provide a genomic background to contexualize the high risk APOE alleles. This genomic background may capture trans expression effects and we will aim to code this using gemomic hashing to represent different anonymized haploblocks that could be used in a federated matter across studies. 

### Download GWAS Data
downloaded from [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST013197)
```
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST013001-GCST014000/GCST013197/GCST013197.tsv.gz
gunzip GCST013197.tsv.gz
```
<img width="927" height="702" alt="image" src="https://github.com/user-attachments/assets/6198d356-a1c6-4c1a-a55f-3596352a1f72" />
# Results

# Conclusions

# References
