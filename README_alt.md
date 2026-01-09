> **CMU × NVIDIA Hackathon — January 7–9, 2026**

<img width="256" height="256" alt="omnigenome_logo" src="https://github.com/user-attachments/assets/5d2e5108-603f-4758-82c5-35e6172110d1" />

# OmniGenome  
### Federated Pangenomes & Genomic Background Modeling

![Python](https://img.shields.io/badge/Python-3.11-blue?logo=python&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-Containerized-2496ED?logo=docker&logoColor=white)
![PGGB](https://img.shields.io/badge/PGGB-Pangenome-green)
![vg](https://img.shields.io/badge/vg-v1.71.0-orange)
![Hackathon](https://img.shields.io/badge/CMU%20x%20NVIDIA-Hackathon%202026-76B900?logo=nvidia&logoColor=white)

---

## 🧠 Problem & Motivation

Pangenomes offer a powerful alternative to single-reference genomes by representing genetic variation across many individuals. However, **current pangenome graphs typically include only tens of samples**, while modern biobanks (UK Biobank, All of Us, HPRC) now contain **thousands of long-read genomes** spanning diverse populations.

A major barrier remains:

> **There is no practical way for independent biobanks to collaboratively construct population-scale pangenomes without sharing raw genomic data.**

This leads to:
- Silos across cohorts  
- Persistent population bias  
- Limited downstream interpretability  

---

## 🎯 Project Vision

**OmniGenome** explores whether **federated learning principles** can be applied to pangenome construction and genetic association analysis.

Instead of sharing raw FASTA files, participating sites:
- Build **local pangenome graphs**
- Share **graph representations only**
- Aggregate information *per chromosome*

We demonstrate this idea through:
1. **Federated pangenome graph construction** (HPRC)
2. **Genomic background hashing** for phenotype association at the **APOE locus**

---

## 🚀 Quick Start (Demo)

This quick start runs a **single local PGGB job** to demonstrate the core building block used throughout the project.

### Prerequisites
- Docker (with Compose)
- ~8 GB RAM
- bgzip-compatible FASTA input

### Run Example

```bash
# Start Docker
systemctl start docker

# Build containers
docker compose build

# Run PGGB on an example FASTA
docker compose run pggb pggb \
  -i /data/input.fa.gz \
  -o /output/my_graph \
  -n 12 \
  -t 8 -p 90 -s 10000
```

> ⚠️ FASTA files **must be bgzip-compressed**, not gzip.

This produces a local **GFA pangenome graph**, which is the unit shared in the federated workflow.

---

## 🧬 Project Overview

### Two Complementary Components

#### 1️⃣ Federated Pangenome Construction  
Simulates how multiple biobanks can collaboratively build pangenomes while keeping raw genomic data local.

#### 2️⃣ Genomic Background Hashing (APOE Case Study)  
Uses localized pangenome context to encode anonymized haploblock structure for federated association analysis.

Together, these form an end-to-end framework for **privacy-preserving, population-scale pangenomics**.

---

## 🧬 High-Level Architecture

```mermaid
flowchart LR

  %% =========================
  %% INPUT FASTA CHUNKS
  %% =========================
  subgraph Input["Chromosome 19 FASTA Chunks"]
    direction TB
    C1["<b>Chunk 1</b><br/>(37 individuals)"]
    C2["<b>Chunk 2</b><br/>(38 individuals)"]
    C3["<b>Chunk 3</b><br/>(32 individuals)"]
    C4["<b>Chunk 4</b><br/>(37 individuals)"]
    C5["<b>Chunk 5</b><br/>(28 individuals)"]
  end

  %% =========================
  %% FEDERATED SITES
  %% =========================
  subgraph Sites["Federated Sites (Subsampled)"]
    direction TB
    S1["<b>Site 1</b><br/>(n = 20)"]
    S2["<b>Site 2</b><br/>(n = 20)"]
    S3["<b>Site 3</b><br/>(n = 20)"]
    S4["<b>Site 4</b><br/>(n = 20)"]
    S5["<b>Site 5</b><br/>(n = 20)"]
  end

  %% =========================
  %% LOCAL GRAPH CONSTRUCTION
  %% =========================
  subgraph LocalGraphs["Local Graph Construction (PGGB)"]
    direction TB
    G1["<b>GFA 1</b>"]
    G2["<b>GFA 2</b>"]
    G3["<b>GFA 3</b>"]
    G4["<b>GFA 4</b>"]
    G5["<b>GFA 5</b>"]
  end

  %% =========================
  %% AGGREGATION
  %% =========================
  subgraph Aggregate["Federated Aggregation (Chromosome 19)"]
    direction TB
    Mega["<b>Megagraph</b><br/>(Chr19)"]
  end

  %% =========================
  %% FORWARD FLOW (LEFT → RIGHT)
  %% =========================
  C1 --> S1 -->|PGGB| G1 --> Mega
  C2 --> S2 -->|PGGB| G2 --> Mega
  C3 --> S3 -->|PGGB| G3 --> Mega
  C4 --> S4 -->|PGGB| G4 --> Mega
  C5 --> S5 -->|PGGB| G5 --> Mega

  %% =========================
  %% FEEDBACK LOOP
  %% =========================
  Mega -. "<b>minigraph feedback</b>" .-> G3
  Mega -. "<b>minigraph feedback</b>" .-> G1
  Mega -. "<b>minigraph feedback</b>" .->  G2
  Mega -. "<b>minigraph feedback</b>" .->  G4
  Mega -. "<b>minigraph feedback</b>" .->  G5

  %% =========================
  %% STYLING (HIGH CONTRAST)
  %% =========================
  classDef input fill:#e6e6e6,stroke:#aaaaaa,stroke-width:1px;
  classDef site  fill:#d6ebff,stroke:#3b82f6,stroke-width:1.5px;
  classDef local fill:#e6fff5,stroke:#10b981,stroke-width:2px;
  classDef mega  fill:#fff0cc,stroke:#f59e0b,stroke-width:3px;

  class C1,C2,C3,C4,C5 input
  class S1,S2,S3,S4,S5 site
  class G1,G2,G3,G4,G5 local
  class Mega mega

```

---

## 🔬 Component 1: Federated Pangenome Graph Construction

### Data
- **Human Pangenome Reference Consortium (HPRC)**
- Chromosomes **19** and **22** (tractable proof-of-principle)

### Core Idea

Each simulated “site”:
- Holds private FASTA data
- Builds a **local pangenome graph** using PGGB
- Shares only the resulting **GFA graph**

Graphs are then:
- Aggregated **per chromosome**
- Used to refine local graphs iteratively

> ⚠️ Graphs are **never combined across chromosomes**

---

### Minimal Pipeline Summary

1. Download HPRC assemblies  
2. Extract chr19 / chr22  
3. Subsample individuals to simulate biobank cohorts  
4. Run PGGB locally  
5. Aggregate graphs with `vg combine`  
6. (Optional) Feedback using `minigraph`

---

### Example: Local Graph Construction

```bash
docker compose run pggb pggb \
  -i /data/input.fa.gz \
  -o /output/run_name \
  -n 20 \
  -t 8 -p 90 -s 10000
```

---

## 🧬 Component 2: Genomic Background Hashing (APOE)

To demonstrate downstream utility, we focus on the **APOE locus**, a major genetic risk factor for Alzheimer’s disease.

### Motivation

GWAS hits often lack genomic context:
- Identical risk alleles can appear in different haplotypic backgrounds
- These backgrounds may influence penetrance or downstream effects

Pangenomes provide a natural representation of this structure.

---

### Workflow

1. Identify APOE-associated loci from GWAS summary statistics  
2. Extract the corresponding pangenome subgraph  
3. Encode **anonymized haploblock structure** as hashes  
4. Use hashes as a federated genomic background feature  

This enables:
- Privacy-preserving association testing  
- Context-aware interpretation of risk alleles  
- Cross-cohort comparison without raw genotype sharing  

---

## 🧰 Key Tools

| Tool | Role |
|-----|-----|
| **PGGB** | Local pangenome construction |
| **vg** | Graph conversion & aggregation |
| **minigraph** | Incremental graph feedback |
| **odgi** | Graph slicing & analysis |
| **Docker** | Reproducibility |

---

## 👥 Team

**CMU × NVIDIA Hackathon 2026**

- Rob Loughnan  
- Adam Kehl  
- Jedrzej Kubica  
- Kumar Koushik Telaprolu  
- **Jeff Winchell**  
- Sanjnaa Sridhar  

---

## 📚 References

- Wightman DP *et al*. *Nature Genetics* (2021)  
- Garrison E *et al*. *Nature Methods* (2024)  
- Garrison E *et al*. *Nature Biotechnology* (2018)  
- Guarracino A *et al*. *Bioinformatics* (2022)  
- Sirén J *et al*. *Science* (2021)  
- Liao W-W *et al*. *Nature* (2023)  
- Li H *et al*. *Genome Biology* (2020)

---

## 📄 License

MIT License
