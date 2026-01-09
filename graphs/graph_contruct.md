═══════════════════════════════════════════════════════════════════════════════
                    FEDERATED PANGENOME CONSTRUCTION WORKFLOW
═══════════════════════════════════════════════════════════════════════════════

INPUT: cleaned chr19 chunks (each chunk = ~30-38 individuals)
───────────────────────────────────────────────────────────────────────────────

/mnt/shared_vol/hprc_mini_fasta/
├── chrom19_chunk1.fa.gz (37 individuals)
├── chrom19_chunk2.fa.gz (38 individuals)
├── chrom19_chunk3.fa.gz (32 individuals)
├── chrom19_chunk4.fa.gz (37 individuals)
└── chrom19_chunk5.fa.gz (28 individuals)


═══════════════════════════════════════════════════════════════════════════════
STEP 0: SUBCHUNKING (Take first 20 individuals from each)
═══════════════════════════════════════════════════════════════════════════════

chrom19_chunk1.fa.gz (37 ind) ──take first 20──► chunk1_sub.fa.gz (20 ind)
chrom19_chunk2.fa.gz (38 ind) ──take first 20──► chunk2_sub.fa.gz (20 ind)
chrom19_chunk3.fa.gz (32 ind) ──take first 20──► chunk3_sub.fa.gz (20 ind)
chrom19_chunk4.fa.gz (37 ind) ──take first 20──► chunk4_sub.fa.gz (20 ind)
chrom19_chunk5.fa.gz (28 ind) ──take first 20──► chunk5_sub.fa.gz (20 ind)
                                                      │
                                                      ▼
                                    /mnt/shared_vol/hprc_mini_fasta/subchunks/


═══════════════════════════════════════════════════════════════════════════════
STEP 1: LOCAL GRAPH CONSTRUCTION (PGGB) - Each site builds its own graph
═══════════════════════════════════════════════════════════════════════════════

  SITE 1              SITE 2              SITE 3              SITE 4              SITE 5
     │                   │                   │                   │                   │
     ▼                   ▼                   ▼                   ▼                   ▼
┌──────────┐       ┌──────────┐       ┌──────────┐       ┌──────────┐       ┌──────────┐
│ chunk1   │       │ chunk2   │       │ chunk3   │       │ chunk4   │       │ chunk5   │
│ _sub.fa  │       │ _sub.fa  │       │ _sub.fa  │       │ _sub.fa  │       │ _sub.fa  │
│ (20 ind) │       │ (20 ind) │       │ (20 ind) │       │ (20 ind) │       │ (20 ind) │
└────┬─────┘       └────┬─────┘       └────┬─────┘       └────┬─────┘       └────┬─────┘
     │ PGGB             │ PGGB             │ PGGB             │ PGGB             │ PGGB
     ▼                   ▼                   ▼                   ▼                   ▼
┌──────────┐       ┌──────────┐       ┌──────────┐       ┌──────────┐       ┌──────────┐
│   GFA1   │       │   GFA2   │       │   GFA3   │       │   GFA4   │       │   GFA5   │
└──────────┘       └──────────┘       └──────────┘       └──────────┘       └──────────┘

                              OUTPUT: /mnt/shared_vol/graphs/chr19_chunkN/


═══════════════════════════════════════════════════════════════════════════════
STEP 2: AGGREGATION (vg combine) - Central server merges all graphs
═══════════════════════════════════════════════════════════════════════════════

    GFA1          GFA2          GFA3          GFA4          GFA5
      │             │             │             │             │
      │  vg convert │  vg convert │  vg convert │  vg convert │  vg convert
      ▼             ▼             ▼             ▼             ▼
    VG1           VG2           VG3           VG4           VG5
      │             │             │             │             │
      └─────────────┴──────┬──────┴─────────────┴─────────────┘
                           │
                           │  vg combine -p
                           ▼
                   ┌───────────────┐
                   │ MEGAGRAPH.vg  │
                   └───────────────┘
                           │
                           │  vg convert -f (VG → GFA)
                           ▼
                   ┌───────────────┐
                   │ MEGAGRAPH.gfa │
                   └───────────────┘


═══════════════════════════════════════════════════════════════════════════════
STEP 3: FEEDBACK (minigraph) - Improve local graphs with global knowledge
═══════════════════════════════════════════════════════════════════════════════

                         MEGAGRAPH.gfa
                              │
        ┌───────────┬─────────┼─────────┬───────────┐
        │           │         │         │           │
        ▼           ▼         ▼         ▼           ▼
    chunk1_sub  chunk2_sub chunk3_sub chunk4_sub chunk5_sub
    .fa.gz      .fa.gz     .fa.gz     .fa.gz     .fa.gz
        │           │         │         │           │
        │ minigraph │         │         │           │
        ▼           ▼         ▼         ▼           ▼
    ┌───────┐   ┌───────┐ ┌───────┐ ┌───────┐   ┌───────┐
    │ GFA1' │   │ GFA2' │ │ GFA3' │ │ GFA4' │   │ GFA5' │
    └───────┘   └───────┘ └───────┘ └───────┘   └───────┘
    improved    improved  improved  improved    improved

    Command: minigraph -cxggs MEGAGRAPH.gfa chunkN_sub.fa.gz > GFAN_prime.gfa


═══════════════════════════════════════════════════════════════════════════════
FINAL OUTPUT STRUCTURE
═══════════════════════════════════════════════════════════════════════════════

/mnt/shared_vol/
├── hprc_mini_fasta/
│   ├── chrom19_chunk1.fa.gz         (original - 37 ind)
│   ├── chrom19_chunk2.fa.gz         (original - 38 ind)
│   ├── chrom19_chunk3.fa.gz         (original - 32 ind)
│   ├── chrom19_chunk4.fa.gz         (original - 37 ind)
│   ├── chrom19_chunk5.fa.gz         (original - 28 ind)
│   │
│   └── subchunks/                   ◄── STEP 0 OUTPUT
│       ├── chunk1_sub.fa.gz         (20 ind)
│       ├── chunk2_sub.fa.gz         (20 ind)
│       ├── chunk3_sub.fa.gz         (20 ind)
│       ├── chunk4_sub.fa.gz         (20 ind)
│       └── chunk5_sub.fa.gz         (20 ind)
│
└── graphs/
    ├── chr19_chunk1/*.gfa           ◄── STEP 1 OUTPUT (GFA1)
    ├── chr19_chunk2/*.gfa           (GFA2)
    ├── chr19_chunk3/*.gfa           (GFA3)
    ├── chr19_chunk4/*.gfa           (GFA4)
    ├── chr19_chunk5/*.gfa           (GFA5)
    │
    ├── MEGAGRAPH.vg                 ◄── STEP 2 OUTPUT
    ├── MEGAGRAPH.gfa
    │
    └── federated/                   ◄── STEP 3 OUTPUT
        ├── chunk1_federated.gfa     (GFA1' - improved)
        ├── chunk2_federated.gfa     (GFA2' - improved)
        ├── chunk3_federated.gfa     (GFA3' - improved)
        ├── chunk4_federated.gfa     (GFA4' - improved)
        └── chunk5_federated.gfa     (GFA5' - improved)


═══════════════════════════════════════════════════════════════════════════════