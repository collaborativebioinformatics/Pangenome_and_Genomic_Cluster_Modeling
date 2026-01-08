# APOE Graph Extraction Status

## ✅ Completed
1. **GWAS Data Downloaded** - GCST013197 Alzheimer's dataset (12.7M SNPs)
2. **APOE Locus Identified** - chr19:45,274,138-45,479,708 (206 kb)
3. **BED File Created** - `gwas/data/apoe_locus.bed` for odgi extract
4. **Coordinates Extracted** - `gwas/data/apoe_locus_interval.tsv`
5. **Visualization Generated** - `gwas/data/apoe_locus_plot.png`
6. **Extraction Script Ready** - `scripts/extract_graph.sh`

## ❌ Blocked - Needs HPRC Data
Cannot complete graph extraction without:
- HPRC chr19 assemblies (multi-sample FASTA)
- Or existing chr19 pangenome graph file

## To Complete:

### Option 1: Download HPRC Data
```bash
# Download HPRC assemblies for chr19
# Requires HPRC download scripts from:
# https://github.com/human-pangenomics/HPP_Year1_Assemblies
```

### Option 2: Use Existing Graph
If you have an existing chr19 graph:
```bash
cd scripts
./extract_graph.sh /path/to/chr19_graph.og
```

### Option 3: Build Test Graph
With any chr19 FASTA data:
```bash
cd docker
docker compose run pggb pggb \
    -i /data/input.fa.gz \
    -o /output/APOE_test \
    -n <haplotypes> -t 8 -p 90 -s 10000
```

## Summary
**Task Status**: Partially complete - GWAS extraction ✅, Graph extraction ❌ (blocked on data)
**Blocker**: Need HPRC chr19 assembly files or existing pangenome graph
**Ready**: All tools, scripts, and coordinates prepared for extraction when data is available
