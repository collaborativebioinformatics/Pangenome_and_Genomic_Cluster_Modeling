#!/bin/bash

SUBCHUNK_DIR="/mnt/shared_vol/hprc_mini_fasta/subchunks"
OUTPUT_DIR="/mnt/shared_vol/graphs"

echo "=== Starting PGGB pipeline at $(date) ==="

for i in 1 2 3 4 5; do
  SUBCHUNK="chr19_chunk${i}_sub20.fa.gz"
  OUTPUT="${OUTPUT_DIR}/chr19_chunk${i}_sub20"
  
  echo ""
  echo "=== Processing ${SUBCHUNK} at $(date) ==="
  mkdir -p ${OUTPUT}
  
  # Get sequence count
  N=$(zcat ${SUBCHUNK_DIR}/${SUBCHUNK} | grep -c '^>')
  echo "Sequences: $N"
  
  # Run PGGB
  docker run --rm \
    -v ${SUBCHUNK_DIR}:/data \
    -v ${OUTPUT}:/output \
    ghcr.io/pangenome/pggb:latest \
    pggb -i /data/${SUBCHUNK} -o /output -n ${N} -t 8 -p 90 -s 10000
  
  echo "=== Done ${SUBCHUNK} at $(date) ==="
done

echo ""
echo "=== All PGGB jobs completed at $(date) ==="
echo "=== Now combining graphs ==="

# Step 2: Combine graphs
echo "Finding GFA files..."
GFA_FILES=$(find ${OUTPUT_DIR}/chr19_chunk*_sub20 -name "*.smooth.final.gfa" | sort)
echo "Found: $GFA_FILES"

echo "=== Pipeline finished at $(date) ==="
