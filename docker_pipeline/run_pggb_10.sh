#!/bin/bash
SUBCHUNK_DIR="/mnt/shared_vol/hprc_mini_fasta/subchunks"
OUTPUT_DIR="/mnt/shared_vol/graphs"

echo "=== Starting PGGB (10 individuals) at $(date) ==="

for i in 1 2 3 4 5; do
  SUBCHUNK="chr19_chunk${i}_sub10.fa.gz"
  OUTPUT="${OUTPUT_DIR}/chr19_chunk${i}_sub10"
  
  echo ""
  echo "=== Processing ${SUBCHUNK} at $(date) ==="
  mkdir -p ${OUTPUT}
  
  N=$(zcat ${SUBCHUNK_DIR}/${SUBCHUNK} | grep -c '^>')
  echo "Sequences: $N"
  
  docker run --rm \
    -v ${SUBCHUNK_DIR}:/data \
    -v ${OUTPUT}:/output \
    ghcr.io/pangenome/pggb:latest \
    pggb -i /data/${SUBCHUNK} -o /output -n ${N} -t 8 -p 90 -s 10000
  
  echo "=== Done ${SUBCHUNK} at $(date) ==="
done

echo "=== All done at $(date) ==="
