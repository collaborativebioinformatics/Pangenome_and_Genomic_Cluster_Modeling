#!/bin/bash
# Extract APOE region from pangenome graph

set -e

APOE_REGION="chr19:45274138-45479708"
INPUT_GRAPH="$1"
OUTPUT_DIR="../results/graphs/APOE_test"

if [ -z "$INPUT_GRAPH" ]; then
    echo "Usage: $0 <input_graph.og>"
    echo "Example: $0 ../results/graphs/chr19_graph/chr19.og"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Extracting APOE region: $APOE_REGION"
echo "From graph: $INPUT_GRAPH"

# Extract APOE region using odgi
docker compose -f ../docker/docker-compose.yml run vg bash -c "odgi extract \
    -i /data/$(basename $INPUT_GRAPH) \
    -r $APOE_REGION \
    -o /output/APOE_test/apoe_extracted.og -P"

# Convert to GFA
docker compose -f ../docker/docker-compose.yml run vg bash -c "odgi view \
    -i /output/APOE_test/apoe_extracted.og \
    -g > /output/APOE_test/apoe_extracted.gfa"

# Create Giraffe index
docker compose -f ../docker/docker-compose.yml run vg autoindex --workflow giraffe \
    -g /output/APOE_test/apoe_extracted.gfa \
    -p /output/APOE_test/apoe_giraffe

echo "APOE graph extraction complete!"
echo "Output: $OUTPUT_DIR/"
