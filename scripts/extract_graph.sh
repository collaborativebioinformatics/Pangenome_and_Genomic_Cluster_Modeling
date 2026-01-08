#!/bin/bash
# Extract APOE region from pangenome graph

set -e

APOE_REGION="chr19:45274138-45479708"
INPUT_GRAPH="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_DIR="${2:-$ROOT_DIR/results/graphs/APOE_test}"

if [ -z "$INPUT_GRAPH" ]; then
    echo "Usage: $0 <input_graph.og> [output_dir]"
    echo "Example: $0 $ROOT_DIR/results/graphs/chr19_graph/chr19.og $ROOT_DIR/results/graphs/APOE_test"
    exit 1
fi

if [ ! -f "$INPUT_GRAPH" ]; then
    echo "Error: Input graph file not found: $INPUT_GRAPH"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Extracting APOE region: $APOE_REGION"
echo "From graph: $INPUT_GRAPH"
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Note: Ensure that the docker-compose.yml mounts the directory containing"
echo "      the input graph to /data/ inside the container."
echo ""

# Extract APOE region using odgi
docker compose -f ../docker/docker-compose.yml run vg bash -c "odgi extract \
    -i /data/$(basename "$INPUT_GRAPH") \
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
