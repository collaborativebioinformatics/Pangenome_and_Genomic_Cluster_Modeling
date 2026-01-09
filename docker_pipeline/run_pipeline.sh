#!/bin/bash
echo "=========================================="
echo "Federated Pangenome Pipeline Runner"
echo "=========================================="
echo ""

# Build Docker image
echo "Building Docker image..."
docker build -t federated_pangenome_pipeline .

echo ""
echo "Starting pipeline in background..."

# Run container in detached mode
docker run -d \
    --name federated_pangenome_pipeline \
    -v /mnt/shared_vol:/mnt/shared_vol \
    -v /var/run/docker.sock:/var/run/docker.sock \
    federated_pangenome_pipeline

echo ""
echo "Pipeline started!"
echo ""
echo "Monitor logs with:"
echo "  docker logs -f federated_pangenome_pipeline"
echo ""
echo "Or check log file:"
echo "  tail -f /mnt/shared_vol/graphs/pipeline_*.log"
echo ""
echo "Check status:"
echo "  docker ps | grep federated"
