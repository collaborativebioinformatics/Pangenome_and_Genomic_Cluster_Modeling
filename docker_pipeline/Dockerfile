FROM python:3.11-slim

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    gzip \
    docker.io \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# Install minigraph
RUN wget https://github.com/lh3/minigraph/releases/download/v0.21/minigraph-0.21_x64-linux.tar.bz2 \
    && tar -xjf minigraph-0.21_x64-linux.tar.bz2 \
    && mv minigraph-0.21_x64-linux/minigraph /usr/local/bin/ \
    && chmod +x /usr/local/bin/minigraph \
    && rm -rf minigraph-0.21_x64-linux*

# Create working directory
WORKDIR /pipeline

# Copy pipeline script
COPY federated_pangenome_pipeline.py /pipeline/

# Make executable
RUN chmod +x /pipeline/federated_pangenome_pipeline.py

# Entry point
ENTRYPOINT ["python3", "/pipeline/federated_pangenome_pipeline.py"]
