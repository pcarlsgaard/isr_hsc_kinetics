FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential pkg-config git-lfs jq curl ca-certificates \
    r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    libhdf5-dev libhts-dev zlib1g-dev \
    samtools bcftools bedtools tabix fastqc \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /home/vscode/.claude /home/vscode/R/library \
    && chown -R vscode:vscode /home/vscode/.claude /home/vscode/R
ENV R_LIBS_USER=/home/vscode/R/library
USER vscode
