# Base image with R 4.5 (requires Bioconductor 3.22)
FROM rocker/shiny:4.5

# Permissions
RUN mkdir -p /.cache && chmod -R 777 /.cache

# Metadata
LABEL maintainer="potapgene@gmail.com"

# Update the system and install CLI tools
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    pipx \
    bash \
    wget \
    curl \
    git \
    gfortran \
    build-essential \
    zlib1g-dev \
    pkg-config \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libglpk40 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install plink1.9
ENV PLINK_VERSION=20231211
RUN mkdir -p /tmp/plink && \
  cd /tmp/plink && \
  curl -fsSL -o plink.zip "http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_${PLINK_VERSION}.zip" && \
  unzip plink.zip && \
  mv plink /bin/plink && \
  cd $HOME && \
  rm -rf /tmp/plink

# Configure pipx and Install Snakemake
ENV PIPX_HOME=/usr/local/pipx
ENV PIPX_BIN_DIR=/usr/local/bin
RUN pipx install snakemake && pipx inject snakemake pulp==2.7

# Set working directory
WORKDIR /pipeline

# Install vcftools
RUN set -ex \
    && wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
    && tar zxf vcftools-0.1.16.tar.gz \
    && cd vcftools-0.1.16 \
    && ./configure --prefix=/app/vcftools-0.1.16 \
    && make \
    && make install \
    && cd .. && rm -rf vcftools-0.1.16*

#=============================================================================
# R PACKAGES WITH PINNED VERSIONS
# Last verified: 2025-01-25
# Bioconductor: 3.22 (for R 4.5)
#=============================================================================

# Install remotes first for version-pinned installations
RUN Rscript -e "install.packages('remotes')"

# Core tidyverse and data manipulation (pinned versions)
RUN Rscript -e " \
    remotes::install_version('ggplot2', version = '3.5.1'); \
    remotes::install_version('dplyr', version = '1.1.4'); \
    remotes::install_version('tidyr', version = '1.3.1'); \
    remotes::install_version('tibble', version = '3.2.1'); \
    remotes::install_version('purrr', version = '1.0.2'); \
    remotes::install_version('stringr', version = '1.5.1'); \
    remotes::install_version('forcats', version = '1.0.0'); \
    remotes::install_version('tidyverse', version = '2.0.0'); \
"

# Data handling packages
RUN Rscript -e " \
    remotes::install_version('data.table', version = '1.16.4'); \
    remotes::install_version('reshape2', version = '1.4.4'); \
    remotes::install_version('qs', version = '0.27.2'); \
"

# Visualization packages
RUN Rscript -e " \
    remotes::install_version('viridis', version = '0.6.5'); \
    remotes::install_version('gridExtra', version = '2.3'); \
    remotes::install_version('scales', version = '1.3.0'); \
    remotes::install_version('cowplot', version = '1.1.3'); \
    remotes::install_version('egg', version = '0.4.5'); \
    remotes::install_version('ggpubr', version = '0.6.0'); \
    remotes::install_version('ggcorrplot', version = '0.1.4.1'); \
    remotes::install_version('ggplotify', version = '0.1.2'); \
    remotes::install_version('see', version = '0.9.0'); \
    remotes::install_version('svglite', version = '2.1.3'); \
    remotes::install_version('CMplot', version = '4.5.1'); \
    remotes::install_version('scattermore', version = '1.2'); \
"

# Spatial and geographic packages
RUN Rscript -e " \
    remotes::install_version('sp', version = '2.1-4'); \
    remotes::install_version('raster', version = '3.6-30'); \
    remotes::install_version('geodata', version = '0.6-2'); \
    remotes::install_version('geosphere', version = '1.5-20'); \
    remotes::install_version('scatterpie', version = '0.2.4'); \
    remotes::install_version('ggnewscale', version = '0.5.0'); \
    remotes::install_version('ggspatial', version = '1.1.9'); \
"

# Shiny packages
RUN Rscript -e " \
    remotes::install_version('shinydashboard', version = '0.7.2'); \
    remotes::install_version('DT', version = '0.33'); \
"

# topr - CRITICAL: version >= 2.0.0 required for custom (non-human) genome builds
RUN Rscript -e "remotes::install_version('topr', version = '2.0.2')"

# BiocManager for Bioconductor packages
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install(version = '3.22', ask = FALSE)"

# Bioconductor packages (pinned to 3.22 versions)
RUN Rscript -e "BiocManager::install('LEA', version = '3.22', ask = FALSE)"
RUN Rscript -e "BiocManager::install('qvalue', version = '3.22', ask = FALSE)"
RUN Rscript -e "BiocManager::install('GenomicRanges', version = '3.22', ask = FALSE)"
RUN Rscript -e "BiocManager::install('WGCNA', version = '3.22', ask = FALSE)"
RUN Rscript -e "BiocManager::install('clusterProfiler', version = '3.22', ask = FALSE)"
RUN Rscript -e "BiocManager::install('AnnotationDbi', version = '3.22', ask = FALSE)"
RUN Rscript -e "BiocManager::install('GO.db', version = '3.22', ask = FALSE)"

# CRAN packages used with Bioconductor workflows
RUN Rscript -e " \
    remotes::install_version('vegan', version = '2.6-8'); \
    remotes::install_version('vcfR', version = '1.15.0'); \
    remotes::install_version('adegenet', version = '2.1.10'); \
    remotes::install_version('poppr', version = '2.9.6'); \
"

# gradientForest from r-forge (version 0.1-37)
# Requires patching for R 4.5 (Calloc/Free -> R_Calloc/R_Free)
RUN git clone --depth 1 https://github.com/r-forge/gradientforest.git /tmp/gf && \
    cd /tmp/gf/pkg/extendedForest/src && \
    sed -i 's/\bCalloc(/R_Calloc(/g; s/\bFree(/R_Free(/g' *.c && \
    cd /tmp/gf/pkg/gradientForest/src && \
    sed -i 's/\bCalloc(/R_Calloc(/g; s/\bFree(/R_Free(/g' *.c 2>/dev/null || true && \
    Rscript -e "install.packages('/tmp/gf/pkg/extendedForest', repos = NULL, type = 'source')" && \
    Rscript -e "install.packages('/tmp/gf/pkg/gradientForest', repos = NULL, type = 'source')" && \
    rm -rf /tmp/gf

# Verify critical package versions
RUN Rscript -e " \
    stopifnot(packageVersion('topr') >= '2.0.0'); \
    stopifnot(packageVersion('ggplot2') >= '3.5.0'); \
    stopifnot(requireNamespace('gradientForest', quietly = TRUE)); \
    cat('All package version checks passed.\n'); \
"
