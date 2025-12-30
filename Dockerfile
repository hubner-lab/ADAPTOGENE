# Base image with R 4.2 pre-installed
FROM rocker/shiny:4.2

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
# Plink install
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
#RUN pipx install snakemake pulp==2.7
RUN pipx install snakemake && pipx inject snakemake pulp==2.7
#ENV PATH="/root/.local/bin:$PATH" 
# to make it accesable via cli (due to pipx)


# Copy Snakemake workflow files
#TODO temporary with absolute path
#COPY /mnt/data/eugene/PopGenPipe/Snakefile /pipeline/Snakefile
#COPY /mnt/data/eugene/PopGenPipe/config.yaml /pipeline/config.yaml
#COPY /mnt/data/eugene/PopGenPipe/scripts/ /pipeline/scripts/
#TODO I don't think we need this in the future
#COPY /mnt/data/eugene/PopGenPipe/data/ /pipeline/data/

# Create important dirs #TODO I think it's no needed to create it in container, it should be created by Snakefile
#RUN mkdir -p "/pipeline/results/plots" "/pipeline/results/tables" "/pipeline/logs"
# Set working directory
WORKDIR /pipeline

# Make scripts executable
#RUN chmod +x /pipeline/scripts/*.sh

# Install R packages
RUN Rscript -e "install.packages(c('ggplot2', 'dplyr', 'data.table', 'BiocManager', 'viridis', 'gridExtra', 'remotes', 'scatterpie', 'geodata', 'raster', 'sp', 'geosphere', 'qs', 'ggnewscale', 'ggspatial', 'see', 'ggpubr', 'forcats', 'CMplot', 'ggcorrplot', 'ggplotify', 'shinydashboard', 'DT', 'topr'))"
RUN Rscript -e "BiocManager::install(c('LEA', 'WGCNA', 'vegan', 'qvalue', 'GenomicRanges', 'poppr', 'vcfR', 'adegenet'))"
# Install cli tools
RUN set -ex \
    && wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
    && tar zxf vcftools-0.1.16.tar.gz \
    && cd vcftools-0.1.16 \
    && ./configure --prefix=/app/vcftools-0.1.16 \
    && make \
    && make install

RUN Rscript -e "install.packages('svglite')"
RUN Rscript -e "remotes::install_github('r-forge/gradientforest/pkg/extendedForest')"
RUN Rscript -e "remotes::install_github('r-forge/gradientforest/pkg/gradientForest')"
#RUN Rscript -e "remotes::install_github('slarge/gradientForest')"
# Define the default command to execute the main pipeline script
#CMD ["/bin/bash", "/pipeline/scripts/main.sh"]

# Default entrypoint to run Snakemake
#ENTRYPOINT ["snakemake"]
