FROM rocker/r-ver:4.3.2

ARG BCFTOOLS_VERSION=1.23

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    libicu-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    ca-certificates \
    wget \
    xz-utils \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# bcftools
RUN wget -q https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && cd bcftools-${BCFTOOLS_VERSION} \
    && make \
    && make install \
    && cd / \
    && rm -rf /tmp/bcftools-${BCFTOOLS_VERSION} /tmp/bcftools-${BCFTOOLS_VERSION}.tar.bz2

WORKDIR /opt/mip-poolfreq

# Afreq CLI
COPY Afreq/ /opt/mip-poolfreq/Afreq/
RUN chmod +x /opt/mip-poolfreq/Afreq/bin/Afreq \
    /opt/mip-poolfreq/Afreq/scripts/compute_allelic_freqs.sh \
    && ln -sf /opt/mip-poolfreq/Afreq/bin/Afreq /usr/local/bin/Afreq


# pak
RUN R -q -e 'site <- "/usr/local/lib/R/site-library"; \
.libPaths(c(site, .libPaths())); \
install.packages("pak", repos="https://cloud.r-project.org", lib=site)'


# Gfreq dependencies
COPY Gfreq/DESCRIPTION Gfreq/DESCRIPTION
COPY Gfreq/NAMESPACE Gfreq/NAMESPACE

RUN R -q -e 'site <- "/usr/local/lib/R/site-library"; \
  .libPaths(c(site, .libPaths())); \
  pak::local_install_deps("Gfreq", dependencies=TRUE, upgrade=FALSE)'


# Gfreq CLI
COPY Gfreq/ Gfreq/
RUN R CMD INSTALL Gfreq \
    && ln -sf /usr/local/lib/R/site-library/Gfreq/bin/Gfreq /usr/local/bin/Gfreq

# Entrypoint
COPY docker/entrypoint.sh /usr/local/bin/container-help
RUN chmod +x /usr/local/bin/container-help

ENTRYPOINT ["/usr/local/bin/container-help"]
CMD ["--help"]
