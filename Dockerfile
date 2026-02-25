FROM rocker/r-ver:4.3.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    libicu-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/mip-poolfreq

# 1) Install pak into the site library
RUN R -q -e 'site <- "/usr/local/lib/R/site-library"; \
  .libPaths(c(site, .libPaths())); \
  install.packages("pak", repos="https://cloud.r-project.org", lib=site)'

# 2) Copy only metadata first (cache-friendly)
COPY Gfreq/DESCRIPTION Gfreq/DESCRIPTION
COPY Gfreq/NAMESPACE  Gfreq/NAMESPACE

# 3) Install dependencies declared in DESCRIPTION
RUN R -q -e 'site <- "/usr/local/lib/R/site-library"; .libPaths(c(site, .libPaths())); pak::local_install_deps("Gfreq", dependencies=TRUE, upgrade=FALSE)'

# 4) Copy the rest of the package and install it
COPY Gfreq/ Gfreq/
RUN R -q -e 'site <- "/usr/local/lib/R/site-library"; .libPaths(c(site, .libPaths())); R.home();' \
 && R CMD INSTALL Gfreq

RUN ln -s /usr/local/lib/R/site-library/Gfreq/bin/Gfreq /usr/local/bin/Gfreq

ENTRYPOINT ["Gfreq"]
CMD ["--help"]
