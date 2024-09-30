# pull base image
FROM rocker/tidyverse:4.2.0

# build image: docker build -t core_bioinf/dockrstudio:4.1.0-v1 -f Rstudio-4.1.0-v1.Dockerfile .

# image based on:
# https://rocker-project.org/
# https://github.com/rocker-org/rocker-versioned2

# who maintains this image
LABEL maintainer="Aleksandr Bykov"
LABEL version="4.2-java"

RUN apt-get update && \
      DEBIAN_FRONTEND=noninteractive \
      apt-get install --assume-yes \
      libglpk40 \
      libcairo2-dev \
      liblzma-dev \
      libxt-dev \
      libcurl4-openssl-dev \
      libcurl4 \
      libxml2-dev \
      libxt-dev \
      openssl \
      libssl-dev \
      wget \
      curl \
      bzip2 \
      libbz2-dev \
      libpng-dev \
      libhdf5-dev \
      pigz \
      libmagick++-dev

RUN apt-get install -y openjdk-11-jdk

RUN export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
RUN apt-get install -y r-cran-rjava
RUN R CMD javareconf
# commonly used R packages 
# installation from CRAN
RUN install2.r --error Cairo XML RCurl R.utils digest optparse

# installation from Github repositiories:
#RUN installGithub.r kassambara/ggpubr

# # installation from Bioconductor:
RUN R -e "BiocManager::install(c('BiocParallel', 'Biostrings'))"

RUN chmod -R a+rw /usr/local/lib/R/site-library # so that everyone can dynamically install more libraries within container

# # Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# If not done above - clean temp dir after install2.r command [https://rocker-project.org/use/extending.html]
# RUN rm -rf /tmp/downloaded_packages

# should already be exposed by rstudio (on top of which tidyverse is build)?
#EXPOSE 8787 

#CMD ["R"]
#CMD ["/init"] # needed only when installing from scratch; e.g. r-ver; and rstudio, shiny

