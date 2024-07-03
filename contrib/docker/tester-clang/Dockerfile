FROM ubuntu:24.04

LABEL maintainer <rene.gassmoeller@mailbox.org>

RUN DEBIAN_FRONTEND=noninteractive apt update && apt upgrade -yq && \
  DEBIAN_FRONTEND=noninteractive apt install -yq --no-install-recommends \
  build-essential ca-certificates file gcc g++ \
  gfortran git libblas-dev liblapack-dev libopenmpi-dev \
  lsb-release ninja-build numdiff openmpi-bin \
  openmpi-common wget zlib1g-dev \
  clang-tidy clang cmake libp4est-dev libsundials-dev \
  mpi-default-dev trilinos-all-dev zlib1g-dev libhdf5-dev libhdf5-mpi-dev

COPY local.cfg /opt/

# Install deal.II with clang using candi:

RUN cd /opt && \
    git clone https://github.com/dealii/candi && \
    cd candi && \
    mv /opt/local.cfg . && \
    CXX=clang++ CC=clang ./candi.sh -p /opt -j4 && \
    rm -rf /opt/tmp

# Set environment variables for this image to be used
# by Github Actions
ENV PATH="/opt/astyle-2.04:$PATH"
ENV DEAL_II_DIR /opt/deal.II-master
ENV OMPI_MCA_btl_base_warn_component_unused=0
ENV OMPI_MCA_mpi_yield_when_idle=1
ENV OMPI_MCA_rmaps_base_oversubscribe=1
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /opt
