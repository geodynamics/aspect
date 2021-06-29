FROM ubuntu:20.04

LABEL maintainer <rene.gassmoeller@mailbox.org>

USER root

RUN DEBIAN_FRONTEND=noninteractive apt update && apt upgrade -yq && \
  DEBIAN_FRONTEND=noninteractive apt install -yq --no-install-recommends \
  build-essential ca-certificates file gcc g++ \
  gfortran git libblas-dev liblapack-dev libopenmpi-dev \
  lsb-release ninja-build numdiff openmpi-bin \
  openmpi-common wget zlib1g-dev \
  texlive-generic-extra texlive-base texlive-latex-recommended \
  texlive-latex-base texlive-fonts-recommended \
  texlive-bibtex-extra lmodern texlive-latex-extra \
  texlive-science graphviz python-pip python-setuptools

RUN pip install cpp-coveralls

# we need a newer version of cmake to support unity builds:
RUN cd $HOME && wget https://github.com/Kitware/CMake/releases/download/v3.17.3/cmake-3.17.3-Linux-x86_64.tar.gz && tar xf cmake*.tar.gz && rm cmake*.tar.gz
ENV PATH="/root/cmake-3.17.3-Linux-x86_64/bin:$PATH"

COPY local.cfg /root/

# Install deal.II with minimal dependencies for ASPECT
RUN cd $HOME && \
    git clone https://github.com/dealii/candi && \
    cd candi && \
    mv /root/local.cfg . && \
    git checkout dealii-9.2 && \
    ./candi.sh --packages="once:hdf5 once:p4est once:trilinos once:dealii" -p $HOME -j64 && \
    rm -rf $HOME/tmp

ENV DEAL_II_DIR /root/deal.II-v9.2.0

ENV OMPI_MCA_btl_base_warn_component_unused=0
ENV OMPI_MCA_mpi_yield_when_idle=1
ENV OMPI_MCA_rmaps_base_oversubscribe=1
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /root
