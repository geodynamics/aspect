FROM ubuntu:20.04

LABEL maintainer <rene.gassmoeller@mailbox.org>

RUN DEBIAN_FRONTEND=noninteractive apt update && apt upgrade -yq && \
  DEBIAN_FRONTEND=noninteractive apt install -yq --no-install-recommends \
  build-essential ca-certificates cmake cmake-curses-gui file gcc g++ \
  gfortran git libblas-dev liblapack-dev libopenmpi-dev \
  lsb-release ninja-build numdiff openmpi-bin \
  openmpi-common wget zlib1g-dev \
  texlive-plain-generic texlive-base texlive-latex-recommended \
  texlive-latex-base texlive-fonts-recommended \
  texlive-bibtex-extra lmodern texlive-latex-extra \
  texlive-science graphviz python3-pip python-setuptools

RUN pip3 install cpp-coveralls

COPY local.cfg /opt/

# Install deal.II with minimal dependencies for ASPECT
RUN cd /opt && \
    git clone https://github.com/dealii/candi && \
    cd candi && \
    mv /opt/local.cfg . && \
    git checkout v9.3.1-r1 && \
    ./candi.sh -p /opt -j64 && \
    rm -rf /opt/tmp && \
    ln -s /opt/astyle-2.04/astyle /usr/bin/astyle

# Set environment variables for this image to be used
# by Github Actions
ENV DEAL_II_DIR /opt/deal.II-v9.3.1
ENV OMPI_MCA_btl_base_warn_component_unused=0
ENV OMPI_MCA_mpi_yield_when_idle=1
ENV OMPI_MCA_rmaps_base_oversubscribe=1
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Set environment variables for derived images, which
# are build by Jenkins
ONBUILD ENV DEAL_II_DIR /opt/deal.II-v9.3.1
ONBUILD ENV OMPI_MCA_btl_base_warn_component_unused=0
ONBUILD ENV OMPI_MCA_mpi_yield_when_idle=1
ONBUILD ENV OMPI_MCA_rmaps_base_oversubscribe=1
ONBUILD ENV OMPI_ALLOW_RUN_AS_ROOT=1
ONBUILD ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /opt
