FROM tacc/tacc-ubuntu18-impi19.0.7-common:latest

LABEL maintainer <rene.gassmoeller@mailbox.org>

RUN apt-get update && apt-get upgrade -y

# we need a newer version of cmake to support unity builds:
RUN apt-get update && apt-get -yq upgrade && apt-get -yq install \
    libblas-dev \
    liblapack-dev \
    git \
    wget \
    cmake \
    python \
    numdiff \
    bc \
    ninja-build

RUN cd /opt && wget https://github.com/Kitware/CMake/releases/download/v3.17.3/cmake-3.17.3-Linux-x86_64.tar.gz && tar xf cmake*.tar.gz && rm cmake*.tar.gz
ENV PATH /opt/cmake-3.17.3-Linux-x86_64/bin:$PATH

RUN cd /opt && \
    git clone https://github.com/dealii/candi.git && \
    cd candi && \
    echo "NATIVE_OPTIMIZATIONS=true" >> local.cfg && \
    echo "STABLE_BUILD=false" >> local.cfg && \
    echo "DEAL_II_VERSION=master" >> local.cfg && \
    echo "DEAL_CONFOPTS=-D DEAL_II_WITH_CXX17=OFF -D DEAL_II_WITH_64BIT_INDICES=ON -D DEAL_II_COMPONENT_EXAMPLES=OFF -DDEAL_II_WITH_COMPLEX_VALUES=OFF -DDEAL_II_WITH_LAPACK=ON" >> local.cfg && \
    ./candi.sh -j 64 --packages="once:p4est once:trilinos once:hdf5 dealii" --platform=deal.II-toolchain/platforms/supported/ubuntu18.platform -p /opt && \
    rm -rf /opt/tmp

# Build aspect, replace git checkout command to create image for release
RUN cd /opt && \
    export DEAL_II_DIR=/opt/deal.II-master && \
    git clone https://github.com/geodynamics/aspect.git ./aspect && \ 
    mkdir aspect/build-release && \
    cd aspect/build-release && \
    git checkout master && \
    /opt/cmake-3.17.3-Linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release \
          .. && \
    make -j64 && \
    make install && \
    make clean

RUN cd /opt/aspect && \
    mkdir build-debug && \
    cd build-debug && \
    /opt/cmake-3.17.3-Linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Debug \
          -DCMAKE_INSTALL_PREFIX=/opt/aspect/debug .. && \
    make -j64 && \
    make install && \
    ln -s /opt/aspect/debug/bin/aspect /usr/local/bin/aspect.debug && \
    make clean
