FROM tacc/tacc-ubuntu18-impi19.0.7-common:latest

LABEL maintainer <rene.gassmoeller@mailbox.org>

# we need a newer version of cmake to support unity builds:
RUN apt-get update && apt-get -yq upgrade && apt-get -yq install \
    libblas-dev \
    liblapack-dev \
    git \
    wget \
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
    echo "DEAL_II_VERSION=v9.4.2" >> local.cfg && \
    echo "DEAL_CONFOPTS=-D DEAL_II_WITH_64BIT_INDICES=ON -D DEAL_II_COMPONENT_EXAMPLES=OFF -DDEAL_II_WITH_COMPLEX_VALUES=OFF -DDEAL_II_WITH_LAPACK=ON" >> local.cfg && \
    ./candi.sh -j 2 --packages="once:cmake once:p4est once:trilinos once:hdf5 dealii" -p /opt && \
    rm -rf /opt/tmp

# Build aspect, replace git checkout command to create image for release
ENV Aspect_DIR /opt/aspect
RUN source /opt/configuration/enable.sh && \
    git clone https://github.com/geodynamics/aspect.git $Aspect_DIR && \
    mkdir ${Aspect_DIR}/build && \
    cd ${Aspect_DIR}/build && \
    cmake -DCMAKE_BUILD_TYPE=DebugRelease -DCMAKE_INSTALL_PREFIX=$Aspect_DIR \
    -DASPECT_INSTALL_EXAMPLES=ON .. && \
    make -j2 && make install && make clean

ENV PATH ${Aspect_DIR}/bin:${PATH}
