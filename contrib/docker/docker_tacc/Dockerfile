FROM tacc/tacc-base:ubuntu22.04-impi19.0.9-common

LABEL maintainer <rene.gassmoeller@mailbox.org>

# As of 2024/12/02 the GPG key in the docker image is not public. To avoid the error
# we need to disable the key check.
RUN echo 'deb [trusted=yes] https://apt.repos.intel.com/mpi all main' > /etc/apt/sources.list.d/intel-mpi.list

RUN apt-get update && apt-get -yq upgrade && apt-get -yq install --no-install-recommends \
    libblas-dev \
    liblapack-dev \
    git \
    zlib1g-dev \
    numdiff \
    bc \
    ninja-build \
    && rm -rf /var/lib/apt/lists/*

# Build deal.II dependencies
RUN cd /opt && \
    git clone https://github.com/dealii/candi.git && \
    cd candi && \
    echo "NATIVE_OPTIMIZATIONS=ON" >> local.cfg && \
    echo "DEAL_II_VERSION=v9.6.0" >> local.cfg && \
    echo "BUILD_EXAMPLES=OFF" >> local.cfg && \
    echo "USE_64_BIT_INDICES=ON" >> local.cfg && \
    echo "DEAL_CONFOPTS='-DDEAL_II_WITH_COMPLEX_VALUES=OFF -DDEAL_II_WITH_LAPACK=ON'" >> local.cfg && \
    ./candi.sh -j 4 --packages="once:cmake once:hdf5 once:netcdf once:sundials once:p4est once:trilinos" -p /opt && \
    rm -rf /opt/tmp

# Build deal.II
RUN source /opt/configuration/enable.sh && \
    cd /opt/candi && \
    ./candi.sh -j 4 --packages="dealii" -p /opt && \
    rm -rf /opt/tmp

# Build aspect
ENV Aspect_DIR /opt/aspect
RUN source /opt/configuration/enable.sh && \
    git clone https://github.com/geodynamics/aspect.git $Aspect_DIR && \
    mkdir ${Aspect_DIR}/build && \
    cd ${Aspect_DIR}/build && \
    cmake -DCMAKE_BUILD_TYPE=DebugRelease -DCMAKE_INSTALL_PREFIX=$Aspect_DIR \
    -DASPECT_INSTALL_EXAMPLES=OFF .. && \
    make -j4 && make install && make clean

ENV PATH ${Aspect_DIR}/bin:${PATH}
