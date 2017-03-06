FROM dealii/dealii:v8.5.pre.4-gcc-mpi-fulldepsmanual-debugrelease

LABEL maintainer <rene.gassmoeller@mailbox.org>

# Build aspect
RUN git clone https://github.com/geodynamics/aspect.git ./aspect && \ 
    mkdir aspect/build && \
    cd aspect/build && \
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DDEAL_II_DIR=$HOME/deal.II-install \
          .. && \
    make -j4 >log 2>&1 && \
    mv aspect $HOME/aspect/aspect-release && \
    cmake -DCMAKE_BUILD_TYPE=Debug \
          -DDEAL_II_DIR=$HOME/deal.II-install \
          .. && \
    make -j4 >log 2>&1 && \
    mv aspect $HOME/aspect/aspect-debug && \
    cd .. && \
    rm -rf build

WORKDIR /home/dealii/aspect
