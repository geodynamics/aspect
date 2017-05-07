FROM dealii/dealii:v8.5.0-gcc-mpi-fulldepscandi-debugrelease

LABEL maintainer <rene.gassmoeller@mailbox.org>

# Build aspect
RUN git clone https://github.com/geodynamics/aspect.git ./aspect && \ 
    mkdir aspect/build-release && \
    cd aspect/build-release && \
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DDEAL_II_DIR=$HOME/deal.II-install \
          .. && \
    make -j4 >log 2>&1 && \
    mv aspect ../aspect-release && \
    make clean && \
    cd .. && \
    mkdir build-debug && \
    cd build-debug && \
    cmake -DCMAKE_BUILD_TYPE=Debug \
          -DDEAL_II_DIR=$HOME/deal.II-install \
          .. && \
    make -j4 >log 2>&1 && \
    mv aspect $HOME/aspect/aspect && \
    make clean

ENV ASPECT_DIR /home/dealii/aspect/build-debug

WORKDIR /home/dealii/aspect
