FROM tjhei/dealii:v9.0.1-full-v9.0.1-r5-gcc5

LABEL maintainer <rene.gassmoeller@mailbox.org>

# Build aspect, replace git checkout command to create image for release
RUN git clone https://github.com/geodynamics/aspect.git ./aspect && \ 
    mkdir aspect/build-release && \
    cd aspect/build-release && \
    git checkout master && \
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DASPECT_PRECOMPILE_HEADERS=ON \
          .. && \
    make -j6 && \
    mv aspect ../aspect-release && \
    make clean && \
    cd .. && \
    mkdir build-debug && \
    cd build-debug && \
    cmake -DCMAKE_BUILD_TYPE=Debug \
          -DASPECT_PRECOMPILE_HEADERS=ON \
          .. && \
    make -j6 && \
    mv aspect $HOME/aspect/aspect && \
    make clean

ENV ASPECT_DIR /home/dealii/aspect/build-debug

WORKDIR /home/dealii/aspect
