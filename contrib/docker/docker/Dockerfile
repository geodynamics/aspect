FROM tjhei/dealii:v9.2.0-full-v9.2.0-r2-gcc5

LABEL maintainer <rene.gassmoeller@mailbox.org>

# we need a newer version of cmake to support unity builds:
RUN cd $HOME/libs && wget https://github.com/Kitware/CMake/releases/download/v3.17.3/cmake-3.17.3-Linux-x86_64.tar.gz && tar xf cmake*.tar.gz && rm cmake*.tar.gz
ENV PATH $HOME/libs/cmake-3.17.3-Linux-x86_64/bin:$PATH

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
