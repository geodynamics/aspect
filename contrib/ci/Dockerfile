FROM tjhei/dealii:v9.2.0-full-v9.2.0-r1-gcc5

LABEL maintainer <rene.gassmoeller@mailbox.org>

USER root

RUN apt-get update && apt-get install -yq --no-install-recommends texlive-generic-extra texlive-base texlive-latex-recommended texlive-latex-base texlive-fonts-recommended texlive-bibtex-extra lmodern texlive-latex-extra texlive-science graphviz python-pip python-setuptools

RUN pip install cpp-coveralls

# we need a newer version of cmake to support unity builds:

RUN cd / && wget https://github.com/Kitware/CMake/releases/download/v3.17.3/cmake-3.17.3-Linux-x86_64.tar.gz && tar xf cmake*.tar.gz && rm cmake*.tar.gz
ENV PATH /cmake-3.17.3-Linux-x86_64/bin:$PATH

USER dealii

ENV DEAL_II_VERSION v9.2.0

RUN cd dealii-${DEAL_II_VERSION}-src && mkdir build && cd build && cmake -DDEAL_II_CXX_FLAGS_DEBUG='-O3' -DDEAL_II_WITH_MPI=ON -DDEAL_II_WITH_COMPLEX_VALUES=OFF -DDEAL_II_COMPONENT_EXAMPLES=OFF -DCMAKE_INSTALL_PREFIX=~/dealii-${DEAL_II_VERSION} -DCMAKE_BUILD_TYPE=Debug -GNinja ../ && ninja install && cp summary.log ~/dealii-${DEAL_II_VERSION}/summary.log && cd .. && rm -rf build
