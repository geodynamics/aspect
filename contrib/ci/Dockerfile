FROM tjhei/dealii:v9.2.0-full-v9.2.0-r2-gcc5

LABEL maintainer <rene.gassmoeller@mailbox.org>

USER root

RUN apt-get update && apt-get install -yq --no-install-recommends texlive-generic-extra texlive-base texlive-latex-recommended texlive-latex-base texlive-fonts-recommended texlive-bibtex-extra lmodern texlive-latex-extra texlive-science graphviz python-pip python-setuptools

RUN pip install cpp-coveralls

# we need a newer version of cmake to support unity builds:

RUN cd / && wget https://github.com/Kitware/CMake/releases/download/v3.17.3/cmake-3.17.3-Linux-x86_64.tar.gz && tar xf cmake*.tar.gz && rm cmake*.tar.gz
ENV PATH /cmake-3.17.3-Linux-x86_64/bin:$PATH

USER dealii
