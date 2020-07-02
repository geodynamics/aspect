FROM tjhei/dealii:v9.2.0-full-v9.2.0-r1-gcc5

LABEL maintainer <rene.gassmoeller@mailbox.org>

USER root

RUN apt-get update && apt-get install -yq --no-install-recommends texlive-generic-extra texlive-base texlive-latex-recommended texlive-latex-base texlive-fonts-recommended texlive-bibtex-extra lmodern texlive-latex-extra texlive-science graphviz python-pip python-setuptools

RUN pip install cpp-coveralls

USER dealii
