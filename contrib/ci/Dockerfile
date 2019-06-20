FROM tjhei/dealii:v9.0.1-full-v9.0.1-r5-gcc5

LABEL maintainer <rene.gassmoeller@mailbox.org>

USER root

RUN apt-get update && apt-get install -yq --no-install-recommends texlive-generic-extra texlive-base texlive-latex-recommended texlive-latex-base texlive-fonts-recommended texlive-bibtex-extra lmodern texlive-latex-extra texlive-science

USER dealii
