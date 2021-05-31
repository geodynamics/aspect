FROM gassmoeller/aspect-tester:focal-9.3

LABEL maintainer <rene.gassmoeller@mailbox.org>

RUN apt-get update && apt-get install -yq --no-install-recommends texlive-plain-generic texlive-base texlive-latex-recommended texlive-latex-base texlive-fonts-recommended texlive-bibtex-extra lmodern texlive-latex-extra texlive-science graphviz python3-pip python-setuptools

RUN pip3 install cpp-coveralls
