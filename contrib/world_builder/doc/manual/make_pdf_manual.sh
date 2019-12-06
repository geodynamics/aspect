mkdir -p pdf && \
pdflatex -output-directory pdf -interaction=batchmode manual.tex && \
bibtex pdf/manual && \
pdflatex -output-directory pdf -interaction=batchmode manual.tex && \
pdflatex -output-directory pdf -interaction=batchmode manual.tex
