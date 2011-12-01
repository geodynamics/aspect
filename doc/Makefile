############################################################
# $Id: Makefile 4165 2011-03-28 13:10:33Z bangerth $
############################################################

DEAL_DIR?=../../../deal.II
D = $(DEAL_DIR)
include $D/common/Make.global_options

aspect-files    = $(shell echo ../include/aspect/*.h ../include/aspect/*/*h modules/*h)

all: aspect.tag manual.pdf

aspect.tag: $(aspect-files) options.dox
	@if test ! -x "`which doxygen`" ; then \
	  echo "-----------------------------------------" ; \
	  echo "Can't find 'doxygen'." ; \
	  echo "-----------------------------------------" ; \
	  false ; \
	 fi
	@echo "=====aspect=================== Generating documentation"
	@echo '@INCLUDE = options.dox'                                > aspect.dox
	@echo 'INCLUDE_PATH           = $D/include/deal.II'          >> aspect.dox
	@echo 'PROJECT_NAME = "ASPECT mantle convection simulator"'  >> aspect.dox
	@echo 'INPUT = $(aspect-files)'                              >> aspect.dox
	@echo 'HTML_OUTPUT      = doxygen'                           >> aspect.dox
	@echo 'GENERATE_TAGFILE = $@'                                >> aspect.dox
	@doxygen aspect.dox || echo "Running doxygen to generate documentation failed."

manual.pdf: manual/manual.tex manual/manual.bib
	@if test -x "`which pdflatex`" -a -x "`which bibtex`" ; then \
	   (cd manual ; \
		pdflatex manual.tex ; \
		bibtex manual       ; \
		pdflatex manual.tex ; \
		pdflatex manual.tex) ; \
	   mv manual/manual.pdf . ; \
         else \
	  echo "-----------------------------------------" ; \
	  echo "Can't find either 'pdflatex' or 'bibtex'." ; \
	  echo "-----------------------------------------" ; \
	  false ; \
	 fi

clean:
	-rm -f aspect.dox aspect.tag
	-rm -rf doxygen
