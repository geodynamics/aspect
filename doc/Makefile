############################################################
# $Id: Makefile 4165 2011-03-28 13:10:33Z bangerth $
############################################################

D           ?= ../../../deal.II
include $D/common/Make.global_options

aspect-files    = $(shell echo ../include/aspect/*.h ../include/aspect/*/*h modules/*h)

aspect.tag: $(aspect-files) options.dox
	@echo "=====aspect=================== Generating documentation"
	@echo '@INCLUDE = options.dox'                                > aspect.dox
	@echo 'INCLUDE_PATH           = $D/include/deal.II'          >> aspect.dox
	@echo 'PROJECT_NAME = "ASPECT mantle convection simulator"'  >> aspect.dox
	@echo 'INPUT = $(aspect-files)'                              >> aspect.dox
	@echo 'HTML_OUTPUT      = doxygen'                           >> aspect.dox
	@echo 'GENERATE_TAGFILE = $@'                                >> aspect.dox
	@doxygen aspect.dox

clean:
	-rm -f aspect.dox aspect.tag
	-rm -rf doxygen
