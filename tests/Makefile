DEAL_DIR?=../../../deal.II
D = $(DEAL_DIR)
include $D/common/Make.global_options

tests = $(shell echo *.prm | sed s/.prm//g)

SUBDIR = x-app-aspect
ID     = `id -un`"@"`hostname`


all: $(tests) Makefile.dep


# a test consists of creating the output (through the dependency
# listed in Makefile.dep) and then comparing the output
#
# the task is made a bit more complicated because screen outputs
# contain run times which we need to screen. fortunately, all
# run times are reported on lines starting with a '|', so
# simply grep them out before running the diff command
$(tests):
	@if test -d $@ ; then \
	   for i in `cd $@; echo *` ; do \
             cat $@/$$i | egrep -v '^\|' > $@/$$i.notime ; \
             cat output-$@/$$i | egrep -v '^\|' > output-$@/$$i.notime ; \
	     if ((diff -q $@/$$i.notime output-$@/$$i.notime 2>&1) > /dev/null) ; then \
		  : ; \
	     else \
		  echo "Files $@/$$i output-$@/$$i differ" ; \
		  failed=true ; \
	     fi ; \
	     rm -f $@/$$i.notime output-$@/$$i.notime ; \
	   done ; \
	 else \
	   failed=true ; \
	 fi ; \
	 if test "x$$failed" = "x" ; then \
	  echo "=====$(SUBDIR)=======================testing============ $< OK" ; \
	  exit 0 ; \
	 else \
	  echo "=====$(SUBDIR)=======================testing============ $< FAILED" ; \
	  exit 1 ; \
	 fi

# creating the output involves running ../lib/aspect; at first, however,
# we make sure that we don't care what the parameter file says
# regarding the output directory, but that we simply set what is
# appropriate. we do so in a way that creates a temporary parameter
# file with a fixed output directory and if the two are different
# overwrite the original one
output-%:
	@echo =====$(SUBDIR)=======================testing============ $<
	@rm -rf $@
	@mkdir $@
	@cat $< | $(PERL) -pi -e 's/(Output directory\s*=).*/\1 $@/g;' \
		> $(<).tmp
	@if (diff -q $< $(<).tmp > /dev/null) ; then \
		rm $(<).tmp ; \
	 else \
		mv $(<).tmp $< ; \
	 fi
	@MPI=`grep "MPI:" $< | perl -pi -e 's/.*: *(\d+).*/\1/g;'` ; \
	    if test -z "$$MPI" ; then \
		../lib/aspect $< > $@/screen-output ; \
	    else \
		mpirun -np $$MPI ../lib/aspect $< > $@/screen-output ; \
	    fi
	@$(PERL) -pi ./normalize.pl $@/*
	@rm -f $@/log


# Creating the report is done by first running all tests (allowing
# running them in parallel) using -k so that we really do get to all
# test, and then re-doing the check
report:
	@-$(MAKE) -k $(tests)
	@-svn info $D | grep '^Revision'
	@echo 'Date: ' `date -u +"%Y %j  %F  %U-%w"`
	@echo 'Id: ' $(ID)-aspect@`svn info . | grep '^Revision' | perl -pi -e 's/Revision: //g;'`
	@for i in $(tests) ; do \
	    if (($(MAKE) $$i 2>&1) > /dev/null) ; then \
	      echo `date -u +"%Y-%m-%d %H:%M"` "  +  $(SUBDIR)/$$i" ; \
	    else \
	      echo `date -u +"%Y-%m-%d %H:%M"` " 0   $(SUBDIR)/$$i" ; \
	    fi ; \
	 done


report+mail:
	@$(MAKE) report | tee testresults
	@cat testresults | mail regression-tests@dealii.org
	@rm testresults


clean:
	rm -rf output-*
	rm -rf Makefile.dep
	rm -rf *~


Makefile.dep: Makefile $(shell echo *.prm)
	@echo =====$(SUBDIR)========================================== $@
	@for i in $(tests) ; do \
		echo "output-$$i: $$i.prm ../lib/aspect" ; \
		echo "$$i: output-$$i" ; \
	done > $@


include Makefile.dep

.PHONY: all report report+mail $(tests) clean
