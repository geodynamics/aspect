.PHONY: all

nproc=1

refinements:=4 5 6 7

# bases_without is "case1a case1b ..."
bases:=$(wildcard base_case*.prm)
bases_without:=$(subst base_,,$(basename $(bases)))

# create a list of "caseXX_refY.stat" names:
files:=$(subst base_,,$(bases))
prmfiles:=$(foreach ref,$(refinements),$(subst .prm,_ref$(ref).prm,$(files)))
statfiles:=$(addsuffix .stat, $(basename $(prmfiles)))

all: $(statfiles)
	./get_statistics *.stat


$(foreach ref,$(refinements),%_ref$(ref).prm): base_%.prm
	@for ref in $(refinements) ; do \
		f="$*_ref$$ref.prm";\
		echo "making $$f";\
		echo "include base_$*.prm" > $$f; \
		echo "set Additional shared libraries= ./plugin/libblankenbach.so" >>$$f; \
		echo "  subsection Mesh refinement" >> $$f; \
		echo "set Initial global refinement=$$ref" >> $$f; \
		echo "end" >> $$f; \
		echo "set Output directory=output-$*_ref$${ref}" >> $$f; \
	done 


%.stat: %.prm
	@echo "** $*.prm -> output-$*/statistics -> $@"
	mpirun -n $(nproc) ./aspect $*.prm >/dev/null
	@cp output-$*/statistics $@
	#@rm -rf output-$*
	@echo created $@.
