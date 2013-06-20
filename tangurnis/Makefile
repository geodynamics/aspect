all: ba tala tala_c

ba: ba/vel_8.csv ba/vel_16.csv ba/vel_32.csv ba/vel_64.csv ba/vel_128.csv
tala: tala/vel_8.csv tala/vel_16.csv tala/vel_32.csv tala/vel_64.csv tala/vel_128.csv
tala_c: tala_c/vel_8.csv tala_c/vel_16.csv tala_c/vel_32.csv tala_c/vel_64.csv tala_c/vel_128.csv



%.csv: ../lib/aspect ba/tan.prm tala/tan.prm tala_c/tan.prm
	@echo "- generating $@"
	@cp $(*D)/tan.prm $@temp.prm
	@ref=`echo $@ | tr -cd '[[:digit:]]' | awk '{print log($$1)/log(2)}'` && (echo  "subsection Mesh refinement"; echo "set Initial global refinement = $$ref ";echo "end";echo "set Output directory=tangurnis/$(*D)") >> $@temp.prm
	@cd ..; ./lib/aspect tangurnis/$@temp.prm
	@cd $(*D); rm -f parameters* statistics *temp.prm

