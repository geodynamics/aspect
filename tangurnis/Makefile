all: ba tala tala_c

ba: ba/vel_8.dat ba/vel_16.dat ba/vel_32.dat ba/vel_64.dat
tala: tala/vel_8.dat tala/vel_16.dat tala/vel_32.dat tala/vel_64.dat
tala_c: tala_c/vel_8.dat tala_c/vel_16.dat tala_c/vel_32.dat tala_c/vel_64.dat


#vel: vel_8.dat vel_16.dat vel_32.dat vel_64.dat vel_128.dat

%.dat: ../lib/aspect tan.prm
	@echo "- generating $@"
	@cp $(*D)/tan.prm $@temp.prm
	@ref=`echo $@ | tr -cd '[[:digit:]]' | awk '{print log($$1)/log(2)}'` && (echo  "subsection Mesh refinement"; echo "set Initial global refinement = $$ref ";echo "end";echo "set Output directory=tangurnis/$(*D)") >> $@temp.prm
	@cd ..; ./lib/aspect tangurnis/$@temp.prm
	@cd $(*D); rm -f parameters* solution* statistics *temp.prm

