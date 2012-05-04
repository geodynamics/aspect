vel: vel_8.dat vel_16.dat vel_32.dat vel_64.dat vel_128.dat

vel_%.dat: ../lib/aspect tan.prm
	@echo "- generating $@"
	@cp tan.prm temp.prm
	@ref=`echo $@ | tr -cd '[[:digit:]]' | awk '{print log($$1)/log(2)}'` && echo -e "subsection Mesh refinement\nset Initial global refinement = $$ref \nend\n" >> temp.prm
	@cd ..; ./lib/aspect tangurnis/temp.prm
	@rm -f temp.prm
	@mv ../output/output.csv $@

