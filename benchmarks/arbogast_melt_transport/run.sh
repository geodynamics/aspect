#!/bin/bash

# global refinement:
filename="convergence_results_may"
echo "writing $filename..."
rm -f $filename
echo "# DoFs L_2_error_velocity L_2_error_fluid_pressure L_2_error_compaction_pressure" >> $filename
for r in "2" "3" "4" "5"
	do
	cp 1d_arbogast_test.prm temp.prm
	echo "refinement $r:"
	echo "subsection Mesh refinement" >>temp.prm
	echo "  set Initial global refinement = $r" >> temp.prm
	echo "end" >> temp.prm
	
	dofs=$((2 ** r * 5))
        # Extract errors and prepend the computed DoFs
        ./aspect-release temp.prm | gawk -v dofs=$dofs '/Errors/ {print dofs, $5, $7, $9}' | sed 's/,//g' >> $filename
	
	grep Solving output/log.txt
	rm -f temp.prm
done

