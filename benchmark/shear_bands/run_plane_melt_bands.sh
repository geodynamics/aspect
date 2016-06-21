#!/bin/bash

# global refinement:
filename="plane_wave_melt_bands_8pi"
rm $filename
echo "# DoFs Growth_rate_error" >> $filename

for r in "3" "4" "5" "6" "7" "8"
do
echo "ref $r:"
cp magmatic_shear_bands.prm temp.prm
echo "subsection Mesh refinement" >>temp.prm
echo "set Initial global refinement = $r" >> temp.prm
echo "end" >> temp.prm
mpirun -np 4 ./aspect temp.prm  | gawk '/freedom/ {printf "%s ", $6}; /Error/ {print $13}' | sed 's/,//g' >> $filename
rm -f temp.prm
done

filename="plane_wave_melt_bands_16pi"
rm $filename
echo "# DoFs Growth_rate_error" >> $filename

for r in "3" "4" "5" "6" "7" "8"
do
echo "ref $r:"
cp magmatic_shear_bands.prm temp.prm
echo "subsection Mesh refinement" >>temp.prm
echo "set Initial global refinement = $r" >> temp.prm
echo "end" >> temp.prm
echo "subsection Compositional initial conditions" >>temp.prm
echo "subsection Plane wave melt bands initial condition" >>temp.prm
echo "set Wave number = 8000" >>temp.prm
echo "end" >>temp.prm
echo "end" >>temp.prm

mpirun -np 4 ./aspect temp.prm  | gawk '/freedom/ {printf "%s ", $6}; /Error/ {print $13}' | sed 's/,//g' >> $filename
rm -f temp.prm
done
