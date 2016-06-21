#!/bin/bash

# global refinement:
filename="plane_wave_melt_bands_angle"
rm $filename
echo "#Angle Analytic_growth_rate Numerical_growth_rate" >> $filename

for a in "0" "14.4775121859" "30" "43.4325365578" "61.0449756281" "69.6358651937" "90" "104.4775121859" "120" "133.4325365578" "151.0449756281" "159.6358651937" "180"
do
echo "angle $a:"
cp magmatic_shear_bands.prm temp.prm
echo "subsection Mesh refinement" >>temp.prm
echo "set Initial global refinement = 7" >> temp.prm
echo "end" >> temp.prm
echo "subsection Compositional initial conditions" >>temp.prm
echo "subsection Plane wave melt bands initial condition" >>temp.prm
echo "set Wave number = 4000" >>temp.prm
echo "set Initial band angle = $a" >>temp.prm
echo "end" >>temp.prm
echo "end" >>temp.prm
mpirun -np 4 ./aspect temp.prm  | gawk '/Error/ {print $10, $11, $12}' | sed 's/,//g' >> $filename
rm -f temp.prm
done

