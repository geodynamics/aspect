#!/bin/bash
cd ../..

file=benchmark/solcx/output.txt

rm -f $file


for visc in 1 1e6
do

echo -e "\nvisc = $visc \n\n" >> $file

for ref in {1..8}
do

cp benchmark/solcx/sol_cx.prm temp.prm

echo "running visc=$visc ref=$ref..."

echo -e "subsection Material model\nsubsection SolCx\nset Viscosity jump = $visc\nend\nend\n" >> temp.prm
echo -e "subsection Mesh refinement\nset Initial global refinement = $ref \nend\n" >> temp.prm


./lib/aspect temp.prm | grep Errors >> $file
#./lib/aspect temp.prm

rm temp.prm
done
done
