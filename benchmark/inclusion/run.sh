#!/bin/bash
cd ../..

file=benchmark/inclusion/output.txt

rm -f $file


for disc in 0 1
do

disc_str="";
if [[ "$disc" == "1" ]]; then
disc_str="disc";
fi

for visc in 1e3
do

echo -e "$disc_str  visc = $visc\n\n" >> $file

for ref in {1..8}
do

cp benchmark/inclusion/inc.prm temp.prm

echo "running $disc_str visc=$visc ref=$ref..."

echo -e "subsection Material model\nsubsection Inclusion\nset Viscosity jump = $visc\nend\nend\n" >> temp.prm
echo -e "subsection Mesh refinement\nset Initial global refinement = $ref \nend\n" >> temp.prm

if [[ "$disc" == "1" ]]; then
echo -e "subsection Discretization\nset Use locally conservative discretization = true\nend" >> temp.prm
fi

./lib/aspect temp.prm | grep Errors >> $file

rm temp.prm
done
done
done
