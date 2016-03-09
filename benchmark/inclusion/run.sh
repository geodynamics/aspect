#!/bin/bash

# global refinement:
echo "--Global Refinement--"
for r in "4" "5" "6" "7" "8" "9"
do
echo "ref $r:"
cp global.prm.base temp.prm
echo "set Output directory = output/global/ref$r" >> temp.prm
echo "subsection Mesh refinement" >> temp.prm
echo "set Initial global refinement = $r" >> temp.prm
echo "end" >> temp.prm
./aspect temp.prm | grep DoFs
rm -f temp.prm
done


# adaptive refinement:
echo "--Adaptive Refinement--"
./aspect adaptive.prm | grep DoFs
