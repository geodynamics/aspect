#!/bin/bash

mkdir output

# global ref:
for r in "3" "4" "5" "6" "7" "8"
do
echo "ref $r:"
cp base.prm temp.prm
echo "subsection Mesh refinement" >>temp.prm
echo "set Initial global refinement = $r" >> temp.prm
echo "end" >> temp.prm
../../aspect temp.prm | grep Error
rm -f temp.prm
done



# adaptive ref:
../../aspect adaptive.prm | egrep "freedom|Error"
