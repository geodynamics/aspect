#!/bin/bash
cd ../..

file=benchmark/solkz/output.txt

rm -f $file


for ref in {1..8}
do

cp benchmark/solkz/sol_kz.prm temp.prm

echo "running ref=$ref..."

echo -e "subsection Mesh refinement\nset Initial global refinement = $ref \nend\n" >> temp.prm


./lib/aspect temp.prm | grep Errors >> $file
#./lib/aspect temp.prm

rm temp.prm
done

