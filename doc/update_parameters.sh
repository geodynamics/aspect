cd ..
make
./aspect doc/manual/empty.prm >/dev/null 2>/dev/null
cp output/parameters.tex doc/manual/
cd doc/manual
echo patching parameters.tex
sed -i 's/LD_LIBRARY_PATH/LD\\_LIBRARY\\_PATH/g' parameters.tex
sed -i 's/tecplot_binary/tecplot\\_binary/g' parameters.tex
cd ..
echo done