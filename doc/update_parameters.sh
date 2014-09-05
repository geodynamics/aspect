cd ..
rm -f output/parameters.tex
./aspect doc/manual/empty.prm >/dev/null 2>/dev/null
cp output/parameters.tex doc/manual/ || echo "ERROR: could not copy parameters.tex"
cd doc/manual
echo patching parameters.tex
sed -i 's/LD_LIBRARY_PATH/LD\\_LIBRARY\\_PATH/g' parameters.tex
sed -i 's/tecplot_binary/tecplot\\_binary/g' parameters.tex
sed -i 's/hyper_shell/hyper\\_shell/g' parameters.tex
sed -i 's/\$ASPECT_SOURCE_DIR/\\\$ASPECT\\_SOURCE\\_DIR/g' parameters.tex
sed -i 's/<depth_average.ext>/$<$depth\\_average.ext$>$/g' parameters.tex
sed -i 's/<myplugin.so>/$<$myplugin.so$>$/g' parameters.tex
sed -i 's/<\.\/myplugin.so>/$<$\.\/myplugin.so$>$/g' parameters.tex
sed -i 's/<\[/\[/g' parameters.tex
sed -i 's/\]>/\]/g' parameters.tex
sed -i 's/dynamic_topography.NNNNN/dynamic\\_topography.NNNNN/g' parameters.tex
sed -i 's/Spline_knots.txt/Spline\\_knots.txt/g' parameters.tex
sed -i 's/melt_fraction/melt\\_fraction/g' parameters.tex
sed -i 's/phi\.%d/phi\.\\%d/g' parameters.tex

grep '[^\\]%' parameters.tex

cd ..
echo done