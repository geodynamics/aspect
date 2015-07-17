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
sed -i 's/box_2d_%s.%d/box\\_2d\\_\\%s.\\%d/g' parameters.tex
sed -i 's/box_2d\.txt/box\\_2d\.txt/g' parameters.tex
sed -i 's/#/\\#/g' parameters.tex

# Process index entries to contain at most three levels (by replacing the
# fourth separator marker ! by /). This is repeated 10 times because only one
# nesting level is removed in each call to sed. The replacement is necessary
# as makeindex only allows for three levels of nesting.
for i in `seq 1 10`; do
  sed -i 's/{\([^!]*\)!\([^!]*\)!\([^!]*\)!\([^}]*\)}/{\1!\2!\3\/\4}/' parameters.tex
done

grep '[^\\]%' parameters.tex

cd ..
echo done
