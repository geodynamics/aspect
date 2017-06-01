#!/bin/bash

# run this script from the doc directory to update parameters.tex. Note that
# you need an in-source build or a symbolic link to the ASPECT binary in the
# main directory.

ASPECT=${1:-"./aspect"}

pushd .
cd ..

if test ! -f $ASPECT ; then
  echo "Please provide the path of the ASPECT executable as the first argument to this script or create a link in the main source directory. "
  exit 1
fi


# first thing: run ASPECT so that it produces the parameters.tex file that
# documents all parameters and that we can use for the manual
echo Creating parameters.tex
rm -f output/parameters.tex
$ASPECT doc/manual/empty.prm >/dev/null 2>/dev/null \
    || (echo "Running ASPECT for parameters.tex failed" ; exit 1)
cp output/parameters.tex doc/manual/ \
    || (echo "ERROR: could not copy parameters.tex" ; exit 1)

echo Patching parameters.tex
cd doc/manual
sed -i 's/LD_LIBRARY_PATH/LD\\_LIBRARY\\_PATH/g' parameters.tex
sed -i 's/tecplot_binary/tecplot\\_binary/g' parameters.tex
sed -i 's/MIN_DOUBLE/MIN\\_DOUBLE/g' parameters.tex
sed -i 's/MAX_DOUBLE/MAX\\_DOUBLE/g' parameters.tex
sed -i 's/hyper_shell/hyper\\_shell/g' parameters.tex
sed -i 's/\$ASPECT_SOURCE_DIR/\\\$ASPECT\\_SOURCE\\_DIR/g' parameters.tex
sed -i 's/<depth_average.ext>/$<$depth\\_average.ext$>$/g' parameters.tex
sed -i 's/<myplugin.so>/$<$myplugin.so$>$/g' parameters.tex
sed -i 's/<\.\/myplugin.so>/$<$\.\/myplugin.so$>$/g' parameters.tex
sed -i 's/<\[/\[/g' parameters.tex
sed -i 's/\]>/\]/g' parameters.tex
sed -i 's/dynamic_topography.NNNNN/dynamic\\_topography.NNNNN/g' parameters.tex
sed -i 's/Spline_knots.txt/Spline\\_knots.txt/g' parameters.tex
sed -i 's/stokes_residuals.txt/stokes\\_residuals.txt/g' parameters.tex
sed -i 's/adiabatic_boundary.txt/adiabatic\\_boundary.txt/g' parameters.tex
sed -i 's/melt_fraction/melt\\_fraction/g' parameters.tex
sed -i 's/phi\.%d/phi\.\\%d/g' parameters.tex
sed -i 's/box_2d_%s.%d/box\\_2d\\_\\%s.\\%d/g' parameters.tex
sed -i 's/box_2d_%s./box\\_2d\\_\\%s./g' parameters.tex
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

cd ../..

# next, generate the output file that is used to create the
# connection graph between all plugins and the core of ASPECT
echo Creating plugin graph
$ASPECT --output-plugin-graph doc/manual/empty.prm >plugin_graph.dot 2>/dev/null \
    || (echo "Running ASPECT for the plugin graph failed" ; exit 1)
neato plugin_graph.dot -Tpdf -o plugin_graph.pdf \
    || (echo "Can't run neato" ; exit 1)
cp plugin_graph.* doc/manual || echo "ERROR: could not copy plugin_graph.*"


popd
echo done
