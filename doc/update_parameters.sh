#!/bin/bash

# run this script from the doc directory to update the parameters for the documentation.
# Note that you need an in-source build or a symbolic link to the ASPECT binary in the
# main directory.

ASPECT=${1:-"./aspect"}

pushd .
cd ..

if test ! -f $ASPECT ; then
  echo "Please provide the absolute path of the ASPECT executable as the first argument to this script or create a link in the main source directory. "
  exit 1
fi


# first thing: run ASPECT so that it produces the parameters.json file that
# documents all parameters and that we can use for the documentation
echo Creating parameters.json
rm -f output/parameters.json
$ASPECT doc/manual/empty.prm >/dev/null 2>/dev/null \
    || { echo "Running ASPECT for parameters.json failed"; exit 1; }

echo Convert parameters to markdown files
./contrib/utilities/jsontomarkdown.py output/parameters.json \
    || { echo "Conversion of parameters to markdown failed"; exit 1; }

# The jsontomarkdown script currently can leave spaces at the ends of lines.  
./contrib/utilities/indent

# next, generate the output file that is used to create the
# connection graph between all plugins and the core of ASPECT
echo Creating plugin graph
$ASPECT --output-plugin-graph doc/manual/empty.prm >plugin_graph.dot 2>/dev/null \
    || { echo "Running ASPECT for the plugin graph failed"; exit 1; }

neato plugin_graph.dot -Tsvg -o plugin_graph.svg \
    || { echo "Can't run neato"; cat plugin_graph.dot; exit 1; }
mv plugin_graph.svg plugin_graph.dot doc/sphinx/user/extending/images/ || echo "ERROR: could not move plugin_graph.*"

popd
echo done
exit 0
