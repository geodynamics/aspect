#!/bin/bash

# run this script from the doc directory to update the plugin graph for the documentation.
# Note that you need an in-source build or a symbolic link to the ASPECT binary in the
# main directory.

ASPECT=${1:-"./aspect"}

pushd .
cd ..

if test ! -f $ASPECT ; then
  echo "Please provide the absolute path of the ASPECT executable as the first argument to this script or create a link in the main source directory. "
  exit 1
fi

# generate the output file that is used to create the
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
