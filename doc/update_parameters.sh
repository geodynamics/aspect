#!/bin/bash

# update the parameters for the documentation.
# This script can be run from any working directory; relative paths to the
# ASPECT executable are interpreted relative to the caller's working directory.
# Note that you need an in-source build or a symbolic link to the ASPECT binary
# in the main directory.

ASPECT_INPUT=${1:-"./aspect"}
CALLER_DIR="$(pwd)"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

case "${ASPECT_INPUT}" in
  /*)
    ASPECT="${ASPECT_INPUT}"
    ;;
  *)
    ASPECT="${CALLER_DIR}/${ASPECT_INPUT}"
    ;;
esac

pushd . >/dev/null
cd "${REPO_ROOT}" || exit 1

if test ! -f $ASPECT ; then
  echo "Please provide the path to the ASPECT executable as the first argument to this script, or create a link in the main source directory."
  exit 1
fi


# run ASPECT so that it produces the parameters.json file that
# documents all parameters and that we can use for the documentation
echo Creating parameters.json
rm -f output/parameters.json
$ASPECT doc/manual/empty.prm >/dev/null 2>/dev/null

if test ! -f output/parameters.json ; then
  echo "Running ASPECT for parameters.json failed"
  exit 1
fi

echo Convert parameters to markdown files
./contrib/utilities/jsontomarkdown.py output/parameters.json \
    || { echo "Conversion of parameters to markdown failed"; exit 1; }

# The jsontomarkdown script currently can leave spaces at the ends of lines.  
./contrib/utilities/indent

popd >/dev/null
echo done
exit 0
