#!/bin/bash
## ---------------------------------------------------------------------
##
##  Copyright (C) 2019 by the ASPECT authors
##
##  This file is part of ASPECT.
##
##  ASPECT is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2, or (at your option)
##  any later version.
##
##  ASPECT is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with ASPECT; see the file LICENSE.  If not see
##  <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

#
# This script runs the clang-tidy tool on the ASPECT code base.
#
# This is based on the deal.II script
# https://github.com/dealii/dealii/blob/master/contrib/utilities/run_clang_tidy.sh
#
# Usage:
# /contrib/utilities/run_clang_tidy.sh SRC_DIR OPTIONAL_CMAKE_ARGS
#   with:
#     SRC_DIR points to an ASPECT source directory
#     OPTIONAL_CMAKE_ARGS are optional arguments to pass to CMake
#  Make sure to run this script in an empty build directory
#
# Requirements:
# Clang 5.0.1+ and have clang, clang++, and run-clang-tidy.py in
# your path.
#
# Run inside docker (assuming a checkout in the current directory):
# docker run --rm -it -v $PWD:/src tjhei/dealii:v9.0.1-full-v9.0.1-r5-clang6 \
#   /src/contrib/utilities/run_clang_tidy.sh /src


# grab first argument and make relative path an absolute one:
SRC=$1
SRC=$(cd "$SRC";pwd)
shift

if test ! -d "$SRC/source" -o ! -d "$SRC/include" -o ! -f "$SRC/CMakeLists.txt" ; then
    echo "Usage:"
    echo "  run_clang_tidy.sh /path/to/ASPECT"
    exit 1
fi
echo "SRC-DIR=$SRC"

# set cmake arguments:
# - export compile commands (so that run-clang-tidy.py works)
# - disable unity/precompiled headers (they confuse clang tidy)
# - append user arguments if present
ARGS=("-D" "CMAKE_EXPORT_COMPILE_COMMANDS=ON" "-D" "ASPECT_UNITY_BUILD=OFF" "-D" "ASPECT_PRECOMPILE_HEADERS=OFF" "$@")

# for a list of checks, see /.clang-tidy
cat "$SRC/.clang-tidy"

if ! [ -x "$(command -v run-clang-tidy.py)" ] || ! [ -x "$(command -v clang++)" ]; then
    echo "make sure clang, clang++, and run-clang-tidy.py (part of clang) are in the path"
    exit 2
fi

CC=clang CXX=clang++ cmake "${ARGS[@]}" "$SRC" || (echo "cmake failed!"; false) || exit 2

# finally run it:
# pipe away stderr (just contains nonsensical "x warnings generated")
# pipe output to output.txt
run-clang-tidy.py -p . -quiet -header-filter="$SRC/include/*" 2>error.txt >output.txt

if grep -E -q '(warning|error): ' output.txt; then
    grep -E '(warning|error): ' output.txt
    exit 4
fi

echo "OK"
exit 0

