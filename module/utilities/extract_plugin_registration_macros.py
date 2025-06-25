#!/usr/bin/python3

## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

# Given a C++ header file as input, read through it and output an
# equivalent interface module unit that exports a module partition
# with the declarations that are part of the header file, and that
# uses 'import' statements in place of existing '#include' directives.
#
# Call this script via
#   python3 contrib/utilities/convert_header_file_to_interface_module_unit.py in.h out.ccm


import sys
import re

interface_header_files = sys.argv[1:]

assert len(interface_header_files)>1, "No input file name given."

match_define_start = re.compile(r"^# *define ASPECT_REGISTER_")


for interface_header_file in interface_header_files:
    # Read the entire header file into memory:
    input = open(interface_header_file, "r")
    lines = input.readlines()
    input.close()

    for i in range(len(lines)):
        if match_define_start.match(lines[i]):
            # We found a line that starts the definition of a
            # registration macro. Copy it into the output, omitting
            # the endline character.
            print(lines[i][0:-1])

            # Then see whether it ends in a backslash, and if so copy
            # successive lines into the output as well:
            while (lines[i][-2] == '\\'):
                i = i+1;
                print(lines[i][0:-1])
