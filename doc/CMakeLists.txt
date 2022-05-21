# Copyright (C) 2022 by the authors of the ASPECT code.
#
# This file is part of ASPECT.
#
# ASPECT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# ASPECT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ASPECT; see the file doc/COPYING.  If not see
# <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

MESSAGE(STATUS "===== Configuring ASPECT documentation =============")

FIND_PACKAGE(Perl)


####################################################
# Generate input files for the manual
####################################################

# Step 1: Find all of the .prm files that will be used in
#         the manual
FILE(GLOB _dot_prm_directories_cookbooks
     ${CMAKE_SOURCE_DIR}/cookbooks/*/doc)
FILE(GLOB _dot_prm_directories_benchmarks
     ${CMAKE_SOURCE_DIR}/benchmarks/*/doc)
SET(_dot_prm_directories ${CMAKE_SOURCE_DIR}/doc/manual)
LIST(APPEND _dot_prm_directories
    ${_dot_prm_directories_cookbooks}
    ${_dot_prm_directories_benchmarks})

FOREACH(_dir ${_dot_prm_directories})
  FILE(GLOB_RECURSE _my_dot_prm_files
       RELATIVE ${CMAKE_SOURCE_DIR}
       ${_dir}/*prm)

  LIST(APPEND _dot_prm_files
       ${_my_dot_prm_files})
ENDFOREACH()


# Step 2: Set up targets to run each of the .prm files
#         through perl and process them
ADD_CUSTOM_TARGET(build_dot_prm_dot_tex)

FOREACH(_dot_prm ${_dot_prm_files})
  # Set up a command to generate the .prm.tex file
  #
  # We do this using the annotate.pl script, and then processing
  # index entries to contain at most three levels using a
  # second perl script.
  ADD_CUSTOM_COMMAND(
    OUTPUT
      ${CMAKE_CURRENT_BINARY_DIR}/${_dot_prm}.tex
    COMMAND
      ${PERL_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/manual/cookbooks/annotate.pl
          ${CMAKE_SOURCE_DIR}/${_dot_prm}
          > ${CMAKE_CURRENT_BINARY_DIR}/${_dot_prm}.tex
    COMMAND
      ${PERL_EXECUTABLE}
         -pi
          ${CMAKE_CURRENT_SOURCE_DIR}/manual/cookbooks/escape_index_entries.pl
          ${CMAKE_CURRENT_BINARY_DIR}/${_dot_prm}.tex
    DEPENDS
      ${CMAKE_SOURCE_DIR}/${_dot_prm}
      ${CMAKE_CURRENT_SOURCE_DIR}/manual/cookbooks/annotate.pl
      ${CMAKE_CURRENT_SOURCE_DIR}/manual/cookbooks/escape_index_entries.pl
    COMMENT
      "Building marked up version of ${_dot_prm} in latex format."
    )

  # Then we also need a target for the operation above. Target
  # names may not contain slashes, so we need to escape those
  STRING(REGEX REPLACE "[^a-zA-Z0-9]" "" _dot_prm_escaped "${_dot_prm}")
  ADD_CUSTOM_TARGET(build_dot_prm_dot_tex_${_dot_prm_escaped}
    DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/${_dot_prm}.tex
    )

  # Finally add the target to the dependencies of the higher-level
  # target
  ADD_DEPENDENCIES(build_dot_prm_dot_tex build_dot_prm_dot_tex_${_dot_prm_escaped})
ENDFOREACH()


# Step 3: Say that creating the .prm.tex files is a dependency
#         for the manual
ADD_CUSTOM_TARGET(manual)
ADD_DEPENDENCIES(manual build_dot_prm_dot_tex)
