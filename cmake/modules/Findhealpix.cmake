# The MIT License (MIT)
#
# Copyright (c) 2024 Benoit Bovy
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#
# FindHealpix
# -----------
#
# Find Healpix include directories and libraries.
#
# This module will set the following variables:
#
#  healpix_FOUND          - System has Healpix
#  healpix_INCLUDE_DIRS   - The Healpix include directories
#  healpix_LIBRARIES      - The libraries needed to use Healpix
#

include(FindPackageHandleStandardArgs)

find_path(healpix_INCLUDE_DIR healpix_cxx
  HINTS
    ENV healpix_ROOT
    ENV healpix_DIR
    ENV CONDA_PREFIX
    ${healpix_ROOT_DIR}
  PATH_SUFFIXES
    include
  )

find_library(healpix_LIBRARY
  NAMES healpix_cxx
  HINTS
    ENV healpix_ROOT
    ENV healpix_DIR
    ENV CONDA_PREFIX
    ${healpix_ROOT_DIR}
  PATH_SUFFIXES
    lib
    libs
    Library
  )

find_package_handle_standard_args(healpix
  REQUIRED_VARS healpix_INCLUDE_DIR healpix_LIBRARY
  )

if(healpix_FOUND)
  set(healpix_INCLUDE_DIRS ${healpix_INCLUDE_DIR})
  set(healpix_LIBRARIES ${healpix_LIBRARY})

  add_library(healpix SHARED IMPORTED)
  set_target_properties(healpix PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${healpix_INCLUDE_DIRS}
    IMPORTED_LOCATION ${healpix_LIBRARIES}
    IMPORTED_IMPLIB ${healpix_LIBRARIES}
    )

  mark_as_advanced(healpix_INCLUDE_DIRS healpix_LIBRARIES)
endif()