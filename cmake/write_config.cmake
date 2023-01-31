# Copyright (C) 2013 - 2021 by the authors of the ASPECT code.
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
# along with ASPECT; see the file LICENSE.  If not see
# <http://www.gnu.org/licenses/>.


SET(_log_detailed "${CMAKE_BINARY_DIR}/detailed.log")
FILE(REMOVE ${_log_detailed})

MACRO(_detailed)
  FILE(APPEND ${_log_detailed} "${ARGN}")
ENDMACRO()

_detailed(
"###
#
#  ASPECT configuration:
#        ASPECT_VERSION:            ${ASPECT_PACKAGE_VERSION}
#        GIT REVISION:              ${ASPECT_GIT_SHORTREV} (${ASPECT_GIT_BRANCH})
#        CMAKE_BUILD_TYPE:          ${CMAKE_BUILD_TYPE}
#
#        DEAL_II_DIR:               ${deal.II_DIR}
#        DEAL_II VERSION:           ${DEAL_II_PACKAGE_VERSION}
#        ASPECT_USE_FP_EXCEPTIONS:  ${ASPECT_USE_FP_EXCEPTIONS}
#        ASPECT_RUN_ALL_TESTS:      ${ASPECT_RUN_ALL_TESTS}
#        ASPECT_USE_SHARED_LIBS:    ${ASPECT_USE_SHARED_LIBS}
#        ASPECT_HAVE_LINK_H:        ${ASPECT_HAVE_LINK_H}
#        ASPECT_WITH_LIBDAP:        ${ASPECT_WITH_LIBDAP}
#        ASPECT_WITH_WORLD_BUILDER: ${ASPECT_WITH_WORLD_BUILDER} ${WORLD_BUILDER_SOURCE_DIR}
#        ASPECT_PRECOMPILE_HEADERS: ${ASPECT_PRECOMPILE_HEADERS}
#        ASPECT_UNITY_BUILD:        ${ASPECT_UNITY_BUILD}
#
#        CMAKE_INSTALL_PREFIX:      ${CMAKE_INSTALL_PREFIX}
#        CMAKE_SOURCE_DIR:          ${CMAKE_SOURCE_DIR}
#        CMAKE_BINARY_DIR:          ${CMAKE_BINARY_DIR}
#        CMAKE_CXX_COMPILER:        ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}
#                                   ${CMAKE_CXX_COMPILER}
#        PARAMETER_GUI_EXECUTABLE:  ${PARAMETER_GUI_EXECUTABLE}
")

IF(CMAKE_C_COMPILER_WORKS)
  _detailed(
"#        CMAKE_C_COMPILER:          ${CMAKE_C_COMPILER}\n")
ENDIF()


IF(DEAL_II_STATIC_EXECUTABLE)
  _detailed(
"#
#        LINKAGE:                   STATIC
")
ELSE()
_detailed(
"#
#        LINKAGE:                   DYNAMIC
")
ENDIF()

GET_PROPERTY(_COMPILE_FLAGS TARGET aspect PROPERTY COMPILE_FLAGS)
_detailed("#
#        COMPILE_FLAGS:             ${_COMPILE_FLAGS}
")

_detailed("#
#        _WITH_CXX14:               ${DEAL_II_WITH_CXX14}
#        _WITH_CXX17:               ${DEAL_II_WITH_CXX17}
#        _MPI_VERSION:              ${DEAL_II_MPI_VERSION}
#        _WITH_64BIT_INDICES:       ${DEAL_II_WITH_64BIT_INDICES}
")


_detailed("#\n###")
