#
# Copyright (C) 2013 by Matthias Maier
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
#
# $Id$



CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
MESSAGE("-- This is CTest ${CMAKE_VERSION}")


#
# And finally submit:
#

MESSAGE("-- Running CTEST_SUBMIT()")
CTEST_SUBMIT(RETURN_VALUE _res)

IF("${_res}" STREQUAL "0")
  MESSAGE("-- Submission successful. Goodbye!")
ENDIF()

