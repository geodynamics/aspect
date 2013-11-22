## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the authors of the ASPECT code.
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
##  along with ASPECT; see the file doc/COPYING.  If not see
##  <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

#
# Dashboard configuration:
#

SET(CTEST_PROJECT_NAME "aspect")

SET(CTEST_DROP_METHOD "http")
SET(CTEST_DROP_SITE "cdash.kyomu.43-1.org")
SET(CTEST_DROP_LOCATION "/submit.php?project=aspect")
SET(CTEST_DROP_SITE_CDASH TRUE)

SET(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS   100)
SET(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 300)

# number of lines to submit before an error:
SET(CTEST_CUSTOM_ERROR_PRE_CONTEXT            5)
# number of lines to submit after an error:
SET(CTEST_CUSTOM_ERROR_POST_CONTEXT          20)

#
# Coverage options:
#

SET(CTEST_EXTRA_COVERAGE_GLOB
  # These files should have executable lines and therefore coverage:
  # source/**/*.cc
  )

SET(CTEST_CUSTOM_COVERAGE_EXCLUDE
  "/tests/"
  )
