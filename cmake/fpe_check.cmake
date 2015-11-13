# Check that we can use feenableexcept. Sets HAVE_FP_EXCEPTIONS.

# The test is a bit more complicated because we want to check that no garbage
# exception is thrown if we convert -std::numeric_limits<double>::max to a
# string. Sadly, this bug only triggers if we link with all the deal.II
# libraries.

INCLUDE (CheckCXXSourceRuns)

SET(_backup_libs ${CMAKE_REQUIRED_LIBRARIES})
SET(_backup_includes ${CMAKE_REQUIRED_LIBRARIES})
SET(_build "RELEASE")
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
IF("${_cmake_build_type}" MATCHES "debug")
  SET(_build "DEBUG")
ENDIF()

LIST(APPEND CMAKE_REQUIRED_LIBRARIES ${DEAL_II_TARGET_${_build}})
LIST(APPEND CMAKE_REQUIRED_INCLUDES ${DEAL_II_INCLUDE_DIRS})

CHECK_CXX_SOURCE_RUNS("
#include <fenv.h>
#include <limits>
#include <sstream>

#include <deal.II/base/utilities.h>

int main()
{
  feenableexcept(FE_DIVBYZERO|FE_INVALID);
  std::ostringstream description;
  const double lower_bound = -std::numeric_limits<double>::max();

  description << lower_bound;
  description << dealii::Utilities::string_to_int (\"1\");

  return 0;
}
" HAVE_FP_EXCEPTIONS)

SET(CMAKE_REQUIRED_LIBRARIES ${_backup_libs})
SET(CMAKE_REQUIRED_INCLUDES ${_backup_includes})
