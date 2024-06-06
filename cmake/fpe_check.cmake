# Check that we can use feenableexcept. Sets HAVE_FP_EXCEPTIONS.

# The test is a bit more complicated because we want to check that no garbage
# exception is thrown if we convert -std::numeric_limits<double>::max to a
# string. Sadly, this bug only triggers if we link with all the deal.II
# libraries.

INCLUDE (CheckCXXSourceRuns)

set(_backup_flags ${CMAKE_REQUIRED_FLAGS})
set(_backup_libs ${CMAKE_REQUIRED_LIBRARIES})
set(_backup_includes ${CMAKE_REQUIRED_INCLUDES})

set(_build "RELEASE")
string(TOLOWER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
if("${_cmake_build_type}" MATCHES "debug")
  set(_build "DEBUG")
endif()

list(APPEND CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}")
list(APPEND CMAKE_REQUIRED_LIBRARIES ${DEAL_II_TARGET_${_build}})
list(APPEND CMAKE_REQUIRED_INCLUDES ${DEAL_II_INCLUDE_DIRS})

check_cxx_source_runs("
#include <fenv.h>
#include <limits>
#include <sstream>

#include <deal.II/base/utilities.h>

int main()
{
  // Some implementations seem to not initialize the FPE bits to zero.
  // Make sure we start from a clean state
  feclearexcept(FE_DIVBYZERO|FE_INVALID);

  // Enable floating point exceptions
  feenableexcept(FE_DIVBYZERO|FE_INVALID);

  std::ostringstream description;
  const double lower_bound = -std::numeric_limits<double>::max();

  description << lower_bound;
  description << dealii::Utilities::string_to_int (\"1\");

  return 0;
}
" HAVE_FP_EXCEPTIONS)


set(CMAKE_REQUIRED_FLAGS ${_backup_flags})
set(CMAKE_REQUIRED_LIBRARIES ${_backup_libs})
set(CMAKE_REQUIRED_INCLUDES ${_backup_includes})
