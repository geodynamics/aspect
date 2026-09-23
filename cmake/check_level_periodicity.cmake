# Detect the public level-periodicity API, including development snapshots
# that share a version number but predate deal.II PR #20212.
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state(RESET)

set(_periodicity_build RELEASE)
if(ASPECT_BUILD_DEBUG)
  set(_periodicity_build DEBUG)
endif()
set(CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_periodicity_build}}")
set(CMAKE_REQUIRED_INCLUDES ${DEAL_II_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${DEAL_II_TARGET_${_periodicity_build}})

# Recheck if this build directory is reconfigured against another deal.II.
unset(ASPECT_HAVE_LEVEL_PERIODICITY_CONSTRAINTS CACHE)
check_cxx_source_compiles("
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>

void check(const dealii::DoFHandler<2> &dofs,
           const dealii::DoFHandler<2>::level_face_iterator &first,
           const dealii::DoFHandler<2>::level_face_iterator &second)
{
  dealii::MGConstrainedDoFs mg;
  mg.initialize(dofs, dealii::MGLevelObject<dealii::IndexSet>(), false);
  dealii::AffineConstraints<double> constraints;
  dealii::DoFTools::make_periodicity_constraints_on_level(first, second, 0, constraints);
}

int main() { return 0; }
" ASPECT_HAVE_LEVEL_PERIODICITY_CONSTRAINTS)

cmake_pop_check_state()
unset(_periodicity_build)
message(STATUS "Rotated multigrid level periodicity: ${ASPECT_HAVE_LEVEL_PERIODICITY_CONSTRAINTS}")
