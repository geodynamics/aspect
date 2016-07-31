#include <aspect/simulator.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/geometry_model/initial_topography_model/prm_polygon.h>
#include <aspect/simulator_access.h>

#include <iostream>
#include <typeinfo>

int f()
{
  using namespace aspect;
  const int dim=3;
  std::string parameters = "100 > 0,0;5,0;5,5;0,5 & -100 > 10,10;10,15;15,15";
  InitialTopographyModel::PrmPolygon<dim> topo;
  ParameterHandler prm;
  topo.declare_parameters(prm);
  prm.enter_subsection("Geometry model");
  prm.enter_subsection ("Initial topography model");
  prm.enter_subsection ("Prm polygon");
  prm.set ("Topography parameters", parameters);
  prm.leave_subsection();
  prm.leave_subsection();
  prm.leave_subsection();

  topo.parse_parameters(prm);

  std::cout << "Testing prm polygon plugin with the following parameters: " << parameters << std::endl;
  std::cout << "Topo at (-1,-1) = " << topo.value(Point<2>(-1,-1)) << std::endl;
  std::cout << "Topo at (0,0) = " << topo.value(Point<2>(0,0)) << std::endl;
  std::cout << "Topo at (0,5) = " << topo.value(Point<2>(0,5)) << std::endl;
  std::cout << "Topo at (5,0) = " << topo.value(Point<2>(5,0)) << std::endl;
  std::cout << "Topo at (5,5) = " << topo.value(Point<2>(5,5)) << std::endl;
  std::cout << "Topo at (5,5.01) = " << topo.value(Point<2>(5,5.01)) << std::endl;
  std::cout << "Topo at (1,1) = " << topo.value(Point<2>(1,1)) << std::endl;
  std::cout << "Topo at (12.5,12) = " << topo.value(Point<2>(12.5,12)) << std::endl;
  std::cout << "Topo at (11.5,12) = " << topo.value(Point<2>(11.5,12)) << std::endl;

  exit(0);
  return 42;
}
// run this function by initializing a global variable by it
int i = f();

