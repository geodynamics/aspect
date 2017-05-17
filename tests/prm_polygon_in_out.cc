#include <aspect/simulator.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/geometry_model/initial_topography_model/prm_polygon.h>
#include <aspect/simulator_access.h>

#include <iostream>

int f()
{
  using namespace aspect;
  const int dim=3;
  std::string parameters_clockwise = "100 > 0,0;0,5;5,5;5,0 & -100 > 10,10;10,15;15,15";
  std::string parameters_anticlockwise = "100 > 0,0;5,0;5,5;0,5 & -100 > 10,10;15,15;10,15";
  InitialTopographyModel::PrmPolygon<dim> topo_clockwise;
  InitialTopographyModel::PrmPolygon<dim> topo_anticlockwise;

  ParameterHandler prm_clockwise;
  topo_clockwise.declare_parameters(prm_clockwise);
  prm_clockwise.enter_subsection("Geometry model");
  prm_clockwise.enter_subsection ("Initial topography model");
  prm_clockwise.enter_subsection ("Prm polygon");
  prm_clockwise.set ("Topography parameters", parameters_clockwise);
  prm_clockwise.leave_subsection();
  prm_clockwise.leave_subsection();
  prm_clockwise.leave_subsection();

  ParameterHandler prm_anticlockwise;
  topo_anticlockwise.declare_parameters(prm_anticlockwise);
  prm_anticlockwise.enter_subsection("Geometry model");
  prm_anticlockwise.enter_subsection ("Initial topography model");
  prm_anticlockwise.enter_subsection ("Prm polygon");
  prm_anticlockwise.set ("Topography parameters", parameters_anticlockwise);
  prm_anticlockwise.leave_subsection();
  prm_anticlockwise.leave_subsection();
  prm_anticlockwise.leave_subsection();

  topo_clockwise.parse_parameters(prm_clockwise);
  topo_anticlockwise.parse_parameters(prm_anticlockwise);

  Point<2> points[] = {Point<2>(-1,-1),Point<2>(0,0),Point<2>(0.001,0),Point<2>(0,0.001),Point<2>(0.001,0.001),
                       Point<2>(0,-0.001),Point<2>(-0.001,-0.001),Point<2>(-0.01,2.5),Point<2>(0,2.5),Point<2>(0.01,2.5),
                       Point<2>(0,4.99),Point<2>(0,5),Point<2>(0.01,5),Point<2>(2.5,5),Point<2>(2.5,-0.01),
                       Point<2>(2.5,0),Point<2>(2.5,0.01),Point<2>(4.99,0),Point<2>(5,0),Point<2>(5,0.01),
                       Point<2>(5,2.5),Point<2>(4.99,5),Point<2>(5,4.99),Point<2>(5,5),Point<2>(5,5.01),
                       Point<2>(1,1),Point<2>(12.5,12),Point<2>(11.5,12),Point<2>(15,10),Point<2>(14,10),
                       Point<2>(15,11),Point<2>(5,6),Point<2>(6,5),Point<2>(0,-1),Point<2>(-1,0),
                       Point<2>(-1,5),Point<2>(0,6),Point<2>(5,-1),Point<2>(6,0),Point<2>(6,-1),
                       Point<2>(16,16),Point<2>(12.5,15),Point<2>(10,12.5),Point<2>(10,15),Point<2>(10,10)
                      };

  std::cout << "Testing prm polygon plugin with the following parameters: (clockwise) " << parameters_clockwise << ", (anticlockwise) " << parameters_anticlockwise << std::endl;
  for (unsigned int i = 0; i < 45; i++)
    {
      std::cout << "Clockwise topo at     (" << points[i] << ") = " << topo_clockwise.value(points[i]) << std::endl;
      std::cout << "Anticlockwise Topo at (" << points[i] << ") = " << topo_anticlockwise.value(points[i]) << std::endl;
    }

  exit(0);
  return 42;
}
// run this function by initializing a global variable by it
int i = f();

