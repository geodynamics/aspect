#include <aspect/simulator.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <iostream>

int f()
{
  using namespace aspect::Utilities;

  const int dim=3;

  // A square polygon
  std::vector< Point<2> > polygon(4);
  polygon[0] = Point<2>(0.0,0.0);
  polygon[1] = Point<2>(1.0,0.0);
  polygon[2] = Point<2>(1.0,1.0);
  polygon[3] = Point<2>(0.0,1.0);

  // A concave polygon
  std::vector< Point<2> > concave_polygon(5);
  concave_polygon[0] = Point<2>(0.0,0.0);
  concave_polygon[1] = Point<2>(1.0,0.0);
  concave_polygon[2] = Point<2>(1.0,1.0);
  concave_polygon[3] = Point<2>(0.5,0.5);
  concave_polygon[4] = Point<2>(0.0,1.0);

  // A selfcrossing polygon
  std::vector< Point<2> > crossing_polygon(4);
  crossing_polygon[0] = Point<2>(0.0,0.0);
  crossing_polygon[1] = Point<2>(1.0,0.0);
  crossing_polygon[2] = Point<2>(1.0,-1.0);
  crossing_polygon[3] = Point<2>(0.0,1.0);


  // Some points inside and outside the polygon
  Point<2> points[] = {Point<2>(0.5,-1), Point<2>(0.5,0.5), Point<2>(0.001,0.2), Point<2>(2.0,2.0), Point<2>(0.25,0.70)};

  std::cout << "Testing distance to polygon function with the following parameters: (polygon 1) "
            << polygon[0] << ", " << polygon[1] << ", " << polygon[2] << ", " << polygon[3] << ", "
            << "(polygon 2) " << concave_polygon[0] << ", " << concave_polygon[1] << ", " << concave_polygon[2] << ", " << concave_polygon[3] << ", " << concave_polygon[4] << ", "
            << "(polygon 3) " << crossing_polygon[0] << ", " << crossing_polygon[1] << ", " << crossing_polygon[2] << ", " << crossing_polygon[3]
            << ", (points) "
            << points[0] << ", " << points[1] << ", " << points[2] << ", " << points[3] << ", " << points[4] << std::endl;

  for (unsigned int i = 0; i < 5; i++)
    {
      std::cout << "Minimal distance of point " << points[i] << " to polygon 1 = " << signed_distance_to_polygon<dim>(polygon,points[i]) << std::endl;
      std::cout << "Minimal distance of point " << points[i] << " to polygon 2 = " << signed_distance_to_polygon<dim>(concave_polygon,points[i]) << std::endl;
      std::cout << "Minimal distance of point " << points[i] << " to polygon 3 = " << signed_distance_to_polygon<dim>(crossing_polygon,points[i]) << std::endl;
    }

  exit(0);
  return 42;
}
// run this function by initializing a global variable by it
int i = f();

