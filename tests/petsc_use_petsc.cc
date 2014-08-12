// make sure we can include deal.II and aspect files
#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <typeinfo>

// create a function that is run upon loading the plugin
// and that produces some output
int f()
{
#ifdef ASPECT_USE_PETSC
  std::cout << "hey petsc!" << std::endl;
#else
  std::cout << "no petsc :-(" << std::endl;
#endif
  return 42;
}

// run this function by initializing a global variable by it
int i = f();
