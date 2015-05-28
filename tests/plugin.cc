// make sure we can include deal.II and aspect files
#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <typeinfo>

// create a function that is run upon loading the plugin
// and that produces some output
int f()
{
  std::cout << typeid(dealii::Triangulation<2>).name() << " "
            << typeid(aspect::Simulator<2>).name() << std::endl;
  return 42;
}

// run this function by initializing a global variable by it
int i = f();
