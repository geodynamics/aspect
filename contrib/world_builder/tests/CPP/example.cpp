
#include "world_builder/wrapper_cpp.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
  // Declare the types which will be needed.
  void * ptr_world =  NULL;
  double temperature = 0.0;
  double x = 120e3;
  double y = 500e3;
  double z = 0;
  double depth = 0;
  double gravity = 10.;
  unsigned int composition_number = 3;
  unsigned int random_number_seed = 1; // use a random number seed larger than zero
  double composition = 0;
  bool has_output_dir = 0; // false
  char output_dir[] = "../../doc/";

   if( argc > 2 ) {
      printf("Too many arguments supplied.\n");
      return 1;
   }
   else if (argc != 2) {
      printf("One argument expected.\n");
      return 1;
   }

  // Show how to call the functions.
  printf("create world \n");

  using namespace wrapper_cpp;
  
  WorldBuilderWrapper world(argv[1],has_output_dir, output_dir, random_number_seed);

  printf("2d temperature: \n");
  printf("temperature in CPP = %f \n", world.temperature_2d(x,z,depth));

  printf("2d temperature (deprecated): \n");
  printf("temperature in CPP = %f \n", world.temperature_2d(x,z,depth, gravity));

  printf("3d temperature: \n");
  printf("temperature in CPP = %f \n", world.temperature_3d(x,y,z,depth));

  printf("3d temperature (deprecated): \n");
  printf("temperature in CPP = %f \n", world.temperature_3d(x,y,z,depth, gravity));

  printf("2d composition: \n");
  printf("composition in CPP = %f \n", world.composition_2d(x,z,depth,composition_number));

  printf("3d composition: \n");
  printf("composition in CPP = %f \n", world.composition_3d(x,y,z,depth,composition_number));

  return 0;
}
