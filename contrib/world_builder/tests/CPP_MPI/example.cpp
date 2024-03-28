
#include "world_builder/world.h"

#include <stdio.h>

// we don't need the c++ MPI wrappers
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX
#include <mpi.h>

int main(int argc, char *argv[]) {
  // Declare the types which will be needed.
  double temperature = 0.0;
  double x = 120e3;
  double y = 500e3;
  double z = 0;
  double depth = 0;
  double gravity = 10;
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
  
  int MPI_RANK = 0;
  int MPI_SIZE = 1;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
  
  std::unique_ptr<WorldBuilder::World> world = std::unique_ptr<WorldBuilder::World>(new WorldBuilder::World(argv[1], has_output_dir, output_dir, random_number_seed)); 


  if(MPI_RANK == 1)
  {
  printf("MPI size %i and %i according to the world builder, and this is written from mpi rank %i and %i according to the world builder\n",MPI_SIZE,world->MPI_SIZE,MPI_RANK,world->MPI_RANK);

  printf("2d temperature: \n");
  std::array<double,2> coords_2d = {{x, z}};
  temperature = world->temperature(coords_2d,depth);
  printf("temperature in C = %f \n", temperature);

  printf("2d temperature (deprecated): \n");
  coords_2d = {{x, z}};
  temperature = world->temperature(coords_2d,depth,gravity);
  printf("temperature in C = %f \n", temperature);

  printf("3d temperature: \n");
  std::array<double,3> coords_3d = {{x, y, y}};
  temperature = world->temperature(coords_3d,depth);
  printf("temperature in C = %f \n", temperature);

  printf("3d temperature (deprecated): \n");
  coords_3d = {{x, y, y}};
  temperature = world->temperature(coords_3d,depth,gravity);
  printf("temperature in C = %f \n", temperature);

  printf("2d composition: \n");
  composition = world->composition(coords_2d,depth,composition_number);
  printf("composition in C = %f \n", composition);

  printf("3d composition: \n");
  composition = world->composition(coords_3d,depth,composition_number);
  printf("composition in C = %f \n", composition);
  }

  MPI_Finalize();
  return 0;
}
