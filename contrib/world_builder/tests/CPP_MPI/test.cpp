#include <stdio.h>

// we don't need the c++ MPI wrappers
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX
#include <mpi.h>

int main(int argc, char *argv[]) {

  // initialize MPI
  int MPI_RANK = 0;
  int MPI_SIZE = 1;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);

  if(MPI_RANK == 1)
    printf(" hello mpi world! There are %i mpi ranks in this test and I am writing from mpi rank %i.\n",MPI_SIZE, MPI_RANK);

  MPI_Finalize();
  return 0;
}
