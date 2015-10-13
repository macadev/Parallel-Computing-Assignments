
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  int npes, rank;
  MPI_Init(&argc, &argv);
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
  int key = 0;
  if (rank == 0) {
    key = 1;
  }
  MPI_Bcast(&key, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank != 0) {
    printf("I am %d and I got the key=%d\n", rank, key);
   }

  MPI_Finalize();
  return 0;
}
