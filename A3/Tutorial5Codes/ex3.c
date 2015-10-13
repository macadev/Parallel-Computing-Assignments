
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  int npes, rank;
  int a[10], b[10];
  MPI_Init(&argc, &argv);
  MPI_Status status;


  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
  if (rank == 0) {
    MPI_Ssend(a, 10, MPI_INT, 1, 1, MPI_COMM_WORLD);
    MPI_Ssend(b, 10, MPI_INT, 1, 2, MPI_COMM_WORLD);
  }
  else if (rank == 1) {
    MPI_Recv(b, 10, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
    MPI_Recv(a, 10, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
   }

  MPI_Finalize();
  return 0;
}
