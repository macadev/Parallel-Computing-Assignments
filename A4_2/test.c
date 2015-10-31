#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include "math.h"

int main(int argc, char * argv[])
{

  int npes, rank;
  double *cp;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int i;

  cp = malloc(255 * sizeof(double));
  if (rank == 0) {
    for (i = 0; i < 255; i++) {
      cp[i] = 12312.222;
    }
  }

  MPI_Bcast(cp, 255, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 1){
    for (i = 0; i < 255; i++) {
      printf("%lf ", cp[i]);
    }
    printf("\n");
  }

  MPI_Finalize();
}