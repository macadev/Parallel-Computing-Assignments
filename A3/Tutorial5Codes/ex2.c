
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int f(int x)
{
  x = x+1;
  return x;
}

int main(int argc, char *argv[])
{
  int v0,v1,sum,rank;
  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (rank == 1) {
    v1 = f(1);
    sleep(5);
    MPI_Send(&v1,1,MPI_INT,0,50,MPI_COMM_WORLD);
  }
  else if (rank == 0) {
    v0=f(0);
    MPI_Recv(&v1,1,MPI_INT, 1,50,MPI_COMM_WORLD,&status);
    sum=v0+v1;
    printf("f(0)+f(1) = %d + %d = %d\n",v0,v1,sum);
  }
  MPI_Finalize();
  return 0;
}
