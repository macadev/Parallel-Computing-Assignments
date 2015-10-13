
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
  int sum,rank, mypart;
  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  mypart = f(rank);
  MPI_Reduce(&mypart,&sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
 
  if (rank ==0){
    printf("Sum is :%d\n",sum);}
  
  MPI_Finalize();
  return 0;
}
