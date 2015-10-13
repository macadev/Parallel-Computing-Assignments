#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
  int myrank;
  int npes;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0) {
    int A, D, AmultB, BdivCplusA, secondMult;
    A = 1;

    // broadcast A to all processes
    MPI_Bcast(&A, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // receive result from (A * B)
    MPI_Recv(&AmultB, 1, MPI_INT, 1, 7, MPI_COMM_WORLD, &status);

    // receive result from B/(C + A)
    MPI_Recv(&BdivCplusA, 1, MPI_INT, 2, 7, MPI_COMM_WORLD, &status);
    
    // receive result from (A - 1)*(A - 2)
    MPI_Recv(&secondMult, 1, MPI_INT, 3, 7, MPI_COMM_WORLD, &status);

    // calculate final result
    D = AmultB + BdivCplusA + secondMult;

    printf("Result form expression is: %d\n", D);
  }

  if (myrank == 1) {
    // first multipication
    int B, receivedA, AmultB;
    B = 2;
    
    // receive A
    MPI_Bcast(&receivedA, 1, MPI_INT, 0, MPI_COMM_WORLD);        
    // Send B
    MPI_Send(&B, 1, MPI_INT, 2, 7, MPI_COMM_WORLD);

    AmultB = B * receivedA;
    // send result of A * B
    MPI_Send(&AmultB, 1, MPI_INT, 0, 7, MPI_COMM_WORLD);

  }

  if (myrank == 2) {
    //adition-division
    int C, receivedA, receivedB, CplusA, BdivCplusA, Aminus1, Aminus2, secondMult, valueD, AmultB;
    C = 3;
    
    // recieve A
    MPI_Bcast(&receivedA, 1, MPI_INT, 0, MPI_COMM_WORLD);    
    // receive B
    MPI_Recv(&receivedB, 1, MPI_INT, 1, 7, MPI_COMM_WORLD, &status);

    CplusA = C + receivedA;

    // result from B/(C + A)
    BdivCplusA = receivedB / CplusA;

    MPI_Send(&BdivCplusA, 1, MPI_INT, 0, 7, MPI_COMM_WORLD);
  }

  if (myrank == 3) {
    int receivedA, secondMult;

    // receive A
    MPI_Bcast(&receivedA, 1 , MPI_INT, 0 , MPI_COMM_WORLD);    

    secondMult = (receivedA - 1) * (receivedA - 2);
    // send result from second multiplication
    MPI_Send(&secondMult, 1, MPI_INT, 0, 7, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;
}
