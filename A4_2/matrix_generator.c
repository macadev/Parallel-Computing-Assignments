#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

void print_linear_matrix(double *matrix, int rows, int columns) {
  int col_counter = 0;
  int i;
  for (i = 0; i < rows*columns; i++) {
    printf("%lf", matrix[i]);
    if (col_counter != columns - 1) {
      printf(" ");
    }
    col_counter++;
    if ( col_counter == columns && i != rows*columns - 1) {
      printf("\n");
      col_counter = 0;
    }
  }
}

int main(int argc, char * argv[])
{

  double *A, *b, *y, *a, *tmp, *final_y;  // var decls
  int i, j, n, row, r;
  double tstart, tfinish, TotalTime;    // timing decls

  n = 2000;

  A = (double *) malloc(n * (n+1) * sizeof(double *));

  /* Coefficient computer: */
  for( i=0; i<n; i++ ) {    // creates a matrix of random
    for( j=0; j<n; j++ ) {    // scalar equations
      r = rand();
      A[i*n + j] = r;
      A[i*n + n] += j*r;
    }
  }
  print_linear_matrix(A, n, n+1);
}