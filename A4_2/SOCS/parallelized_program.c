#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include "math.h"

#define EPSILON 0.000001


/*
______               _ _      _ _             _  ____________ ___________ 
| ___ \             | | |    | (_)           | | | ___ \ ___ \  ___|  ___|
| |_/ /_ _ _ __ __ _| | | ___| |_ _______  __| | | |_/ / |_/ / |__ | |_   
|  __/ _` | '__/ _` | | |/ _ \ | |_  / _ \/ _` | |    /|    /|  __||  _|  
| | | (_| | | | (_| | | |  __/ | |/ /  __/ (_| | | |\ \| |\ \| |___| |    
\_|  \__,_|_|  \__,_|_|_|\___|_|_/___\___|\__,_| \_| \_\_| \_\____/\_|    
                                                                          
Daniel Macario. 
ID: 260503662
*/

/* 
=================================================
Functions for parallel algorithm implementation
=================================================
*/ 

double * read_user_matrix_from_file_and_store_linearly(char *filename, int *rows, int *columns);
void linear_matrix_write_clicking_probabilities_to_file(double *cp, int rows, int cols);
void parallelized_RREF_helper(double *matrix, int rows, int columns, int rank, int rowsPerProcess, int npes);
void parallelized_RREF(double *linear_A, double* matrix_chunk, int columns, int rowsPerProcess, int rowsForLastProcess, int rank, int npes, int *data_division, int *displacements);
void print_best_acceptance_threshold(double *, int);
void print_linear_matrix(double *matrix, int rows, int columns);
void find_and_divide_by_max_element(double *matrix, int rank, int npes, int rows, int columns);
void print_best_acceptance_threshold_parallel(double *matrix_chunk, int rows, int columns, int rank);


/* 
=================================================
Functions for serial algorithm implementation
=================================================
*/

int equals(double, double);
double ** allocate_matrix(int, int);
double ** read_user_matrix_from_file(char*, int*, int*);
void input_clicking_probabilities(double **, int, int, double *);
void write_clicking_probabilities_to_file(double *, int);
void print_matrix(double **, int, int);
void free_matrix(double **, int);
void RREF(double **, int, int);
void divide_by_max(double **, int, int);

struct timeval start, end;

int main(int argc, char * argv[])
{

  int npes, rank;
  int rows, columns, rowsPerProcess, rowsForLastProcess;
  double **A;
  double *cp;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double *linear_A;
  double *matrix_chunk;
  int data_for_parallelism[4];

  if (rank == 0) {
    /* setup */
    if (argc != 2) printf("please provide a user matrix!");
    
    if (npes == 1) {
      A = read_user_matrix_from_file(argv[1], &rows, &columns);
    } else {
      linear_A = read_user_matrix_from_file_and_store_linearly(argv[1], &rows, &columns);
    }
    cp = malloc(columns * sizeof(double)); // clicking probabilities

    rowsPerProcess = rows / npes;
    rowsForLastProcess = rows - (npes * rowsPerProcess) + rowsPerProcess;

    // Communicate data obtained from reading the file to other processes
    data_for_parallelism[0] = rows;
    data_for_parallelism[1] = columns; 
    data_for_parallelism[2] = rowsPerProcess;
    data_for_parallelism[3] = rowsForLastProcess;
    MPI_Bcast(data_for_parallelism, 4, MPI_INT, 0, MPI_COMM_WORLD);

  } else {

    // Receive data from root process about the file read
    MPI_Bcast(data_for_parallelism, 4, MPI_INT, 0, MPI_COMM_WORLD);
    rows = data_for_parallelism[0];
    columns = data_for_parallelism[1]; 
    rowsPerProcess = data_for_parallelism[2];
    rowsForLastProcess = data_for_parallelism[3];
  }

  /* computation */

  if (npes == 1) {
    /* Serial Execution */
    RREF(A, rows, columns);
    divide_by_max(A, rows, columns);
    input_clicking_probabilities(A, rows, columns, cp);
    print_best_acceptance_threshold(cp, rows);
    write_clicking_probabilities_to_file(cp, rows);
    free_matrix(A, rows);
    free(cp);
  } else {
    /* Parallel Execution */

    int data_division[npes]; // Number of elements each process receives
    int displacements[npes]; // Starting point of the data corresponding to each process
    int i, mynum;

    // Prepare arrays needed to scatter the matrix
    for (i = 0; i < npes - 1; i++) {
      data_division[i] = rowsPerProcess * columns;
      displacements[i] = i * (columns * rowsPerProcess);
    }
    data_division[npes - 1] = rowsForLastProcess * columns;
    displacements[i] = i * (columns * rowsPerProcess);

    if (rank != (npes - 1)) {
      // data for all processes except the last
      matrix_chunk = malloc((rowsPerProcess * columns) * sizeof(double));
    } else {
      // data for last process - it might be larger than the rest
      matrix_chunk = malloc((rowsForLastProcess * columns) * sizeof(double));
    }

    mynum = data_division[rank];
    
    // Scatter the user data among all process in blocks of rows
    MPI_Scatterv(linear_A, data_division, displacements, MPI_DOUBLE, matrix_chunk, mynum, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    parallelized_RREF(linear_A, matrix_chunk, columns, rowsPerProcess, rowsForLastProcess, rank, npes, data_division, displacements);
    if (rank == 0) {
      linear_matrix_write_clicking_probabilities_to_file(linear_A, rows, columns);  
      free(linear_A);
    }
    free(matrix_chunk);
  }

  MPI_Finalize();
  return 0;
}

/*
Writes the last column of the matrix to the clicking_probabilities.txt file
*/
void linear_matrix_write_clicking_probabilities_to_file(double *matrix, int rows, int cols) {
  FILE *output_file;
  int row, last_col_index;
  last_col_index = cols - 1;
  output_file = fopen("clicking_probabilities.txt","w");
  for (row = 0; row < rows; row++) {
    fprintf(output_file, "%lf\n", matrix[row*cols + last_col_index]);
  }

  fclose(output_file);
}

/*
Coordinates the parallelized row reduction of the matrix read from the file. The algorithm
starts by having the first process reduce its matrix chunk, and then communicate the reduced
rows to the other processes so that they may perform a reduction themselves. After, the second
process repeats the exact same procedure, and this goes on until the last chunk of the matrix.
At this point, the matrix is in upper triangle form. The algorithm then reduces the matrix further
with a similar procedure to the one discussed above to obtain RREF form of the matrix.

After, the max element is found in parallel and used to calculate the probabilities. Finally,
the max threshold is calculated in parallel and the results are printed to the probabilites file. 
*/
void parallelized_RREF(double *linear_A, double* matrix_chunk, int columns, int rowsPerProcess, int rowsForLastProcess, 
  int rank, int npes, int *data_division, int *displacements) {

  int i, j, row, col, pivot, correct_rows_value;
  double temp[columns];
  double pivot_val;
  MPI_Status status;

  if (rank == npes - 1) {
    correct_rows_value = rowsForLastProcess;
  } else {
    correct_rows_value = rowsPerProcess;
  }

  for (i = 0; i < rank * rowsPerProcess; i++) {
    MPI_Recv(temp, columns, MPI_DOUBLE, i/rowsPerProcess, 7, MPI_COMM_WORLD, &status);

    for (j = 0; j < columns; j++) {
      if (temp[j] == 1) {
        pivot = j;
        break;
      }
    }
    
    for (row = 0; row < correct_rows_value; row++) {

      pivot_val = matrix_chunk[row*columns + pivot];

      for (col = 0; col < columns; col++) {
        matrix_chunk[row*columns + col] = matrix_chunk[row*columns + col] - pivot_val*temp[col];
      }
    }
  }

  // All rows updated with data from rows above
  // Begin calculation of the RREF for the current matrix chunk
  
  if (rank < npes - 1) {
    parallelized_RREF_helper(matrix_chunk, rowsPerProcess, columns, rank, rowsPerProcess, npes);
  } else {
    parallelized_RREF_helper(matrix_chunk, rowsForLastProcess, columns, rank, rowsPerProcess, npes);
  }

  /* Reduce each individual matrix chunk in parallel */

  int target_row;
  for (row = correct_rows_value - 1; row > 0; row--) {

    pivot = rowsPerProcess * rank + row;

    for (target_row = row - 1; target_row >= 0; target_row--) {
      for (col = 0; col < columns; col++) {
        if (col == pivot) continue;
        matrix_chunk[target_row*columns + col] = matrix_chunk[target_row*columns + col] - matrix_chunk[target_row*columns + pivot]*matrix_chunk[row*columns + col];
      }
      matrix_chunk[target_row*columns + pivot] = 0;
    }
  }  

  /* Reduce upper triangle */

  int ranks_to_receive_from, current_sender_rank, rows_to_receive;
  for (ranks_to_receive_from = npes - rank - 1; ranks_to_receive_from > 0; ranks_to_receive_from--) {
    current_sender_rank = ranks_to_receive_from + rank;

    // receiving extra rows from the last process
    if (current_sender_rank == npes - 1) {
      rows_to_receive = rowsForLastProcess;
    } else {
      rows_to_receive = rowsPerProcess;
    }


    for (i = 0; i < rows_to_receive; i++) {
      MPI_Recv(temp, columns, MPI_DOUBLE, current_sender_rank, 7, MPI_COMM_WORLD, &status);

      pivot = rowsPerProcess * current_sender_rank + rows_to_receive - i - 1;

      // reduce matrix_chunk with the received row
      for (row = 0; row < rowsPerProcess; row++) {
        for (col = 0; col < columns; col++) {
          if (col == pivot) continue;
          matrix_chunk[row*columns + col] = matrix_chunk[row*columns + col] - matrix_chunk[row*columns + pivot]*temp[col];
        }
        matrix_chunk[row*columns + pivot] = 0;
      }
    }
  }

  int destination_rank;
  if (rank != 0) {
    for (row = correct_rows_value - 1; row >= 0; row--) {
      for (j = 0; j < columns; j++) {
        temp[j] = matrix_chunk[row*columns + j];
      }
      for (destination_rank = rank - 1; destination_rank >= 0; destination_rank--) {
        MPI_Send(temp, columns, MPI_DOUBLE, destination_rank, 7, MPI_COMM_WORLD);
      }
    }
  }

  /* Finish Processing */
  
  find_and_divide_by_max_element(matrix_chunk, rank, npes, correct_rows_value, columns);
  print_best_acceptance_threshold_parallel(matrix_chunk, correct_rows_value, columns, rank);

  if (rank == 0) {
    MPI_Gatherv(matrix_chunk, columns*rowsPerProcess, MPI_DOUBLE, linear_A, data_division, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } else if (rank == npes - 1) {
    MPI_Gatherv(matrix_chunk, columns*rowsForLastProcess, MPI_DOUBLE, NULL, data_division, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } else {
    MPI_Gatherv(matrix_chunk, columns*rowsPerProcess, MPI_DOUBLE, NULL, data_division, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

}

/*
Performs the reduction of the matrix argument, and it then broadcasts all the rows from
the matrix so that the other processes may reduce their matrix chunks with them.
*/
void parallelized_RREF_helper(double *matrix, int rows, int columns, int rank, int rowsPerProcess, int npes) {
  int row, i, j;
  double temp[columns];

  int pivot = rowsPerProcess * rank;

  for (row = 0; row < rows; row++) {
    for (j = pivot + row + 1; j < columns; j++) {
      matrix[row*columns + j] = matrix[row*columns + j] / matrix[row*columns + pivot + row];
    }
    
    matrix[row*columns + pivot + row] = 1;
    
    // send your row's calculated data
    for (i = 0; i < columns; i++) {
      temp[i] = matrix[row*columns + i];
    }
    
    for (i = rank + 1; i < npes; i++) {
      MPI_Send(temp, columns, MPI_DOUBLE, i, 7, MPI_COMM_WORLD);
    }

    // update lower rows
    for (i = row + 1; i < rows; i++)
    {
      for (j = pivot + row + 1; j < columns; j++) {
        matrix[i*columns + j] = matrix[i*columns + j] - matrix[i*columns + row + pivot] * temp[j];
      }
      matrix[i*columns + row + pivot] = 0;
    }

  }
}

/*
The max element of the matrix is found in parallel, and used to calculate the probabilities.
*/
void find_and_divide_by_max_element(double *matrix, int rank, int npes, int rows, int columns) {
  MPI_Status status;
  int row, last_col_index, i;
  double max, rec_max, element;
  last_col_index = columns - 1;
  
  max = fabs(matrix[last_col_index]);
  for (row = 1; row < rows; row++) {
    element = fabs(matrix[row*columns + last_col_index]);
    if (element > max) {
      max = element;
    }
  }

  if (rank == 0) {
    for (i = 1; i < npes; i++) {
      MPI_Recv(&rec_max, 1, MPI_DOUBLE, i, 7, MPI_COMM_WORLD, &status);
      if (rec_max > max) {
        max = rec_max;
      }
    }
    MPI_Bcast(&max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } else {
    MPI_Send(&max, 1, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
  }

  // Give every process the max
  MPI_Bcast(&max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (row = 0; row < rows; row++) {
    matrix[row*columns + last_col_index] = fabs(matrix[row*columns + last_col_index] / max);
  }

}

/*
Finds and prints the best acceptance threshold by using MPI_Reduce.
*/
void print_best_acceptance_threshold_parallel(double *matrix_chunk, int rows, int columns, int rank) {
  double thresholds[] = {0.2, 0.4, 0.6, 0.8, 1.0};
  int last_col_index, i, row;
  double best_threshold, local_profit, global_profit, best_profit, probability;
  
  last_col_index = columns - 1;
  for (i = 0; i < 5; i++) {
    local_profit = 0;
    global_profit = 0;
    for (row = 0; row < rows; row++) {
      probability = matrix_chunk[row*columns + last_col_index];
      if (probability > thresholds[i]) {
        local_profit += (1.0 * probability - 2.0 * (1 - probability));
      }
    }

    MPI_Reduce(&local_profit, &global_profit, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
      if (i == 0 || global_profit > best_profit) {
        best_profit = global_profit;
        best_threshold = thresholds[i];
      }
    }
  }

  if (rank == 0) {
    printf("%1.1f\n", best_threshold);
  }
}

/*
The user matrix is read into a one-dimensional array.
*/
double *read_user_matrix_from_file_and_store_linearly(char *filename, int *rows, int *columns) {
  FILE *file;
  file = fopen(filename, "r");

  /* get number of rows and columns*/
  *rows = 1;
  *columns = 1;
  char c;
  int columns_known = 0;
  while(!feof(file)) {
    c = fgetc(file);
    if (c == ' ') {
      if (!columns_known) (*columns)++;
    } 

    if (c == '\n') {
      (*rows)++;
      columns_known = 1;
      continue;
    }
  }

  /* read values into array */
  rewind(file);
  int i;
  double * matrix = (double *) malloc((*rows) * (*columns) * sizeof(double *));
  for (i = 0; i < (*rows)*(*columns); i++) {
    fscanf(file,"%lf",&matrix[i]);
  } 
  fclose(file);

  return matrix;
}

/*
Prints the matrix
*/
void print_linear_matrix(double *matrix, int rows, int columns) {
  int col_counter = 0;
  int i;
  for (i = 0; i < rows*columns; i++) {
    printf("%lf ", matrix[i]);
    col_counter++;
    if ( col_counter == columns) {
      printf("\n");
      col_counter = 0;
    }
  }
}

/*
====================================
== Code Used for Serial Execution ==
====================================
*/

double **read_user_matrix_from_file(char *filename, int *rows, int *columns) {
  FILE *file;
  file = fopen(filename, "r");

  /* get number of rows and columns*/
  *rows = 1;
  *columns = 1;
  char c;
  int columns_known = 0;
  while(!feof(file)) {
    c = fgetc(file);
    if (c == ' ') {
      if (!columns_known) (*columns)++;
    } 

    if (c == '\n') {
      (*rows)++;
      columns_known = 1;
      continue;
    }
  }

  /* read values into array */
  rewind(file);
  int i, j;
  double **matrix = allocate_matrix(*rows, *columns);
  for (i = 0; i < *rows; i++) {
    for (j = 0; j < *columns; j++) {
      fscanf(file,"%lf",&matrix[i][j]);
    }
  } 
  fclose(file);

  return matrix;
}

void RREF(double **matrix, int rows, int columns) {
  /* Gaussian elimination */
  int src_row, dest_row, row, row2, column;
  double pivot;
  for (src_row = 0; src_row < rows; src_row++) {
    for (dest_row = 0; dest_row < rows; dest_row++) {
      if (dest_row == src_row) continue;

      pivot = matrix[dest_row][src_row] / matrix[src_row][src_row];

      // Update all the numbers in a row after reducing!
      for (column = src_row; column < columns; column++) {
        matrix[dest_row][column] = matrix[dest_row][column] - pivot*matrix[src_row][column];
      }
    }
  }

  /* Back-substitution */
  for (row = rows-1; row >= 0; row--) {
    // turn main diagonal elements into 1's, find results
    matrix[row][columns-1] = matrix[row][columns-1] / matrix[row][row];
    matrix[row][row] = 1;
    for (row2 = row-1; row2 >= 0; row2--) {
      matrix[row2][columns-1] += matrix[row2][row]*matrix[row][columns-1];
      matrix[row2][row] = 0;
    }
  }
}

void input_clicking_probabilities(double **matrix, int rows, int columns, double *cp) {
  int row;
  for (row = 0; row < rows; row++) {
    cp[row] = matrix[row][columns-1];
  }
}

void write_clicking_probabilities_to_file(double *cp, int rows) {
  /* write clicking probabilities to file */ 
  FILE *output_file;
  int row;
  output_file = fopen("clicking_probabilities.txt","w");
  for (row = 0; row < rows; row++) {
    fprintf(output_file, "%lf\n", cp[row]);
  }

  fclose(output_file);
}

void free_matrix(double ** matrix, int rows)
{
  int i;
  for (i = 0; i < rows; i++) free(matrix[i]);
  free(matrix);
}

void print_best_acceptance_threshold(double *cp, int rows) {
	double thresholds[] = {0.2, 0.4, 0.6, 0.8, 1.0};
	int i, row;
	double best_threshold, profit, best_profit, probability;
	
	for (i = 0; i < 5; i++) {
	  profit = 0;
	  for (row = 0; row < rows; row++) {
	    probability = cp[row];
	    if (probability > thresholds[i]) {
	      profit += (1.0 * probability - 2.0 * (1 - probability));
	    }
	  }
	  
    if (i == 0 || profit > best_profit) {
      best_profit = profit;
      best_threshold = thresholds[i];
    }
	}
	printf("%1.1f\n", best_threshold);
}

void divide_by_max(double **matrix, int rows, int columns) {
  double max = 0; 
  int row;

  /* get max so we can divide by this later to get probabilities */
  for (row = 0; row < rows; row++) {    
    if (max < fabs(matrix[row][columns-1])) max = fabs(matrix[row][columns-1]);
  }

  /* divide by max and take abs */
  for (row = 0; row < rows; row++) {
    /* check for division by zero */
    if (equals(max,0)) {
      matrix[row][columns-1] = 0;
    } else {
      matrix[row][columns-1] = fabs (matrix[row][columns-1]) / max;
    }
  }
}

int equals(double a, double b) {
  if (fabs(a-b) < EPSILON) return 1;
  else return 0;
}

void print_matrix(double **matrix, int rows, int columns) 
{
  int row, column;
  for (row = 0; row < rows; row++) {
    for (column = 0; column < columns; column++) {
      printf("%lf ",matrix[row][column]);
    }
    printf("\n");
  } 
}

double ** allocate_matrix(int rows, int cols)
{
  int i = 0;
  double ** matrix = (double **) malloc(rows * sizeof(double *));

  for (i = 0; i < rows; i++) {
    matrix[i] = (double *) malloc(cols * sizeof(double));
  }

  return matrix;
}

void print_RREF(double **matrix, int rows, int columns) {
  int row, col;
  for (row = 0; row < rows; row++) {
    for (col = 0; col < columns; col++) {
      printf("%f ", matrix[row][col]);
    }
    printf("\n");
  }
}