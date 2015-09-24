#include <stdio.h>
#include <stdlib.h>

int *arrivalQueue;

int main(int argc, char *argv[]) {
     int num_threads;
     int rc;

     if (argc == 1) {
          printf("Please specify the number of threads as a parameter to the program\n");
          exit(1);
     }

     num_threads = atoi(argv[1]);
     printf("Number of threads to span: %d\n", num_threads);
     arrivalQueue = malloc(20);

     arrivalQueue[0] = 12;
     arrivalQueue[19] = 1231;
     arrivalQueue[23] = 123;

     printf("%d\n", arrivalQueue[0]);
     printf("%d\n", arrivalQueue[19]);
     printf("%d\n", arrivalQueue[23]);

     int n = sizeof(arrivalQueue)/sizeof(arrivalQueue[0]);
     printf("sizel %d\n", n);

     exit(0);
}