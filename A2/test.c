#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char *argv[]) 
{
  #pragma omp parallel
  {
    printf("%d\n", omp_get_num_threads());
  }
}