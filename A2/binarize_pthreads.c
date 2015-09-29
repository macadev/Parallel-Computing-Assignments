#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#define THRESHOLD 150

unsigned char *image, *new_image;
unsigned width, height;

typedef struct ThreadData {
  int end_cell;
  int start_cell;
} ThreadData;

void *worker_thread(void *arg) {
  ThreadData *td = (ThreadData *) arg;
  unsigned char value;
  int i, j;
  for (i = 0; i < height; i++) {
    for (j = td->start_cell; j < td->end_cell; j++) {

      if (image[4*width*i + 4*j] > THRESHOLD) {
        value = 255;
      } else {
        value = 0;
      }

      new_image[4*width*i + 4*j] = value;
      new_image[4*width*i + 4*j + 1] = value;
      new_image[4*width*i + 4*j + 2] = value;
      new_image[4*width*i + 4*j + 3] = 255;
    }
    
  }
  pthread_exit(NULL);
}

void binarize(char* input_filename, char* output_filename, int thread_count)
{
  unsigned error;

  // load image from PNG into C array
  error = lodepng_decode32_file(&image, &width, &height, input_filename);
  if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
  new_image = malloc(width * height * 4 * sizeof(unsigned char));

  struct timeval start, end; // struct used to compute execution time
  gettimeofday(&start, NULL);  // set starting point

  int region_size = width / thread_count;

  pthread_t threads[thread_count];
  ThreadData *td;
  int i, rc;
  
  for (i = 0 ; i < thread_count; i++) {
    td = malloc(sizeof(ThreadData));
    td->start_cell = (i * region_size);
    if (i == thread_count - 1) {
      td->end_cell = width;
    } else {
      td->end_cell = td->start_cell + region_size;
    }

    rc = pthread_create(&threads[i], NULL, worker_thread, (void *) &(*td));
    if (rc != 0) {
      printf("Error creating threads\n");
      exit(1);
    }
  }

  // Join the threads.
  for (i = 0; i < thread_count; i++) {
    rc = pthread_join(threads[i], NULL);
    if (rc != 0) {
      printf("Error joining threads\n");
      exit(1);
    }
  }

  gettimeofday(&end, NULL);
  printf("\n\nAlgorithm's computational part duration : %ld\n", \
               ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));


  lodepng_encode32_file(output_filename, new_image, width, height);

  free(image);
  free(new_image);
}

int main(int argc, char *argv[])
{
  char* input_filename = argv[1];
  char* output_filename = argv[2];
  int thread_count = atoi(argv[3]);

  binarize(input_filename, output_filename, thread_count);

  return 0;
}
