#include <stdio.h>
#include <stdlib.h>
#include <semaphore.h>
#include <pthread.h>
#include <unistd.h>
#include <time.h>

static sem_t thread_finish_count_semaphore;
int thread_finished_count = 0;
int num_threads;

// Struct used to associate threads with an ID corresponding
// to the poisition in the sequence of creation.
typedef struct PositionData {
     int position;
} PositionData;

// Sleep the thread for a random number between 0 and 10000 us
void sleepMicroseconds(int us) {
     // usleep sleeps for a specified amount of microseconds.
     if (usleep(us) == -1) {
          printf("Error sleeping thread\n");
          exit(2);
     }
}

// Executing thread sleeps a random number of miliseconds
void *thread_function( void *arg ) {
     PositionData *pos = (PositionData *) arg;
     int position = pos->position;
     sleepMicroseconds(rand() % 10001);
     
     // The semaphore is used for the sole purpose of eliminating the comma
     // at the end of the output string. E.g instead of getting "0,1,2," the program
     // outputs "0,1,2" 
     if (sem_wait(&thread_finish_count_semaphore) == -1) exit(2);
     thread_finished_count += 1;
     printf("%d", position);
     if (thread_finished_count != num_threads) {
         printf(",");
     }
     if (sem_post(&thread_finish_count_semaphore) == -1) exit(2);
}

int main(int argc, char *argv[]) {
     int rc;

     // Function call to be able to generate random numbers
     srand(time(NULL));

     if (argc == 1) {
          printf("Please specify the number of threads as a parameter to the program\n");
          exit(1);
     }

     num_threads = atoi(argv[1]);
     pthread_t threads[num_threads];

     if (sem_init(&thread_finish_count_semaphore, 0, 1) == -1) {
         printf("Error, could initialize semaphore\n");
         exit(1);
     }

     // Span the number of threads specified by the user
     PositionData *posData;
     int i;
     for (i = 0; i < num_threads; i++) {
          posData = malloc(sizeof(PositionData));
          (*posData).position = i;
          rc = pthread_create(&threads[i], NULL, thread_function, (void*) posData);
          if (rc != 0) {
               printf("Error creating threads\n");
               exit(1);
          }
     }

     for (i = 0; i < num_threads; i++) {
          rc = pthread_join(threads[i], NULL);
          if (rc != 0) {
               printf("Error joining threads\n");
               exit(1);
          }
     }

     printf("\n");
     exit(0);
}