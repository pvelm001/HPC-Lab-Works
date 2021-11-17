/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   unsigned long int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int    first;        /* Index of first multiple */
   unsigned long int    prime_first;        /* Index of first multiple */
   unsigned long int    global_count = 0; /* Global prime count */
   unsigned long long int    high_value;   /* Highest value on this proc */
   unsigned long int    i;
   int    id;           /* Process ID number */
   unsigned long int    index;        /* Index of current prime */
   unsigned long long int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 3,...,'n' */
   char  *prime_marked;       /* Portion of 3,...,'n' */
   unsigned long long int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long int    proc0_size;   /* Size of proc 0's subarray */
   unsigned long int    prime;        /* Current prime */
   unsigned long int    size;         /* Elements in 'marked' */

   MPI_Init(&argc, &argv); 

   /* Start the timer */

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit(1);
   }

   n = atoll(argv[1]);
   n = (n / 2); //Odd Numbers 

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

   low_value = (id * (n - 1) / p) * 2 + 3;
   high_value = ((id + 1) * (n - 1) / p) * 2 + 1;
   size = (high_value - low_value) / 2 + 1;

   /* Bail out if all the primes used for sieving are
      not all held by process 0 */

   proc0_size = (n - 1) / p;

   if ((proc0_size) * 2 + 3 < (int) sqrt((double) n * 2)) {
      if (!id) printf("Too many processes\n");
      MPI_Finalize();
      exit(1);
   }

   proc0_size = (int) sqrt((double) n * 2);

   /* Allocate this process's share of the array. */
   marked = (char *) malloc(size);
   prime_marked = (char *) malloc(proc0_size);
   if (marked == NULL || prime_marked == NULL) {
      printf("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit(1);
   }

   /* Initialize the marked array */
   for (i = 0; i < size; i++) marked[i] = 0;
   for (i = 0; i < proc0_size; i++) prime_marked[i] = 0;

   index = 0;
   prime = 3;
   do {
      //Finding the first index of the prime multiple
      first = (prime * prime - 3) / 2;
      //Marking the prime multiples
      for (i = first; i < proc0_size; i += prime) prime_marked[i] = 1;
      while (prime_marked[++index]);
      prime = index * 2 + 3;
   } while (prime * prime <= n * 2);
   
   int b = 1<<17; //Local Prime Block Size
   int currb_iter = 0;

   int B = 1<<17; //Global Index Block Size
   int currB_iter = 0;

   //Blocking Local Prime Array
   do {
      //Blocking Global Marked Array
      do {
         index = currb_iter;
         prime = index * 2 + 3;
         do {
            //Finding the first index of the prime multiple
            if (prime * prime > low_value)
               first = (prime * prime - low_value) / 2;
            else {
               if (!(low_value % prime)) first = 0;
               else { 
                  first = prime - (low_value % prime);
                  if (!((first + low_value) % 2)) first = first + prime;
                  first = first / 2; 
               }
            }
            //Marking the prime multiples
            for (i = first + currB_iter; i < MIN(currB_iter + B, size); i += prime) marked[i] = 1;
            while(prime_marked[++index]);
            prime = index * 2 + 3;

         } while (prime * prime <= high_value && index < MIN(currb_iter + b, proc0_size));

         low_value = low_value + B * 2;
         currB_iter = currB_iter + B;

      } while (currB_iter < size);

      low_value = (id * (n - 1) / p) * 2 + 3;
      currB_iter = 0;
      currb_iter = currb_iter + b;
   
   } while (currb_iter < proc0_size);

   /* Count the number of unmarked numbers i.e, Primes in the local processor */
   count = 0;
   for (i = 0; i < size; i++)
      if (!marked[i]) count++;

   /* Reduce the local counts into global count of Processor 0 */
   if (p > 1) MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      
   /* Stop the timer */
   elapsed_time += MPI_Wtime();

   /* Print the results */
   if (!id) printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count + 1, elapsed_time, p);
   MPI_Finalize();
   return 0;
}

