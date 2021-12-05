/* circuitSatifiability.c solves the Circuit Satisfiability
 *  Problem using a brute-force sequential solution.
 *
 *   The particular circuit being tested is "wired" into the
 *   logic of function 'checkCircuit'. All combinations of
 *   inputs that satisfy the circuit are printed.
 *
 *   16-bit version by Michael J. Quinn, Sept 2002.
 *   Extended to 32 bits by Joel C. Adams, Sept 2013.
 */
#include <stdio.h>    
#include <limits.h>   
#include <mpi.h>

int checkCircuit (int, long);

int main (int argc, char *argv[]){
  long i;                               // loop variable (32 bits)
  int count = 0;                        // number of solutions

  // MPI Variable 
  int my_rank ,comm_size ,my_count = 0; //my_rank: my current rank comm_size: num of process

  // tree-structured communication
  int   remain_process, receive_count;          // receive_count: receive total numbers of solutions from other process 

  double startTime = 0.0, totalTime = 0.0;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	startTime = MPI_Wtime();

  // Calculate the circuit
  for (i = my_rank; i <= USHRT_MAX; i += comm_size){
    my_count += checkCircuit (my_rank, i);
  }

  /*every step split the processes half. The former process receives my_count from the latter process. 
  Then, devide remain_process in half. when remain_process == 1,leave while loop*/
  
  remain_process = comm_size; //remain_process: records the number of the proceess needs to send or receive the value.

  while(remain_process != 1 && my_rank < remain_process){
    if(remain_process % 2 == 1){
      if(my_rank < remain_process / 2){
        MPI_Recv(&receive_count , 1 , MPI_INT , my_rank + remain_process / 2 + 1 , 1 , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
        my_count += receive_count;
      }
      else if(my_rank > remain_process / 2){
        MPI_Send(&my_count , 1, MPI_INT, my_rank - (remain_process / 2) - 1, 1, MPI_COMM_WORLD);
      }
    }
    else{
      if(my_rank < remain_process / 2){
        MPI_Recv(&receive_count, 1, MPI_INT, my_rank + remain_process / 2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        my_count += receive_count;
      }
      else if(my_rank >= remain_process / 2){                
        MPI_Send(&my_count, 1, MPI_INT, my_rank - (remain_process / 2), 1, MPI_COMM_WORLD);
      }
    }

    remain_process = remain_process / 2 + ((remain_process % 2 == 1) ? 1 : 0);
  }
   
  // Print out the result
  if(my_rank == 0){
    printf ("Process %d finished.\n", my_rank);
    fflush(stdout);
    totalTime = MPI_Wtime() - startTime;
    count = my_count;
	  printf("\nA total of %d solutions were found.\n\n", count);
    fflush(stdout);
    printf("Process %d finished in time %f secs.\n", my_rank, totalTime);
    fflush(stdout);
  }
   
  MPI_Finalize();
  return 0;
}


/* EXTRACT_BIT is a macro that extracts the ith bit of number n.
 *
 * parameters: n, a number;
 *             i, the position of the bit we want to know.
 *
 * return: 1 if 'i'th bit of 'n' is 1; 0 otherwise 
 */

#define EXTRACT_BIT(n,i) ( (n & (1<<i) ) ? 1 : 0)


/* checkCircuit() checks the circuit for a given input.
 * parameters: id, the id of the process checking;
 *             bits, the (long) rep. of the input being checked.
 *
 * output: the binary rep. of bits if the circuit outputs 1
 * return: 1 if the circuit outputs 1; 0 otherwise.
 */

#define SIZE 16

int checkCircuit (int id, long bits) {
   int v[SIZE];        /* Each element is a bit of bits */
   int i;

   for (i = 0; i < SIZE; i++)
     v[i] = EXTRACT_BIT(bits,i);

   if (  (v[0] || v[1]) && (!v[1] || !v[3]) && (v[2] || v[3])
       && (!v[3] || !v[4]) && (v[4] || !v[5])
       && (v[5] || !v[6]) && (v[5] || v[6])
       && (v[6] || !v[15]) && (v[7] || !v[8])
       && (!v[7] || !v[13]) && (v[8] || v[9])
       && (v[8] || !v[9]) && (!v[9] || !v[10])
       && (v[9] || v[11]) && (v[10] || v[11])
       && (v[12] || v[13]) && (v[13] || !v[14])
       && (v[14] || v[15])  )
  {
      printf ("%d) %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d \n", id,
         v[15],v[14],v[13],v[12],
         v[11],v[10],v[9],v[8],v[7],v[6],v[5],v[4],v[3],v[2],v[1],v[0]);
      fflush (stdout);
      return 1;
   } else {
      return 0;
   }
}