#include <stdio.h>  
#include <stdlib.h> 
#include <mpi.h> 

int main(int argc, char *argv[]){
  long long int num_of_toss ,toss_i; // num_of_toss: total number of tosses
  double x, y;                          // x,y value for the random coordinate
  double z, pi_value;                   // z: Used to check if sqrt(x^2+y^2)<=1
  long long int num_in_cir = 0;         // total numbers of tosses in circle
    
  // MPI Variable
  int comm_size ,my_rank;               // comm_size: number of process, my_rank: my current rank
  long long int my_num_in_circle = 0;   // my total number of tosses in circle

  // tree structure communication
  int remain_process;                        // while loop
  long long int recv_num_in_circle;     // receive total numbers of tosses in circle from other process
  double startTime=0.0, totalTime=0.0;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  srand(time(NULL));

  if(my_rank == 0){
    printf("Enter the number of tosses: "); // enter the number of tosses
    fflush(stdout);
    scanf("%lld", &num_of_toss);
    startTime = MPI_Wtime();
  }    

  // Broadcast the number of tosses to every process
  MPI_Bcast(&num_of_toss , 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

  // Calculate the number of tosses in circle
  for(toss_i = 0 + (long long int) my_rank ; toss_i < num_of_toss ; toss_i += (long long int) comm_size){
    x = ((double)rand()) / RAND_MAX;            //gets a random from x coordinate
    y = ((double)rand()) / RAND_MAX;            //gets a random from y coordinate
    z = sqrt(x * x + y * y);    
    if(z <= 1.0){
      my_num_in_circle++;
    }
  }

  /*every step split the processes half. The former process receives my_num_in_circle from the latter process. 
  Then, devide remain_process in half. when remain_process == 1,leave while loop*/
  
  remain_process = comm_size; //remain_process: records the number of the proceess needs to send or receive the value.
    
  while(remain_process != 1 && my_rank < remain_process){
    if(remain_process % 2 == 1){
      if(my_rank < remain_process / 2){
         MPI_Recv(&recv_num_in_circle, 1, MPI_LONG_LONG_INT, my_rank + remain_process / 2 + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         my_num_in_circle += recv_num_in_circle;
      }
      else if(my_rank > remain_process / 2){
        MPI_Send(&my_num_in_circle, 1, MPI_LONG_LONG_INT, my_rank - (remain_process / 2) - 1, 1, MPI_COMM_WORLD);
      }
    }
    else{
      if(my_rank < remain_process / 2){
        MPI_Recv(&recv_num_in_circle, 1, MPI_LONG_LONG_INT, my_rank + remain_process / 2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        my_num_in_circle += recv_num_in_circle;
      }
      else if(my_rank >= remain_process / 2){                
        MPI_Send(&my_num_in_circle, 1, MPI_LONG_LONG_INT, my_rank - (remain_process / 2), 1, MPI_COMM_WORLD);
      }
    }

    remain_process = remain_process / 2 + ((remain_process % 2 == 1) ? 1 : 0);
  }

    if(my_rank == 0){
	    num_in_cir = my_num_in_circle;
	    pi_value = 4.0 * ((double) num_in_cir) / ((double) num_of_toss);
	    printf("\nPI: %f\n", pi_value);
	    fflush(stdout);
	    totalTime = MPI_Wtime() - startTime;
	    printf("Process %d finished in time %f secs.\n", my_rank, totalTime);
	    fflush(stdout);
	  }
	    
	  MPI_Finalize();
	  return 0;
}
