#include <mpi.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

/* Constant  */
enum mergeArray {MERGE_SMALL = 0, MERGE_LARGE = 1};
enum num_limit {MAX_NUM = 20};
enum boolean {FALSE = 0, TRUE = 1};
enum process {PROCESS_ROOT = 0};

//function
int arrayCompare(const void *a, const void *b); // compare two value 
void swap(int *a, int *b);                      // swap two value 
void mergeArray(int *local_array, int local_array_size, int *recv_array, int recv_array_size, enum mergeArray state); //merge two array 

int main(void){
  /* Basic Variables  */
  int total_num; //Total number of the array

  /* Global Array  */
  int *global_num_array = NULL;              // Global array 
  
  int *global_num_array_size = NULL;         // Array size of each process in global array  
  int *global_num_array_displacement = NULL; // Displacement of each process in global array 

  /* MPI Variable (basic)  */
  int comm_size, local_rank;  // comm_size :number of process , local_rank: local current rank

  /* MPI Variable (local Array)  */
  int local_remainder_check;  // Checks whether the current process need to add one more value in the array or not
  int *recv_num_array = NULL; /* number array of the recviver */
  int recv_num_array_size;    /* size of the recviver's number array */
  int *local_num_array;       /* local number array */
  int local_num_array_size;   /* size of the local number array */

  /* MPI Variable (time)  */
  double start_time, total_time;
  size_t loop_iter;

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);

  // User input  
  if (local_rank == PROCESS_ROOT){
    char input_str[MAX_NUM];
    printf("Enter how many number you want to calculate:");
    fflush(stdout);
    fgets(input_str, MAX_NUM, stdin);
    total_num = strtol(input_str, NULL, 10);
    start_time = MPI_Wtime();
  }

  srand(local_rank * 100);

  /* Broadcast 'total_num' to each process  */
  MPI_Bcast(&total_num, 1, MPI_INT, PROCESS_ROOT, MPI_COMM_WORLD);

  /* Check the exact size of the array and create */
  local_remainder_check = (((total_num % comm_size) >= (local_rank + 1)) ? 1 : 0);
  local_num_array_size = total_num / comm_size + local_remainder_check;

  local_num_array = (int *)malloc(sizeof(int) * (total_num / comm_size + 1));
  recv_num_array = (int *)malloc(sizeof(int) * (total_num / comm_size + 1));

  /* Randomly create a list of number  */
  for (loop_iter = 0; loop_iter < local_num_array_size; ++loop_iter){
    local_num_array[loop_iter] = rand() % UINT_MAX;
  }

  /* sort local array  */
  qsort(local_num_array, local_num_array_size, sizeof(int), arrayCompare);

  /* Create the array for process root  */
  if (local_rank == PROCESS_ROOT){
    global_num_array = (int *)malloc(sizeof(int) * total_num);
    global_num_array_size = (int *)malloc(sizeof(int) * comm_size);
    global_num_array_displacement = (int *)malloc(sizeof(int) * comm_size);
  }

  /* Gather the array length into 'global_num_array_size'  */
  MPI_Gather(&local_num_array_size, 1, MPI_INT, global_num_array_size, 1, MPI_INT, PROCESS_ROOT, MPI_COMM_WORLD);

  /* Check the displacement of each array  */
  if (local_rank == PROCESS_ROOT){
    global_num_array_displacement[0] = 0;

    for (loop_iter = 1; loop_iter < comm_size; ++loop_iter){
      global_num_array_displacement[loop_iter]
          = global_num_array_displacement[loop_iter - 1] +
            global_num_array_size[loop_iter - 1];
    }
  }

  /* Gather the number array  */
  MPI_Gatherv(local_num_array, local_num_array_size, MPI_INT, global_num_array, global_num_array_size,
              global_num_array_displacement, MPI_INT, PROCESS_ROOT, MPI_COMM_WORLD);

  /* Print out each local array  */
  if (local_rank == PROCESS_ROOT){
    int output_rank = -1;  
    size_t no_counter = 1; // print out position of the value in local array
    printf("\nLocal array:\n");
    for (loop_iter = 0; loop_iter < total_num; ++loop_iter){
      if (loop_iter == global_num_array_displacement[output_rank + 1]){
        no_counter = 1;
        ++output_rank;
        printf("\nprocess %d local array:\n", output_rank);
        fflush(stdout);
      }

      printf("Rank: %3d , No. %5d , %d\n", output_rank, no_counter, global_num_array[loop_iter]);
      fflush(stdout);
      no_counter++;
    }
  }

  /* odd even sort*/
  for (loop_iter = 0; loop_iter < comm_size; ++loop_iter){
    int local_partner; /* The partner of the local process  */

    /* Find the partner of the current phase  */
    if (loop_iter % 2 == 1){ /* Odd phase  */
      local_partner = local_rank + ((local_rank % 2 == 0) ? -1 : 1);
    }
    else{ /* Even phase  */      
      local_partner = local_rank + ((local_rank % 2 == 0) ? 1 : -1);
    }
    /* local process has no partner  */
    if ((local_partner < 0) || (local_partner >= comm_size)){
      continue;
    }

    /* local_rank or local_partner excess total number of array*/
    if ((local_rank >= total_num) || (local_partner >= total_num)){
      continue;
    }

    /* Communicate the partner*/
    if (local_rank % 2 == 0){
      /* check the size of the array  */
      MPI_Send(&local_num_array_size, 1, MPI_INT, local_partner, 0, MPI_COMM_WORLD);
      MPI_Recv(&recv_num_array_size, 1, MPI_INT, local_partner, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      /* If one of the array size is 0, don't need to send and receive*/
      if (local_num_array_size == 0 || recv_num_array_size == 0){
        break;
      }

      /* pass the whole array  */
      MPI_Send(local_num_array, local_num_array_size, MPI_INT, local_partner, 0, MPI_COMM_WORLD);
      MPI_Recv(recv_num_array, recv_num_array_size, MPI_INT, local_partner, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else{
      /*check the size of the array  */
      MPI_Recv(&recv_num_array_size, 1, MPI_INT, local_partner, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&local_num_array_size, 1, MPI_INT, local_partner, 0, MPI_COMM_WORLD);

       /* If one of the array size is 0, don't need to send and receive*/
      if (local_num_array_size == 0 || recv_num_array_size == 0){
        break;
      }

      /*pass the whole array  */
      MPI_Recv(recv_num_array, recv_num_array_size, MPI_INT, local_partner, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(local_num_array, local_num_array_size, MPI_INT, local_partner, 0, MPI_COMM_WORLD);
    }

    /* The smaller rank stored the smaller value , the larger rank stored the larger value  */
    if (local_rank < local_partner){
      mergeArray(local_num_array, local_num_array_size, recv_num_array, recv_num_array_size, MERGE_SMALL);
    }
    else{
      mergeArray(local_num_array, local_num_array_size, recv_num_array, recv_num_array_size, MERGE_LARGE);
    }
  }

  /* finish Odd even sort*/

  /* Gatherv the sorted array  */
  MPI_Gatherv(local_num_array, local_num_array_size, MPI_INT, global_num_array, global_num_array_size,
              global_num_array_displacement, MPI_INT, PROCESS_ROOT, MPI_COMM_WORLD);

  /* Print out the global array  */
  if (local_rank == PROCESS_ROOT){
    printf("\nGlobal Array:\n");

    for (loop_iter = 0; loop_iter < total_num; loop_iter++){
      printf("No. %5d , %d\n", loop_iter + 1, global_num_array[loop_iter]);
      fflush(stdout);
    }

    total_time = MPI_Wtime() - start_time;
    printf("\nThe execute Time: %f\n", total_time);
  }

  free(global_num_array);
  free(global_num_array_size);
  free(local_num_array);
  free(recv_num_array);
  MPI_Finalize();
  return 0;
}

// compare 2 value
int arrayCompare(const void *a, const void *b)
{
  if (*(int *)a <  *(int *)b) {return -1;}
  if (*(int *)a == *(int *)b) {return  0;}
  if (*(int *)a >  *(int *)b) {return  1;}
}

//swap 2 value
void swap(int *a, int *b)
{
  int temp = *a;
  *a = *b;
  *b = temp;
}

//merge 2 array
void mergeArray(int *local_array, int local_array_size,
                int *recv_array, int recv_array_size,
                enum mergeArray state){
  int return_array[local_array_size]; 
  size_t local_iter, recv_iter; 
  int return_iter; 
  if (state == MERGE_SMALL){    
    local_iter = recv_iter = 0;
    return_iter = 0;
    while (return_iter < local_array_size){
      if (local_array[local_iter] <= recv_array[recv_iter]){
        return_array[return_iter] = local_array[local_iter];
        ++return_iter;
        ++local_iter;
      }
      else{
        return_array[return_iter] = recv_array[recv_iter];
        ++return_iter;
        ++recv_iter;
      }
    }
  }
  else{
    local_iter = local_array_size - 1;
    recv_iter = recv_array_size - 1;
    return_iter = local_array_size - 1;
    while (return_iter >= 0){
      if (local_array[local_iter] > recv_array[recv_iter]){
        return_array[return_iter] = local_array[local_iter];
        --return_iter;
        --local_iter;
      }
      else{
        return_array[return_iter] = recv_array[recv_iter];
        --return_iter;
        --recv_iter;
      }
    }
  }
  for (local_iter = 0; local_iter < local_array_size; ++local_iter){
    local_array[local_iter] = return_array[local_iter];
  }
  return;
}