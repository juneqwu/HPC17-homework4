/******************************************************************************
* FILE: mpi_bug5.c
* DESCRIPTION: 
*   This is an "unsafe" program. It's behavior varies depending upon the
*   platform and MPI library
* AUTHOR: Blaise Barney 
* LAST REVISED: 01/24/09
* Hint: If possible, try to run the program on two different machines,
* which are connected through a network. You should see uneven timings;
* try to understand/explain them.
******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define MSGSIZE 2000

int main (int argc, char *argv[])
{
int        numtasks, rank, i, tag=111, dest=1, source=0, count=0;
char       data[MSGSIZE];
double     start, end, result;
MPI_Status status;
MPI_Request req;

MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

if (rank == 0) {
  printf ("mpi_bug5 has started...\n");
  if (numtasks > 2) 
    printf("INFO: Number of tasks= %d. Only using 2 tasks.\n", numtasks);
  }

/******************************* Send task **********************************/
if (rank == 0) {

  /* Initialize send data */
  for(i=0; i<MSGSIZE; i++)
     data[i] =  'x';

  start = MPI_Wtime();
  
  /*stop the while loop at some point) */
  while (count < 200) {
    MPI_Isend(data, MSGSIZE, MPI_BYTE, dest, tag, MPI_COMM_WORLD, &req);
    count++;
    if (count % 10 == 0) {
      end = MPI_Wtime();
      printf("Count= %d  Time= %f sec.\n", count, end-start);
      start = MPI_Wtime();
      }
    }
  }


/****************************** Receive task ********************************/

if (rank == 1) {
  count = 0;
  while (count < 200) {
    MPI_Irecv(data, MSGSIZE, MPI_BYTE, source, tag, MPI_COMM_WORLD, &req);
    /* Do some work  - at least more than the send task */
    result = 0.0;
    for (i=0; i < 1000000; i++) 
      result = result + (double)random();
    count ++;
   }
  }

MPI_Finalize();
 
return 0;
}

/* The program dies on my computer some times, but runs with uneven timing
 on crunchy1 and crunchy3 using -perhost 1 commend, e.g. 
 Count= 530  Time= 0.000029 sec.
 Count= 540  Time= 0.155222 sec.
 The reason is that the receiving side is doing much more work than the sending end,
 but the sending side is still keeping sending messages to the other side regardless 
 of the status. The receving side of the buffer got flooded if it doesn't finish its
 work fast enough to receive.  Fix: switch to non-blocking Isend/Irecv. Another fix 
 is to put a MPI_barrier at the end of each Send/Recv.
 */
