/* Parallel sample sort
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>


static int compare(const void *a, const void *b)
{
  int *da = (int *)a;
  int *db = (int *)b;

  if (*da > *db)
    return 1;
  else if (*da < *db)
    return -1;
  else
    return 0;
}

int main( int argc, char *argv[])
{
  int rank, size, count;
  int i, j, N;
  int *vec, *new_vec, *sample, *all_sample, *splitter;
  int *scounts, *sdispls, *rcounts, *rdispls;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Number of random numbers per processor (this should be increased
   * for actual tests or could be passed in through the command line */
  if (argc != 2){
    N = 100;
  }
  else{
    sscanf(argv[1], "%d", &N);
  }

  vec = calloc(N, sizeof(int));
  /* seed random number generator differently on every core */
  srand((unsigned int) (rank + 393919));

  /* fill vector with random integers */
  for (i = 0; i < N; ++i) {
    vec[i] = rand();
  }
  printf("rank: %d, first entry: %d\n", rank, vec[0]);

  /* sort locally */
  qsort(vec, N, sizeof(int), compare);

  /* randomly sample s entries from vector or select local splitters,
   * i.e., every N/P-th entry of the sorted vector */
  sample = calloc(N/size, sizeof(int));
  for (i = 0; i < N/size; i++){
    sample[i] = vec[1 + size*i];
  }

  /* every processor communicates the selected entries
   * to the root processor; use for instance an MPI_Gather */
  if (rank == 0) {
    all_sample = calloc(N, sizeof(int));
  }
  MPI_Gather(sample, N/size, MPI_INT, all_sample, N/size, MPI_INT, 0,
           MPI_COMM_WORLD);

  /* root processor does a sort, determinates splitters that
   * split the data into P buckets of approximately the same size */
  splitter = calloc(size-1, sizeof(int));
  if (rank == 0) {
    qsort(all_sample, N, sizeof(int), compare);
    for (i = 0; i < size-1; i++){
      splitter[i] = all_sample[N/size*(i+1)];
    }
  }

  /* root process broadcasts splitters */
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(splitter, size-1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  /* every processor uses the obtained splitters to decide
   * which integers need to be sent to which other processor (local bins) */
  count = 1;
  j = 0;
  scounts = calloc(size, sizeof(int));
  sdispls = calloc(size, sizeof(int));
  rcounts = calloc(size, sizeof(int));
  rdispls = calloc(size, sizeof(int));
  for (i = 0; i < N; i++){
    if (vec[i] > splitter[j] && j < size){
      scounts[j] = count;
      count = 1;
      j += 1;
    }
    else
      count += 1;
   }
  scounts[size-1] = count;

  sdispls[0] = 0;
  for (i = 1; i < size; i++){
    sdispls[i] = 0;
    for (j = 0; j < i; j++){
      sdispls[i] = sdispls[i] + scounts[j];
    }
  }

  /* send and receive: either you use MPI_AlltoallV, or
   * (and that might be easier), use an MPI_Alltoall to share
   * with every processor how many integers it should expect,
   * and then use MPI_Send and MPI_Recv to exchange the data */
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, MPI_COMM_WORLD);
  int sum = 0;
  rdispls[0] = 0;
  for (i = 0; i < size; i++){
    sum = sum + rcounts[i];
    rdispls[i] = sum - rcounts[i];
  }
  new_vec = calloc(N*N, sizeof(int));
  MPI_Alltoallv(vec, scounts, sdispls , MPI_INT, new_vec,  rcounts, rdispls, MPI_INT, MPI_COMM_WORLD);
  
  /* do a local sort */
  qsort(new_vec, sum, sizeof(int), compare);
  
  /* every processor writes its result to a file */
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "output%02d.txt", rank);
    fd = fopen(filename,"w+");

    if(NULL == fd)
    {
      printf("Error opening file \n");
      return 1;
    }

    fprintf(fd, "rank %d contains:\n", rank);
    for(i = 0; i < sum; ++i)
      fprintf(fd, "%d\n", new_vec[i]);

    fclose(fd);
    
  /* clean up*/
    if (rank == 0){
      free(all_sample);
    }
  free(new_vec);
  free(scounts);
  free(rcounts);
  free(sdispls);
  free(rdispls);
  free(splitter);
  free(sample);
  free(vec);
  MPI_Finalize();
  return 0;
}
