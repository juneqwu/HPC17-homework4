/* MPI-parallel 2D Jacobi smoothing to solve -u''=f
 * Global vector has N unknowns, each processor works with its
 * part, which has lN = N/p unknowns.
 * Author: June Wu
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "util.h"
#include <string.h>

/* compuate global residual, assuming ghost values are updated */
double compute_residual(double *lu, int lN, double invhsq)
{
  int i,j;
  double tmp, gres = 0.0, lres = 0.0;

  for (i = 1; i <= lN; i++){
    for (j = 1; j <= lN; j++){
      tmp = ((4.0* lu[j+(lN+2)*i] - lu[j+(lN+2)*i+1] - lu[j+(lN+2)*i-1] - lu[j+(lN+2)*(i-1)] - lu[j+(lN+2)*(i+1)]) * invhsq - 1);
      lres += tmp * tmp;
    }
  }
  /* use allreduce for convenience; a reduce would also be sufficient */
  MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sqrt(gres);
}


int main(int argc, char * argv[])
{
  int mpirank, i, j, p, N, lN, iter, max_iters;
  MPI_Status status;
  MPI_Request request_out1, request_in1;
  MPI_Request request_out2, request_in2;
  MPI_Request request_out3, request_in3;
  MPI_Request request_out4, request_in4;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  /* get name of host running MPI process */
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);

  sscanf(argv[1], "%d", &N);
  sscanf(argv[2], "%d", &max_iters);

  /* compute number of unknowns handled by each process */
  lN = N / p;
  if ((N % p != 0) && mpirank == 0 ) {
    printf("N: %d, local N: %d\n", N, lN);
    printf("Exiting. N must be a multiple of p\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  /* timing */
  MPI_Barrier(MPI_COMM_WORLD);
  timestamp_type time1, time2;
  get_timestamp(&time1);

  /* Allocation of vectors, including left and right ghost points */
  double * lu    = (double *) calloc(sizeof(double), pow((lN + 2), 2));
  double * lunew = (double *) calloc(sizeof(double), pow((lN + 2), 2));
  double * lutemp;

  double h = 1.0 / (N + 1);
  double hsq = h * h;
  double invhsq = 1./hsq;
  double gres, gres0, tol = 1e-5;

  /* initial residual */
  gres0 = compute_residual(lu, lN, invhsq);
  gres = gres0;

  for (iter = 0; iter < max_iters && gres/gres0 > tol; iter++) {
    /* interleaf computation and communication: compute the first
     * and last value, which are communicated with non-blocking
     * send/recv. During that communication, do all the local work */

    /* Jacobi step for the bdry points */
    for (j = 1; j<= lN; j++){
      lunew[j+(lN+2)]  = 0.25 * (hsq + lu[j+(lN+2)+1] + lu[j+(lN+2)-1] + lu[j+(lN+2)*2]);
      lunew[j+(lN+2)*lN]  = 0.25 * (hsq + lu[j+(lN+2)*lN+1] + lu[j+(lN+2)*lN-1] + lu[j+(lN+2)*(lN-1)]);
    }

    for (i = 1; i <= lN; i++){
      lunew[1+(lN+2)*i]  = 0.25 * (hsq + lu[1+(lN+2)*i+1] + lu[1+(lN+2)*i-1] + lu[1+(lN+2)*(i-1)] + lu[1+(lN+2)*(i+1)]);
      lunew[lN+(lN+2)*i]  = 0.25 * (hsq + lu[lN+(lN+2)*i+1] + lu[lN+(lN+2)*i-1] + lu[lN+(lN+2)*(i-1)] + lu[lN+(lN+2)*(i+1)]);
    }

    /* communicate ghost values */
    if (mpirank % lN < p - 1) {
      /* If not the right-most process, send/recv bdry values to the right */
      for (i = 1; i<= lN; i++){
	MPI_Isend(&(lunew[(lN+2)*(i+1)-2]), 1, MPI_DOUBLE, mpirank+1, 124, MPI_COMM_WORLD, &request_in1);
	MPI_Irecv(&(lunew[(lN+2)*i]), 1, MPI_DOUBLE, mpirank+1, 123, MPI_COMM_WORLD, &request_out1);
      }
    }
    if (mpirank % lN > 0) {
      /* If not the left-most process, send/recv bdry values to the left */
      for (i = 1; i <= lN; i++){
        MPI_Isend(&(lunew[(lN+2)*i+1]), 1, MPI_DOUBLE, mpirank-1, 123, MPI_COMM_WORLD, &request_in2);
        MPI_Irecv(&(lunew[(lN+2)*(i+1)-1]), 1, MPI_DOUBLE, mpirank-1, 124, MPI_COMM_WORLD, &request_out2);
      }
    }
    if (mpirank - lN*(lN-1) >= lN-1) {
      /* If not the top-most process, send/recv bdry values to the top */
      for (i = 1; i <= lN; i++){
        MPI_Isend(&(lunew[i]), 1, MPI_DOUBLE, mpirank+lN, 123, MPI_COMM_WORLD, &request_in3);
        MPI_Irecv(&(lunew[(lN+2)*lN+i]), 1, MPI_DOUBLE, mpirank-lN, 124, MPI_COMM_WORLD, &request_out3);
      }
    }
     if (mpirank >= lN) {
      /* If not the bottom-most process, send/recv bdry values to the bottom */
      for (i = 1; i <= lN; i++){
        MPI_Isend(&(lunew[(lN+2)*(lN-1)+i]), 1, MPI_DOUBLE, mpirank - lN, 123, MPI_COMM_WORLD, &request_in4);
        MPI_Irecv(&(lunew[i]), 1, MPI_DOUBLE, mpirank + lN, 124, MPI_COMM_WORLD, &request_out4);
      }
    }

    /* Jacobi step for the inner points */
    for (i = 2; i < lN; i++){
      for (j = 2; j< lN; j++){
	lunew[j+(lN+2)*i]  = 0.25 * (hsq + lu[j+(lN+2)*i+1] + lu[j+(lN+2)*i-1] + lu[j+(lN+2)*(i-1)] + lu[j+(lN+2)*(i+1)]);
      }
    }

    /* check if Isend/Irecv are done */
    if (mpirank < p - 1) {
      MPI_Wait(&request_out1, &status);
      MPI_Wait(&request_in1, &status);
    }
    if (mpirank > 0) {
      MPI_Wait(&request_out2, &status);
      MPI_Wait(&request_in2, &status);
    }

    /* copy newu to u using pointer flipping */
    lutemp = lu; lu = lunew; lunew = lutemp;
    if (0 == (iter % 10)) {
      gres = compute_residual(lu, lN, invhsq);
      if (0 == mpirank) {
	printf("Iter %d: Residual: %g\n", iter, gres);
      }
    }
  }

  /* Clean up */
  free(lu);
  free(lunew);

  /* timing */
  MPI_Barrier(MPI_COMM_WORLD);
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  if (0 == mpirank) {
    printf("Time elapsed is %f seconds.\n", elapsed);
  }
  MPI_Finalize();
  return 0;
}
