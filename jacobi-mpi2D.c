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
  int i, j;
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
  int mpirank, i, j, p, N, lN, iter, max_iters, lp;
  MPI_Status status, status1;

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
  lN = (int) N / sqrt(p);
  lp = (int) sqrt(p);
  if ((N % (int) sqrt(p) != 0 ) && mpirank == 0 ) {
    printf("N: %d, local N: %d\n", N, lN);
    printf("Exiting. N must be a multiple of the square root of p\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  /* timing */
  MPI_Barrier(MPI_COMM_WORLD);
  timestamp_type time1, time2;
  get_timestamp(&time1);

  /* Allocation of vectors, including left/upper and right/lower ghost points */
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

    /* Jacobi step for local points */
    for (i = 1; i <= lN; i++){
      for (j = 1; j<= lN; j++){
	lunew[j+(lN+2)*i]  = 0.25 * (hsq + lu[j+(lN+2)*i+1] + lu[j+(lN+2)*i-1] + lu[j+(lN+2)*(i-1)] + lu[j+(lN+2)*(i+1)]);
      }
    }

    /* communicate ghost values */
    if (mpirank % lp < (lp - 1) ) {
      /* If not the right-most process, send/recv bdry values to the right */
      for (i = 1; i<=lN; i++){
	MPI_Send(&(lunew[(lN+2)*(i+1)-2]), 1, MPI_DOUBLE, mpirank+1, 124, MPI_COMM_WORLD);
	MPI_Recv(&(lunew[(lN+2)*i]), 1, MPI_DOUBLE, mpirank+1, 123, MPI_COMM_WORLD, &status);
      }
    }
    if (mpirank % lp > 0 ) {
      /* If not the left-most process, send/recv bdry values to the left */
      for (i = 1; i <= lN; i++){
        MPI_Send(&(lunew[(lN+2)*i+1]), 1, MPI_DOUBLE, mpirank-1, 123, MPI_COMM_WORLD);
        MPI_Recv(&(lunew[(lN+2)*(i+1)-1]), 1, MPI_DOUBLE, mpirank-1, 124, MPI_COMM_WORLD, &status1);
      }
    }
    if (lp*(lp-1) - mpirank > 0) {
      /* If not the top-most process, send/recv bdry values to the top */
      for (i = 1; i <= lN; i++){
        MPI_Send(&(lunew[i]), 1, MPI_DOUBLE, mpirank+lp, 125, MPI_COMM_WORLD);
        MPI_Recv(&(lunew[(lN+2)*lN+i]), 1, MPI_DOUBLE, mpirank + lp, 126, MPI_COMM_WORLD, &status1);
      }
    }
     if (mpirank >= lp ) {
      /* If not the bottom-most process, send/recv bdry values to the bottom */
      for (i = 1; i <= lN; i++){
        MPI_Send(&(lunew[(lN+2)*(lN-1)+i]), 1, MPI_DOUBLE, mpirank - lp, 126, MPI_COMM_WORLD);
        MPI_Recv(&(lunew[i]), 1, MPI_DOUBLE, mpirank - lp, 125, MPI_COMM_WORLD, &status1);
      }
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
