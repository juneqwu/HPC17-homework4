## Spring 2017: High Performance Computing 

### Assignment 4 

Problem 1: Modifications and comments are in the fixed code.

Problem 2: See plots for scalability study. When running jacobi-mpi2D, you need to add the second and the third argument. The second argument is to specify the number of grid points N (N by N matrix), where N must be a multiple of the square root of the number of MPI tasks (processes). The number of MPI tasks should be the power of 4. This elaborate set up is to make sure that each process get a local square matrix with identical size. The third argument specifies the number of maximum iterations.

Problem 3: Parallel sample sort. Used MPI collective communications, e.g. MPI_Alltoall, MPI_Alltoallv, MPI_Bcast. 

\* *Problem 2 and Problem 3 report timing ran on [STAMPEDE](https://www.tacc.utexas.edu/stampede/)*.



