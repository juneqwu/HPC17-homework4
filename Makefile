CC = mpicc
FLAGS=-O3
LDFLAGS = -lm
EXECS = mpi_solved1 mpi_solved2 mpi_solved3 mpi_solved4 mpi_solved5 mpi_solved6 mpi_solved7 jacobi-mpi2D ssort

all: ${EXECS}

mpi_solved1: mpi_solved1.c
	${CC} ${FLAGS} $^ -o mpi_solved1 ${LDFLAGS}

mpi_solved2: mpi_solved2.c
	${CC} ${FLAGS} $^ -o mpi_solved2 ${LDFLAGS}

mpi_solved3: mpi_solved3.c
	${CC} ${FLAGS} $^ -o mpi_solved3 ${LDFLAGS}

mpi_solved4: mpi_solved4.c
	${CC} ${FLAGS} $^ -o mpi_solved4 ${LDFLAGS}

mpi_solved5: mpi_solved5.c
	${CC} ${FLAGS} $^ -o mpi_solved5 ${LDFLAGS}

mpi_solved6: mpi_solved6.c
	${CC} ${FLAGS} $^ -o mpi_solved6 ${LDFLAGS}

mpi_solved7: mpi_solved7.c
	${CC} ${FLAGS} $^ -o mpi_solved7 ${LDFLAGS}

jacobi-mpi2D: jacobi-mpi2D.c
	${CC} ${FLAGS} $^ -o jacobi-mpi2D ${LDFLAGS}

ssort: ssort.c
	${CC} ${FLAGS} $^ -o ssort ${LDFLAGS}

clean:
	rm -f ${EXECS}
