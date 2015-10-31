all: prog

prog: program.o 
	mpicc -o cp_mpi parallelized_program.o

program.o: parallelized_program.c
	mpicc -c parallelized_program.c