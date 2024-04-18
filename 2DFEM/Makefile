CC = gcc
CFLAG = -Ofast -fopenmp

inputdata.o: inputdata.c
	gcc -o inputdata.o -c inputdata.c
construction.o: construction.c
	gcc -o construction.o -c construction.c
solver.o: solver.c
	gcc -o solver.o -c solver.c
solver_mp.o: solver.c
	gcc $(CFLAGS) -o solver.o -c solver.c
main.o: main.c
	gcc -o main.o -c main.c
2DFEM: main.o inputdata.o construction.o solver.o
	gcc -o 2DEFEM main.o inputdata.o construction.o solver.o
2DFEM_mp: main.o inputdata.o construction.o solver_mp.o
	gcc -o 2DFEM_mp main.o inputdata.o construction.o solver_mp.o
clean:
	del -rf *.o