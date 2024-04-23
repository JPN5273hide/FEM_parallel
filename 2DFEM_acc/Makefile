CC = gcc
CFLAG = -std=c11 -Ofast -fopenmp 

inputdata.o: inputdata.c
	$(CC) -o inputdata.o -c inputdata.c
construction.o: construction.c
	$(CC) -o construction.o -c construction.c
solver.o: solver.c
	$(CC) -o solver.o -c solver.c
solver_acc.o: solver.c
	$(CC) $(CFLAG) -o solver_acc.o -c solver.c -lm
main.o: main.c
	$(CC) -o main.o -c main.c
2DFEM: main.o inputdata.o construction.o solver.o
	$(CC) -o 2DEFEM main.o inputdata.o construction.o solver.o
2DFEM_acc: main.o inputdata.o construction.o solver_mp.o
	$(CC) $(CFLAG) -o 2DFEM_acc main.o inputdata.o construction.o solver_acc.o -lm
clean:
	del -rf *.o