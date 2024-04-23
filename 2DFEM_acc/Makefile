CC = gcc
CFLAG = -std=c11

inputdata.o: inputdata.c
	$(CC) -o inputdata.o -c inputdata.c
construction.o: construction.c
	$(CC) -o construction.o -c construction.c
solver.o: solver.c
	$(CC) -o solver.o -c solver.c
solver_mp.o: solver.c
	$(CC) $(CFLAG) -o solver_mp.o -c solver.c
main.o: main.c
	$(CC) -o main.o -c main.c
2DFEM: main.o inputdata.o construction.o solver.o
	$(CC) -o 2DEFEM main.o inputdata.o construction.o solver.o
2DFEM_mp: main.o inputdata.o construction.o solver_mp.o
	$(CC) $(CFLAG) -o 2DFEM_mp main.o inputdata.o construction.o solver_mp.o
clean:
	del -rf *.o