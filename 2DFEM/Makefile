CC = gcc
CFLAG = -fopenmp

inputdata.o: inputdata.c
	$(CC) -o inputdata.o -c inputdata.c
construction.o: construction.c
	$(CC) -o construction.o -c construction.c
solver.o: solver.c
	$(CC) $(CFLAG) -o solver.o -c solver.c
main.o: main.c
	$(CC) -o main.o -c main.c
2DFEM: main.o inputdata.o construction.o solver.o
	$(CC) $(CFLAG) -o 2DFEM main.o inputdata.o construction.o solver.o
clean:
	del -rf *.o