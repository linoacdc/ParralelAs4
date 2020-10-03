CC=gcc
CFLAGS=-fopenmp

sequential:
	$(CC) $(FLAGS) -o Sequential.exe Sequential.c

openmp:
	$(CC) $(CFLAGS) -o ompProject.exe ompProject.c

clean:
	rm -rf *.exe
