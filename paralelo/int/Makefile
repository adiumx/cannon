PROGNAME := check_memory

main: memory.h memory_parallel.h
	mpicc $(PROGNAME).c memory.c memory_parallel.c -o $(PROGNAME) -lm

clean:
	rm -f $(PROGNAME) *.o
