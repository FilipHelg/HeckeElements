CC = gcc
CFLAGS = -std=c11 -Ofast -march=native -Wall
LDFLAGS = -lm
OMPFLAGS = -fopenmp

.PHONY: all clean

all: hecke kl_structure_poc kl_structure_bulk kl_structure_bulk_opt kl_structure_bulk_opt_parallel kl_runtime_sample_estimate

# Build the hecke binary linking against lehmer.o
hecke: hecke.o lehmer.o laurent.o
	$(CC) -o hecke hecke.o lehmer.o laurent.o $(LDFLAGS)

kl_structure_poc: kl_structure_poc.o lehmer.o laurent.o
	$(CC) -o kl_structure_poc kl_structure_poc.o lehmer.o laurent.o $(LDFLAGS)

kl_structure_bulk: kl_structure_bulk.o lehmer.o laurent.o
	$(CC) -o kl_structure_bulk kl_structure_bulk.o lehmer.o laurent.o $(LDFLAGS)

kl_structure_bulk_opt: kl_structure_bulk_opt.o lehmer.o
	$(CC) -o kl_structure_bulk_opt kl_structure_bulk_opt.o lehmer.o $(LDFLAGS)

kl_structure_bulk_opt_parallel: kl_structure_bulk_opt_parallel.o lehmer.o
	$(CC) -o kl_structure_bulk_opt_parallel kl_structure_bulk_opt_parallel.o lehmer.o $(LDFLAGS) $(OMPFLAGS)

kl_runtime_sample_estimate: kl_runtime_sample_estimate.o lehmer.o
	$(CC) -o kl_runtime_sample_estimate kl_runtime_sample_estimate.o lehmer.o $(LDFLAGS) $(OMPFLAGS)

hecke.o: hecke.c lehmer.h laurent.h
	$(CC) $(CFLAGS) -c hecke.c -o hecke.o

kl_structure_poc.o: kl_structure_poc.c lehmer.h laurent.h
	$(CC) $(CFLAGS) -c kl_structure_poc.c -o kl_structure_poc.o

kl_structure_bulk.o: kl_structure_bulk.c lehmer.h laurent.h
	$(CC) $(CFLAGS) -c kl_structure_bulk.c -o kl_structure_bulk.o

kl_structure_bulk_opt.o: kl_structure_bulk_opt.c lehmer.h
	$(CC) $(CFLAGS) -c kl_structure_bulk_opt.c -o kl_structure_bulk_opt.o

kl_structure_bulk_opt_parallel.o: kl_structure_bulk_opt_parallel.c lehmer.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -c kl_structure_bulk_opt_parallel.c -o kl_structure_bulk_opt_parallel.o

kl_runtime_sample_estimate.o: kl_runtime_sample_estimate.c lehmer.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -c kl_runtime_sample_estimate.c -o kl_runtime_sample_estimate.o

laurent.o: laurent.c laurent.h
	$(CC) $(CFLAGS) -c laurent.c -o laurent.o

lehmer.o: lehmer.c lehmer.h
	$(CC) $(CFLAGS) -c lehmer.c -o lehmer.o

clean:
	@echo Cleaning build artifacts (Windows)...
	-del /Q hecke.exe 2>nul || if exist hecke del /Q hecke 2>nul
	-del /Q kl_structure_poc.exe 2>nul || if exist kl_structure_poc del /Q kl_structure_poc 2>nul
	-del /Q kl_structure_bulk.exe 2>nul || if exist kl_structure_bulk del /Q kl_structure_bulk 2>nul
	-del /Q kl_structure_bulk_opt.exe 2>nul || if exist kl_structure_bulk_opt del /Q kl_structure_bulk_opt 2>nul
	-del /Q kl_structure_bulk_opt_parallel.exe 2>nul || if exist kl_structure_bulk_opt_parallel del /Q kl_structure_bulk_opt_parallel 2>nul
	-del /Q kl_runtime_sample_estimate.exe 2>nul || if exist kl_runtime_sample_estimate del /Q kl_runtime_sample_estimate 2>nul
	-del /Q lehmer.exe 2>nul || if exist lehmer del /Q lehmer 2>nul
	-del /Q *.o *.leh 2>nul || echo Files not found
	-del /Q *.o *.leh 2>nul || echo Files not found