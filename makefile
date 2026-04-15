CC = gcc
CFLAGS = -std=c11 -Ofast -march=native -Wall
LDFLAGS = -lm

.PHONY: all clean

all: hecke kl_structure_poc kl_structure_bulk kl_structure_bulk_opt

# Build the hecke binary linking against lehmer.o
hecke: hecke.o lehmer.o laurent.o
	$(CC) -o hecke hecke.o lehmer.o laurent.o $(LDFLAGS)

kl_structure_poc: kl_structure_poc.o lehmer.o laurent.o
	$(CC) -o kl_structure_poc kl_structure_poc.o lehmer.o laurent.o $(LDFLAGS)

kl_structure_bulk: kl_structure_bulk.o lehmer.o laurent.o
	$(CC) -o kl_structure_bulk kl_structure_bulk.o lehmer.o laurent.o $(LDFLAGS)

kl_structure_bulk_opt: kl_structure_bulk_opt.o lehmer.o
	$(CC) -o kl_structure_bulk_opt kl_structure_bulk_opt.o lehmer.o $(LDFLAGS)

hecke.o: hecke.c lehmer.h laurent.h
	$(CC) $(CFLAGS) -c hecke.c -o hecke.o

kl_structure_poc.o: kl_structure_poc.c lehmer.h laurent.h
	$(CC) $(CFLAGS) -c kl_structure_poc.c -o kl_structure_poc.o

kl_structure_bulk.o: kl_structure_bulk.c lehmer.h laurent.h
	$(CC) $(CFLAGS) -c kl_structure_bulk.c -o kl_structure_bulk.o

kl_structure_bulk_opt.o: kl_structure_bulk_opt.c lehmer.h
	$(CC) $(CFLAGS) -c kl_structure_bulk_opt.c -o kl_structure_bulk_opt.o

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
	-del /Q lehmer.exe 2>nul || if exist lehmer del /Q lehmer 2>nul
	-del /Q *.o *.leh 2>nul || echo Files not found
	-del /Q *.o *.leh 2>nul || echo Files not found