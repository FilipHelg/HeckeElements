CC = gcc
CFLAGS = -std=c11 -Ofast -march=native -Wall
LDFLAGS = -lm

.PHONY: all clean

all: hecke

# Build the hecke binary linking against lehmer.o
hecke: hecke.o lehmer.o laurent.o
	$(CC) -o hecke hecke.o lehmer.o laurent.o $(LDFLAGS)

hecke.o: hecke.c lehmer.h laurent.h
	$(CC) $(CFLAGS) -c hecke.c -o hecke.o

laurent.o: laurent.c laurent.h
	$(CC) $(CFLAGS) -c laurent.c -o laurent.o

lehmer.o: lehmer.c lehmer.h
	$(CC) $(CFLAGS) -c lehmer.c -o lehmer.o

clean:
	@echo Cleaning build artifacts (Windows)...
	-del /Q hecke.exe 2>nul || if exist hecke del /Q hecke 2>nul
	-del /Q lehmer.exe 2>nul || if exist lehmer del /Q lehmer 2>nul
	-del /Q *.o *.leh 2>nul || echo Files not found
	-del /Q *.o *.leh 2>nul || echo Files not found