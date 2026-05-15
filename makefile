CC = gcc
CFLAGS = -std=c11 -Ofast -march=native -Wall
LDFLAGS = -lm
OMPFLAGS = -fopenmp

.PHONY: all clean

all: hecke kl_structure_poc kl_structure_bulk kl_structure_bulk_opt kl_structure_bulk_opt_parallel dual_kl_poc dual_kl_bulk dual_kl_reader kl_triple_check kl_triple_check_cached kl_triple_check_nonzero kl_triple_check_q_tableau kl_triple_check_reverse_mult validate_old_products kl_involution_check kl_involution_check_P kl_involution_check_Q_right right_preorder_involutions left_cells basis_product #kl_runtime_sample_estimate

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

kl_triple_check: kl_triple_check.o lehmer.o
	$(CC) -o kl_triple_check kl_triple_check.o lehmer.o $(LDFLAGS) $(OMPFLAGS)

kl_triple_check_cached: kl_triple_check_cached.o lehmer.o
	$(CC) -o kl_triple_check_cached kl_triple_check_cached.o lehmer.o $(LDFLAGS) $(OMPFLAGS)

kl_triple_check_nonzero: kl_triple_check_nonzero.o lehmer.o
	$(CC) -o kl_triple_check_nonzero kl_triple_check_nonzero.o lehmer.o $(LDFLAGS) $(OMPFLAGS)

kl_involution_check: kl_involution_check.o lehmer.o
	$(CC) -o kl_involution_check kl_involution_check.o lehmer.o $(LDFLAGS) $(OMPFLAGS)

kl_involution_check_P: kl_involution_check_P.c lehmer.o
	$(CC) -std=c11 -Ofast -march=native -Wall -fopenmp -o kl_involution_check_P kl_involution_check_P.c lehmer.o $(LDFLAGS)

kl_involution_check_Q_right: kl_involution_check_Q_right.c lehmer.o
	$(CC) -std=c11 -Ofast -march=native -Wall -fopenmp -o kl_involution_check_Q_right kl_involution_check_Q_right.c lehmer.o $(LDFLAGS)

right_preorder_involutions: right_preorder_involutions.c lehmer.o
	$(CC) -std=c11 -Ofast -march=native -Wall -fopenmp -o right_preorder_involutions right_preorder_involutions.c lehmer.o $(LDFLAGS)

left_cells: left_cells.o lehmer.o
	$(CC) -o left_cells left_cells.o lehmer.o $(LDFLAGS)

basis_product: basis_product.c lehmer.o
	$(CC) -std=c11 -Ofast -march=native -Wall $(OMPFLAGS) -o basis_product basis_product.c lehmer.o $(LDFLAGS)

kl_triple_check_q_tableau: kl_triple_check_q_tableau.o lehmer.o
	$(CC) -o kl_triple_check_q_tableau kl_triple_check_q_tableau.o lehmer.o $(LDFLAGS) $(OMPFLAGS)

kl_triple_check_reverse_mult: kl_triple_check_reverse_mult.o lehmer.o
	$(CC) -o kl_triple_check_reverse_mult kl_triple_check_reverse_mult.o lehmer.o $(LDFLAGS) $(OMPFLAGS)

validate_old_products: validate_old_products.o lehmer.o
	$(CC) -o validate_old_products validate_old_products.o lehmer.o $(LDFLAGS) $(OMPFLAGS)

dual_kl_poc: dual_kl_poc.c lehmer.c laurent.c lehmer.h laurent.h
	$(CC) -std=c11 -O0 -g -Wall -o dual_kl_poc dual_kl_poc.c lehmer.c laurent.c $(LDFLAGS)

dual_kl_bulk: dual_kl_bulk.c lehmer.h
	$(CC) -std=c11 -Ofast -march=native -Wall $(OMPFLAGS) -o dual_kl_bulk dual_kl_bulk.c $(LDFLAGS)

dual_kl_reader: dual_kl_reader.c lehmer.c lehmer.h
	$(CC) -std=c11 -O0 -g -Wall -o dual_kl_reader dual_kl_reader.c lehmer.c $(LDFLAGS)

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

kl_triple_check.o: kl_triple_check.c lehmer.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -c kl_triple_check.c -o kl_triple_check.o

kl_triple_check_cached.o: kl_triple_check_cached.c lehmer.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -c kl_triple_check_cached.c -o kl_triple_check_cached.o

kl_triple_check_q_tableau.o: kl_triple_check_q_tableau.c lehmer.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -c kl_triple_check_q_tableau.c -o kl_triple_check_q_tableau.o

validate_old_products.o: validate_old_products.c lehmer.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -c validate_old_products.c -o validate_old_products.o

kl_triple_check_nonzero.o: kl_triple_check_nonzero.c lehmer.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -c kl_triple_check_nonzero.c -o kl_triple_check_nonzero.o

kl_involution_check.o: kl_involution_check.c lehmer.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -c kl_involution_check.c -o kl_involution_check.o

left_cells.o: left_cells.c lehmer.h
	$(CC) $(CFLAGS) -c left_cells.c -o left_cells.o

dual_kl_poc.o: dual_kl_poc.c lehmer.h laurent.h
	$(CC) $(CFLAGS) -c dual_kl_poc.c -o dual_kl_poc.o

dual_kl_reader.o: dual_kl_reader.c lehmer.h
	$(CC) $(CFLAGS) -c dual_kl_reader.c -o dual_kl_reader.o

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
	-del /Q kl_involution_check_P.exe 2>nul || if exist kl_involution_check_P del /Q kl_involution_check_P 2>nul
	-del /Q kl_triple_check_q_tableau.exe 2>nul || if exist kl_triple_check_q_tableau del /Q kl_triple_check_q_tableau 2>nul
	-del /Q kl_triple_check.exe 2>nul || if exist kl_triple_check del /Q kl_triple_check 2>nul
	-del /Q kl_triple_check_cached.exe 2>nul || if exist kl_triple_check_cached del /Q kl_triple_check_cached 2>nul
	-del /Q kl_triple_check_nonzero.exe 2>nul || if exist kl_triple_check_nonzero del /Q kl_triple_check_nonzero 2>nul	-del /Q kl_triple_check_reverse_mult.exe 2>nul || if exist kl_triple_check_reverse_mult del /Q kl_triple_check_reverse_mult 2>nul	-del /Q dual_kl_poc.exe 2>nul || if exist dual_kl_poc del /Q dual_kl_poc 2>nul
	-del /Q left_cells.exe 2>nul || if exist left_cells del /Q left_cells 2>nul
	-del /Q kl_triple_check_nonzero.exe 2>nul || if exist kl_triple_check_nonzero del /Q kl_triple_check_nonzero 2>nul	-del /Q kl_triple_check_reverse_mult.exe 2>nul || if exist kl_triple_check_reverse_mult del /Q kl_triple_check_reverse_mult 2>nul	-del /Q dual_kl_poc.exe 2>nul || if exist dual_kl_poc del /Q dual_kl_poc 2>nul
	-del /Q dual_kl_bulk.exe 2>nul || if exist dual_kl_bulk del /Q dual_kl_bulk 2>nul
	-del /Q dual_kl_reader.exe 2>nul || if exist dual_kl_reader del /Q dual_kl_reader 2>nul
	-del /Q lehmer.exe 2>nul || if exist lehmer del /Q lehmer 2>nul
	-del /Q *.o *.leh 2>nul || echo Files not found
	-del /Q *.o *.leh 2>nul || echo Files not found