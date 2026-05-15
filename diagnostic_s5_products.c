#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Forward declarations from lehmer.h */
typedef struct {
    int supportSize;
    int *support;
    long long *coeff;
} SparseHecke;

typedef struct {
    int capacity;
    int *supportSize;
    int **support;
    long long **coeff;
} DenseAccum;

typedef struct {
    int n;
    int count;
    int maxLength;
    int *lengths;
    /* ... other fields */
} FastContext;

/* Include necessary headers */
#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

#include <csv.h>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <csv_file> [num_samples]\n", argv[0]);
        return 1;
    }

    const char *csvPath = argv[1];
    int numSamples = (argc > 2) ? atoi(argv[2]) : 5;

    FILE *csv = fopen(csvPath, "r");
    if (!csv) {
        fprintf(stderr, "Cannot open %s\n", csvPath);
        return 1;
    }

    /* Skip header */
    char line[256];
    if (!fgets(line, sizeof(line), csv)) {
        fprintf(stderr, "Cannot read header\n");
        fclose(csv);
        return 1;
    }

    /* Read sample triples */
    printf("Sampling %d triples from %s:\n\n", numSamples, csvPath);
    
    int sampleCount = 0;
    while (fgets(line, sizeof(line), csv) && sampleCount < numSamples) {
        int w, x, y;
        if (sscanf(line, "%d,%d,%d", &w, &x, &y) != 3) {
            continue;
        }

        printf("Triple: w=%d, x=%d, y=%d\n", w, x, y);
        printf("  Product 1: dkH_%d * kH_%d\n", w, x);
        printf("  Product 2: dkH_%d * kH_%d\n", w, y);
        printf("\n");

        sampleCount++;
    }

    fclose(csv);
    printf("Total samples read: %d\n", sampleCount);
    return 0;
}
