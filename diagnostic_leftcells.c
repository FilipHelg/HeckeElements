
#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

#include <stdio.h>
#include <string.h>

typedef struct {
    int index;
    unsigned char P[64];
} DiagEntry;

static int CompareByP(const void *a, const void *b) {
    const DiagEntry *ea = (const DiagEntry *)a;
    const DiagEntry *eb = (const DiagEntry *)b;
    return memcmp(ea->P, eb->P, 64);
}

int main() {
    int n = 5;
    int count = 120;
    
    DiagEntry *entries = malloc(count * sizeof(DiagEntry));
    
    for (int i = 0; i < count; i++) {
        char P[8][8] = {{0}};
        char Q[8][8] = {{0}};
        char shape[8] = {0};
        RSTableaux(n, i, P, Q, shape);
        entries[i].index = i;
        memcpy(entries[i].P, P, 64);
    }
    
    qsort(entries, count, sizeof(DiagEntry), CompareByP);
    
    // Count unique P-tableaux (left cells)
    int cellCount = 1;
    for (int i = 1; i < count; i++) {
        if (memcmp(entries[i].P, entries[i-1].P, 64) != 0) {
            cellCount++;
        }
    }
    
    printf("S5: %d permutations, %d left cells\n", count, cellCount);
    
    // Count buckets with size > 1
    int multiBuckets = 0;
    int totalPairs = 0;
    for (int i = 0; i < count; ) {
        int j = i + 1;
        while (j < count && memcmp(entries[i].P, entries[j].P, 64) == 0) j++;
        int bucketSize = j - i;
        if (bucketSize > 1) {
            multiBuckets++;
            totalPairs += bucketSize * (bucketSize - 1) / 2;
        }
        i = j;
    }
    
    printf("  Buckets with >1 element: %d\n", multiBuckets);
    printf("  Total (x,y) pairs: %d\n", totalPairs);
    
    free(entries);
    return 0;
}
