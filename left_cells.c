#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "lehmer.h"

typedef struct {
    int index;
    unsigned char key[64];
} CellEntry;

static int CompareCellEntry(const void *a, const void *b) {
    const CellEntry *ea = (const CellEntry *)a;
    const CellEntry *eb = (const CellEntry *)b;
    int cmp = memcmp(ea->key, eb->key, sizeof(ea->key));
    if (cmp != 0) {
        return cmp;
    }
    return ea->index - eb->index;
}

static int BuildCellKey(int n, int index, unsigned char key[64]) {
    char P[8][8] = {{0}};
    char Q[8][8] = {{0}};
    char shape[8] = {0};

    if (!RSTableaux(n, index, P, Q, shape)) {
        return 0;
    }

    memset(key, 0, 64);
    memcpy(key, Q, 64);
    return 1;
}

static void PrintPermutation(FILE *out, int n, int index) {
    char perm[8];
    IndexToPerm(n, index, perm);
    for (int i = 0; i < n; i++) {
        fputc((int)perm[i] + '0', out);
    }
}

static int ParseInt(const char *text, int *outValue) {
    char *end = NULL;
    long value = strtol(text, &end, 10);
    if (!text[0] || !end || *end != '\0') {
        return 0;
    }
    if (value < 0 || value > 1000000L) {
        return 0;
    }
    *outValue = (int)value;
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "Usage: %s n [csv_out]\n", argv[0]);
        fprintf(stderr, "Writes one CSV row per permutation, grouped by left cell.\n");
        return 1;
    }

    int n = 0;
    if (!ParseInt(argv[1], &n)) {
        fprintf(stderr, "Invalid n: %s\n", argv[1]);
        return 1;
    }

    if (n < 1 || n > 8) {
        fprintf(stderr, "This program supports 1 <= n <= 8. Received n=%d.\n", n);
        return 1;
    }

    int count = fac(n);
    const char *outputPath = (argc >= 3) ? argv[2] : NULL;
    char defaultPath[64];
    if (!outputPath) {
        snprintf(defaultPath, sizeof(defaultPath), "S%d_left_cells.csv", n);
        outputPath = defaultPath;
    }

    CellEntry *entries = (CellEntry *)calloc((size_t)count, sizeof(CellEntry));
    if (!entries) {
        fprintf(stderr, "Out of memory while allocating %d entries.\n", count);
        return 1;
    }

    for (int i = 0; i < count; i++) {
        entries[i].index = i;
        if (!BuildCellKey(n, i, entries[i].key)) {
            fprintf(stderr, "Failed to build left-cell key for permutation %d.\n", i);
            free(entries);
            return 1;
        }
    }

    qsort(entries, (size_t)count, sizeof(CellEntry), CompareCellEntry);

    FILE *out = fopen(outputPath, "w");
    if (!out) {
        fprintf(stderr, "Could not open output file %s\n", outputPath);
        free(entries);
        return 1;
    }

    fprintf(out, "left_cell_id,cell_size,element_index,permutation\n");

    int cellCount = 0;
    int singletonCells = 0;
    int multiElementCells = 0;
    int minCellSize = count;
    int maxCellSize = 0;
    int64_t totalElements = 0;

    for (int i = 0; i < count; ) {
        int j = i + 1;
        while (j < count && memcmp(entries[i].key, entries[j].key, sizeof(entries[i].key)) == 0) {
            j++;
        }

        int cellSize = j - i;
        if (cellSize == 1) {
            singletonCells++;
        } else {
            multiElementCells++;
        }
        if (cellSize < minCellSize) {
            minCellSize = cellSize;
        }
        if (cellSize > maxCellSize) {
            maxCellSize = cellSize;
        }

        for (int k = i; k < j; k++) {
            fprintf(out, "%d,%d,%d,", cellCount, cellSize, entries[k].index);
            PrintPermutation(out, n, entries[k].index);
            fputc('\n', out);
            totalElements++;
        }

        cellCount++;
        i = j;
    }

    fclose(out);
    free(entries);

    printf("Left-cell CSV written to: %s\n", outputPath);
    printf("S_%d summary\n", n);
    printf("  Total permutations: %d\n", count);
    printf("  Left cells: %d\n", cellCount);
    printf("  Singleton cells: %d\n", singletonCells);
    printf("  Multi-element cells: %d\n", multiElementCells);
    printf("  Smallest cell size: %d\n", minCellSize);
    printf("  Largest cell size: %d\n", maxCellSize);
    printf("  Elements written: %lld\n", (long long)totalElements);

    return 0;
}