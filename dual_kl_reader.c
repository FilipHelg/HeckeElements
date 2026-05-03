#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lehmer.h"

typedef struct {
    char magic[8];
    uint32_t version;
    uint32_t n;
    uint32_t count;
    uint32_t involutionCount;
    uint32_t recordFormat;
    uint32_t reserved[3];
} DualBulkHeader;

typedef struct {
    int8_t exponent;
    int32_t coeff;
} PolyTerm;

static int FacLocal(int n) {
    int value = 1;
    for (int i = 2; i <= n; i++) {
        value *= i;
    }
    return value;
}

static int ReadU32(FILE *in, uint32_t *value) {
    return fread(value, sizeof(*value), 1, in) == 1;
}

static int ReadI32(FILE *in, int32_t *value) {
    return fread(value, sizeof(*value), 1, in) == 1;
}

static int ReadI8(FILE *in, int8_t *value) {
    return fread(value, sizeof(*value), 1, in) == 1;
}

static void PrintPermutation(int n, int index) {
    char perm[8] = {0};
    IndexToPerm(n, index, perm);
    for (int i = 0; i < n; i++) {
        printf("%d", (int)perm[i]);
    }
}

static void PrintUsage(const char *prog) {
    printf("Usage:\n");
    printf("  %s file.bin [--x INDEX]\n", prog);
    printf("\n");
    printf("Prints the header and decoded dual KL records from a bulk binary cache.\n");
    printf("Use --x to print only one target record.\n");
}

static int ParseArgs(int argc, char *argv[], const char **path, int *filterX) {
    if (argc < 2 || argc > 4) {
        return 0;
    }

    *path = argv[1];
    *filterX = -1;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--x") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--x requires an index.\n");
                return 0;
            }
            *filterX = atoi(argv[++i]);
            continue;
        }

        fprintf(stderr, "Unknown option: %s\n", argv[i]);
        return 0;
    }

    return 1;
}

static int DumpRecord(FILE *in, uint32_t n, uint32_t count, int wantedX) {
    (void)count;
    uint32_t x = 0;
    uint32_t supportCount = 0;

    if (!ReadU32(in, &x) || !ReadU32(in, &supportCount)) {
        return 0;
    }

    int printRecord = (wantedX < 0) || ((uint32_t)wantedX == x);
    if (printRecord) {
        printf("x=%u ", x);
        printf("perm=");
        PrintPermutation((int)n, (int)x);
        printf(" support=%u\n", supportCount);
    }

    for (uint32_t s = 0; s < supportCount; s++) {
        uint32_t z = 0;
        uint8_t termCount = 0;
        if (!ReadU32(in, &z) || fread(&termCount, sizeof(termCount), 1, in) != 1) {
            return 0;
        }

        PolyTerm *terms = NULL;
        if (termCount > 0) {
            terms = (PolyTerm *)calloc(termCount, sizeof(PolyTerm));
            if (!terms) {
                return 0;
            }

            for (uint8_t t = 0; t < termCount; t++) {
                if (!ReadI8(in, &terms[t].exponent) || !ReadI32(in, &terms[t].coeff)) {
                    free(terms);
                    return 0;
                }
            }
        }

        if (printRecord) {
            printf("  z=%u ", z);
            PrintPermutation((int)n, (int)z);
            printf(" : ");
            if (termCount == 0) {
                printf("0");
            } else {
                for (uint8_t t = 0; t < termCount; t++) {
                    if (t > 0) {
                        printf(" + ");
                    }
                    printf("%d*v^%d", terms[t].coeff, (int)terms[t].exponent);
                }
            }
            printf("\n");
        }

        free(terms);
    }

    return 1;
}

int main(int argc, char *argv[]) {
    const char *path = NULL;
    int filterX = -1;
    if (!ParseArgs(argc, argv, &path, &filterX)) {
        PrintUsage(argv[0]);
        return 1;
    }

    FILE *in = fopen(path, "rb");
    if (!in) {
        fprintf(stderr, "Could not open %s\n", path);
        return 1;
    }

    DualBulkHeader header;
    if (fread(&header, sizeof(header), 1, in) != 1) {
        fprintf(stderr, "Failed to read header.\n");
        fclose(in);
        return 1;
    }

    printf("magic=%.7s\n", header.magic);
    printf("version=%u\n", header.version);
    printf("n=%u\n", header.n);
    printf("count=%u\n", header.count);
    printf("involutionCount=%u\n", header.involutionCount);
    printf("recordFormat=%u\n", header.recordFormat);

    if (memcmp(header.magic, "DKBULK1", 7) != 0) {
        fprintf(stderr, "Warning: unexpected magic string.\n");
    }

    if (header.n < 3 || header.n > 8) {
        fprintf(stderr, "Warning: n looks out of range for this cache.\n");
    }

    int totalCount = FacLocal((int)header.n);
    if ((int)header.count != totalCount) {
        fprintf(stderr, "Warning: count=%u does not match n!=%d.\n", header.count, totalCount);
    }

    printf("--- Records ---\n");
    for (uint32_t i = 0; i < header.involutionCount; i++) {
        if (!DumpRecord(in, header.n, header.count, filterX)) {
            fprintf(stderr, "Failed while decoding record %u.\n", i);
            fclose(in);
            return 1;
        }
    }

    fclose(in);
    return 0;
}