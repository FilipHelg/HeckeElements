/* Build the required KL context implementation in this TU, then include the
    nonzero triple checker implementation under a renamed main so we can call
    it and post-process its CSV output. */

/* Bring in core context implementation (rename its main to avoid collision). */
#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

/* When including the nonzero checker, skip its internal re-include of the
    context (we already compiled it above), and rename its main so we can call
    it directly. */
#define KL_TRIPLE_CHECK_NONZERO_SKIP_INCLUDE
#define main kl_triple_check_nonzero_included_main
#include "kl_triple_check_nonzero.c"
#undef main
#undef KL_TRIPLE_CHECK_NONZERO_SKIP_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void Usage(const char *prog) {
    printf("Usage: %s n [--dual-bin PATH] [--kl-data PATH]\n", prog);
}

int main(int argc, char *argv[]) {
    if (argc < 2) { Usage(argv[0]); return 1; }

    int n = atoi(argv[1]);
    if (n < 3 || n > 8) { fprintf(stderr, "n must be between 3 and 8.\n"); return 1; }

    /* First, run the included nonzero checker to produce the CSV file. */
    if (kl_triple_check_nonzero_included_main(argc, argv) != 0) {
        fprintf(stderr, "kl_triple_check_nonzero failed.\n");
        return 1;
    }

    /* Determine output CSV path (same default as the original tool). */
    char defaultOut[128];
    snprintf(defaultOut, sizeof(defaultOut), "kl_triple_matches_nonzero_n%d.csv", n);
    const char *outPath = defaultOut;
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) { outPath = argv[i+1]; break; }
    }

    /* Build same-left-cell pairs. */
    int count = fac(n);
    TableauEntry *entries = NULL;
    int64_t pairCount = 0;
    if (!BuildSameCellEntries(n, count, &entries, &pairCount)) {
        fprintf(stderr, "Failed to build same-cell buckets.\n");
        return 1;
    }

    /* Read CSV and collect covered pairs (x,y). */
    FILE *in = fopen(outPath, "r");
    if (!in) {
        fprintf(stderr, "Could not open %s\n", outPath);
        free(entries);
        return 1;
    }

    /* simple fixed-size hash for pairs using malloc'd array of flags is OK for small n */
    int max = count;
    unsigned char *covered = (unsigned char *)calloc((size_t)max * (size_t)max, 1);
    if (!covered) { fclose(in); free(entries); return 1; }

    char line[512];
    /* skip header */
    if (!fgets(line, sizeof(line), in)) { fclose(in); free(entries); free(covered); return 1; }
    while (fgets(line, sizeof(line), in)) {
        int w,x,y;
        if (sscanf(line, "%d,%d,%d", &w, &x, &y) >= 3) {
            if (x>y) { int t=x; x=y; y=t; }
            covered[x * max + y] = 1;
        }
    }
    fclose(in);

    /* Verify every same-cell pair has at least one covered triple */
    int missing = 0;
    for (int i = 0; i < count; ) {
        int j = i + 1;
        while (j < count && memcmp(entries[i].key.bytes, entries[j].key.bytes, sizeof(entries[i].key.bytes)) == 0) j++;
        int bsize = j - i;
        if (bsize > 1) {
            for (int a = i; a < j; a++) {
                for (int b = a+1; b < j; b++) {
                    int x = entries[a].index;
                    int y = entries[b].index;
                    if (x>y) { int t=x; x=y; y=t; }
                    if (!covered[x * max + y]) {
                        if (missing == 0) fprintf(stderr, "Missing pairs:\n");
                        fprintf(stderr, "%d,%d\n", x, y);
                        missing++;
                    }
                }
            }
        }
        i = j;
    }

    if (missing == 0) {
        printf("All %lld same-left-cell pairs have at least one nonzero matching involution.\n", (long long)pairCount);
    } else {
        printf("%d same-left-cell pairs are missing a nonzero matching involution.\n", missing);
    }

    free(covered);
    free(entries);
    return (missing == 0) ? 0 : 2;
}
