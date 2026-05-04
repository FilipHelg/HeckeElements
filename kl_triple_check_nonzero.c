#ifndef KL_TRIPLE_CHECK_NONZERO_SKIP_INCLUDE
#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main
#endif

#include <time.h>

typedef struct {
    unsigned char bytes[72];
} TableauKey;

typedef struct {
    int index;
    TableauKey key;
} TableauEntry;

static int IndexToLengthSafeLocal(int n, int index) {
    return IndexToLength(n, index);
}

typedef struct {
    int index;
    Poly57 poly;
} CompactTerm;

typedef struct {
    int supportSize;
    int *support;
    Poly57 *coeff;
} CompactHecke;

typedef struct {
    SparseHecke *records;
    unsigned char *present;
    int count;
} DualCache;

typedef struct {
    long buildKLElemCount;
    double buildKLElemTime;
    long multiplyCount;
    double multiplyTime;
    long compactFromSparseCount;
    double compactFromSparseTime;
    long zeroProductSkipCount;
    double zeroProductSkipTime;
    long compactEqualsCount;
    double compactEqualsTime;
    long writeTripleCount;
    double writeTripleTime;
} CompareBenchStats;

static void BuildTableauKey(int n, int index, TableauKey *key) {
    memset(key->bytes, 0, sizeof(key->bytes));

    char P[8][8] = {{0}};
    char Q[8][8] = {{0}};
    char shape[8] = {0};
    (void)RSTableaux(n, index, P, Q, shape);

    memcpy(key->bytes, P, sizeof(P));
    memcpy(key->bytes + sizeof(P), shape, sizeof(shape));
}

static int CompareTableauEntry(const void *a, const void *b) {
    const TableauEntry *ea = (const TableauEntry *)a;
    const TableauEntry *eb = (const TableauEntry *)b;
    int cmp = memcmp(ea->key.bytes, eb->key.bytes, sizeof(ea->key.bytes));
    if (cmp != 0) {
        return cmp;
    }
    return ea->index - eb->index;
}

static int BuildSameCellEntries(int n, int count, TableauEntry **entriesOut, int64_t *pairCountOut) {
    TableauEntry *entries = (TableauEntry *)calloc((size_t)count, sizeof(TableauEntry));
    if (!entries) {
        return 0;
    }

    for (int i = 0; i < count; i++) {
        entries[i].index = i;
        BuildTableauKey(n, i, &entries[i].key);
    }

    qsort(entries, (size_t)count, sizeof(TableauEntry), CompareTableauEntry);

    int64_t pairCount = 0;
    for (int i = 0; i < count; ) {
        int j = i + 1;
        while (j < count && memcmp(entries[i].key.bytes, entries[j].key.bytes, sizeof(entries[i].key.bytes)) == 0) {
            j++;
        }

        int bucketSize = j - i;
        if (bucketSize > 1) {
            pairCount += (int64_t)bucketSize * (int64_t)(bucketSize - 1) / 2;
        }

        i = j;
    }

    *entriesOut = entries;
    *pairCountOut = pairCount;
    return 1;
}

static void CompactFree(CompactHecke *h) {
    free(h->support);
    free(h->coeff);
    h->support = NULL;
    h->coeff = NULL;
    h->supportSize = 0;
}

static int CompareCompactTerm(const void *a, const void *b) {
    const CompactTerm *ta = (const CompactTerm *)a;
    const CompactTerm *tb = (const CompactTerm *)b;
    return ta->index - tb->index;
}

static int CompactFromSparse(const SparseHecke *src, CompactHecke *dst) {
    dst->supportSize = src->supportSize;
    dst->support = NULL;
    dst->coeff = NULL;

    if (src->supportSize <= 0) {
        return 1;
    }

    CompactTerm *terms = (CompactTerm *)calloc((size_t)src->supportSize, sizeof(CompactTerm));
    if (!terms) {
        return 0;
    }

    for (int i = 0; i < src->supportSize; i++) {
        int idx = src->support[i];
        terms[i].index = idx;
        terms[i].poly = src->coeff[idx];
    }

    qsort(terms, (size_t)src->supportSize, sizeof(CompactTerm), CompareCompactTerm);

    dst->support = (int *)calloc((size_t)src->supportSize, sizeof(int));
    dst->coeff = (Poly57 *)calloc((size_t)src->supportSize, sizeof(Poly57));
    if (!dst->support || !dst->coeff) {
        free(terms);
        CompactFree(dst);
        return 0;
    }

    for (int i = 0; i < src->supportSize; i++) {
        dst->support[i] = terms[i].index;
        dst->coeff[i] = terms[i].poly;
    }

    free(terms);
    return 1;
}

static int CompactEquals(const CompactHecke *a, const CompactHecke *b) {
    if (a->supportSize != b->supportSize) {
        return 0;
    }

    for (int i = 0; i < a->supportSize; i++) {
        if (a->support[i] != b->support[i]) {
            return 0;
        }
        if (memcmp(&a->coeff[i], &b->coeff[i], sizeof(Poly57)) != 0) {
            return 0;
        }
    }

    return 1;
}

static void DualCacheFree(DualCache *cache) {
    if (!cache) {
        return;
    }

    if (cache->records && cache->present) {
        for (int i = 0; i < cache->count; i++) {
            if (cache->present[i]) {
                SparseFree(&cache->records[i]);
            }
        }
    }

    free(cache->records);
    free(cache->present);
    cache->records = NULL;
    cache->present = NULL;
    cache->count = 0;
}

typedef struct {
    char magic[8];
    uint32_t version;
    uint32_t n;
    uint32_t count;
    uint32_t involutionCount;
    uint32_t recordFormat;
    uint32_t reserved[3];
} DualBulkHeader;

static int ReadU32(FILE *in, uint32_t *value) {
    return fread(value, sizeof(*value), 1, in) == 1;
}

static int ReadI32(FILE *in, int32_t *value) {
    return fread(value, sizeof(*value), 1, in) == 1;
}

static int ReadI8(FILE *in, int8_t *value) {
    return fread(value, sizeof(*value), 1, in) == 1;
}

static int LoadDualCache(const char *path, int count, DualCache *cache, int *loadedNOut) {
    FILE *in = fopen(path, "rb");
    if (!in) {
        fprintf(stderr, "Could not open dual KL cache %s\n", path);
        return 0;
    }

    DualBulkHeader header;
    if (fread(&header, sizeof(header), 1, in) != 1) {
        fprintf(stderr, "Failed to read dual KL header from %s\n", path);
        fclose(in);
        return 0;
    }

    if (memcmp(header.magic, "DKBULK1", 7) != 0) {
        fprintf(stderr, "Warning: unexpected dual KL cache magic.\n");
    }

    if ((int)header.count != count) {
        fprintf(stderr, "Dual KL cache count=%u does not match n!=%d.\n", header.count, count);
        fclose(in);
        return 0;
    }

    cache->records = (SparseHecke *)calloc((size_t)count, sizeof(SparseHecke));
    cache->present = (unsigned char *)calloc((size_t)count, sizeof(unsigned char));
    cache->count = count;
    if (!cache->records || !cache->present) {
        fprintf(stderr, "Out of memory while loading dual KL cache.\n");
        DualCacheFree(cache);
        fclose(in);
        return 0;
    }

    for (uint32_t rec = 0; rec < header.involutionCount; rec++) {
        uint32_t x = 0;
        uint32_t supportCount = 0;
        if (!ReadU32(in, &x) || !ReadU32(in, &supportCount)) {
            fprintf(stderr, "Failed while reading dual KL record %u.\n", rec);
            DualCacheFree(cache);
            fclose(in);
            return 0;
        }

        if (x >= (uint32_t)count) {
            fprintf(stderr, "Dual KL record index %u out of range.\n", x);
            DualCacheFree(cache);
            fclose(in);
            return 0;
        }

        if (cache->present[x]) {
            fprintf(stderr, "Duplicate dual KL record for index %u.\n", x);
            DualCacheFree(cache);
            fclose(in);
            return 0;
        }

        if (!SparseInit(&cache->records[x], count)) {
            fprintf(stderr, "Out of memory allocating sparse record for %u.\n", x);
            DualCacheFree(cache);
            fclose(in);
            return 0;
        }
        cache->present[x] = 1;

        for (uint32_t s = 0; s < supportCount; s++) {
            uint32_t z = 0;
            uint8_t termCount = 0;
            if (!ReadU32(in, &z) || fread(&termCount, sizeof(termCount), 1, in) != 1) {
                fprintf(stderr, "Failed while reading dual KL support term.\n");
                DualCacheFree(cache);
                fclose(in);
                return 0;
            }

            Poly57 poly;
            PolyZero(&poly);
            for (uint8_t t = 0; t < termCount; t++) {
                int8_t exponent = 0;
                int32_t coeff = 0;
                if (!ReadI8(in, &exponent) || !ReadI32(in, &coeff)) {
                    fprintf(stderr, "Failed while reading dual KL polynomial term.\n");
                    DualCacheFree(cache);
                    fclose(in);
                    return 0;
                }

                int pos = 28 + (int)exponent;
                if (pos >= 0 && pos < 57) {
                    poly.coeff[pos] += (int)coeff;
                    if (poly.minPos > pos) poly.minPos = pos;
                    if (poly.maxPos < pos) poly.maxPos = pos;
                }
            }
            PolyNormalize(&poly);
            SparseAddPoly(&cache->records[x], (int)z, &poly);
        }
    }

    if (loadedNOut) {
        *loadedNOut = (int)header.n;
    }

    fclose(in);
    return 1;
}

static void PrintTripleUsage(const char *prog) {
    printf("Usage:\n");
    printf("  %s n [--dual-bin PATH] [--kl-data PATH] [--out PATH]\n", prog);
    printf("\n");
    printf("Checks dKH_w * kH_x = dKH_w * kH_y for same-left-cell pairs x,y and involutions w.\n");
    printf("Only reports triples where the product is NONZERO.\n");
    printf("Default dual cache: dual_kl_bulk_n<n>.bin\n");
    printf("Default KL data: S<n>.txt\n");
    printf("Default output: kl_triple_matches_nonzero_n<n>.csv\n");
    printf("Optional: --bench prints phase timings and extra progress info.\n");
}

static int ParseArgs(
    int argc,
    char *argv[],
    int *nOut,
    const char **dualPathOut,
    const char **klDataPathOut,
    const char **outputPathOut,
    int *benchOut
) {
    if (argc < 2) {
        return 0;
    }

    *nOut = atoi(argv[1]);
    *dualPathOut = NULL;
    *klDataPathOut = NULL;
    *outputPathOut = NULL;
    if (benchOut) *benchOut = 0;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--dual-bin") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--dual-bin requires a path.\n");
                return 0;
            }
            *dualPathOut = argv[++i];
            continue;
        }

        if (strcmp(argv[i], "--kl-data") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--kl-data requires a path.\n");
                return 0;
            }
            *klDataPathOut = argv[++i];
            continue;
        }

        if (strcmp(argv[i], "--out") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--out requires a path.\n");
                return 0;
            }
            *outputPathOut = argv[++i];
            continue;
        }

        if (strcmp(argv[i], "--bench") == 0) {
            if (benchOut) *benchOut = 1;
            continue;
        }

        fprintf(stderr, "Unknown option: %s\n", argv[i]);
        return 0;
    }

    return 1;
}

static int WriteTriple(FILE *out, int w, int x, int y, const char *permTable, int stride) {
    const char *pw = permTable + (size_t)w * (size_t)stride;
    const char *px = permTable + (size_t)x * (size_t)stride;
    const char *py = permTable + (size_t)y * (size_t)stride;
    return fprintf(out, "%d,%d,%d,%s,%s,%s\n", w, x, y, pw, px, py) > 0;
}

static int CompareBucketProducts(
    FILE *out,
    const FastContext *ctx,
    int w,
    const SparseHecke *dw,
    const TableauEntry *bucketEntries,
    int bucketSize,
    const char *permTable,
    int permStride,
    DenseAccum *accum,
    LeftMulCache *leftCache,
    SparseHecke *tmp1,
    SparseHecke *tmp2,
    SparseHecke *kx,
    SparseHecke *product
, CompareBenchStats *bench) {
    CompactHecke *products = (CompactHecke *)calloc((size_t)bucketSize, sizeof(CompactHecke));
    if (!products) {
        return 0;
    }

    for (int i = 0; i < bucketSize; i++) {
        int x = bucketEntries[i].index;

        if (bench) {
            clock_t _s = clock();
            SparseClear(kx);
            BuildKLElementSparse(ctx, x, kx);
            clock_t _e = clock();
            bench->buildKLElemCount++;
            bench->buildKLElemTime += (double)(_e - _s) / (double)CLOCKS_PER_SEC;
        } else {
            SparseClear(kx);
            BuildKLElementSparse(ctx, x, kx);
        }

        if (bench) {
            clock_t _s2 = clock();
            SparseClear(product);
            MultiplyHeckeSparse(ctx, x, dw, kx, product, tmp1, tmp2, accum, NULL, 0);
            clock_t _e2 = clock();
            bench->multiplyCount++;
            bench->multiplyTime += (double)(_e2 - _s2) / (double)CLOCKS_PER_SEC;
        } else {
            SparseClear(product);
            MultiplyHeckeSparse(ctx, x, dw, kx, product, tmp1, tmp2, accum, NULL, 0);
        }

        /* Early termination: if product is zero, mark it as such without expensive CompactFromSparse */
        if (product->supportSize <= 0) {
            memset(&products[i], 0, sizeof(CompactHecke));
            if (bench) {
                bench->zeroProductSkipCount++;
            }
            continue;
        }

        if (bench) {
            clock_t _s3 = clock();
            if (!CompactFromSparse(product, &products[i])) {
                for (int k = 0; k < i; k++) {
                    CompactFree(&products[k]);
                }
                free(products);
                return 0;
            }
            clock_t _e3 = clock();
            bench->compactFromSparseCount++;
            bench->compactFromSparseTime += (double)(_e3 - _s3) / (double)CLOCKS_PER_SEC;
        } else {
            if (!CompactFromSparse(product, &products[i])) {
                for (int k = 0; k < i; k++) {
                    CompactFree(&products[k]);
                }
                free(products);
                return 0;
            }
        }
    }

    for (int i = 0; i < bucketSize; i++) {
        /* Skip zero products */
        if (products[i].supportSize <= 0) {
            continue;
        }

        for (int j = i + 1; j < bucketSize; j++) {
            /* Skip zero products */
            if (products[j].supportSize <= 0) {
                continue;
            }

            if (bench) {
                bench->compactEqualsCount++;
            }

            if (!CompactEquals(&products[i], &products[j])) {
                continue;
            }

            int x = bucketEntries[i].index;
            int y = bucketEntries[j].index;
            if (x > y) {
                int tmp = x;
                x = y;
                y = tmp;
            }

            if (bench) {
                bench->writeTripleCount++;
            }

            if (!WriteTriple(out, w, x, y, permTable, permStride)) {
                for (int k = 0; k < bucketSize; k++) {
                    CompactFree(&products[k]);
                }
                free(products);
                return 0;
            }
        }
    }

    for (int k = 0; k < bucketSize; k++) {
        CompactFree(&products[k]);
    }
    free(products);
    (void)leftCache;
    return 1;
}

int main(int argc, char *argv[]) {
    int n = 0;
    const char *dualPath = NULL;
    const char *klDataPath = NULL;
    const char *outputPath = NULL;
    int bench = 0;
    time_t phaseStart = 0;
    time_t phaseEnd = 0;
    int64_t processedTasks = 0;
    int64_t totalTasks = 0;
    int totalBuckets = 0;
    int totalInvolutions = 0;

    if (!ParseArgs(argc, argv, &n, &dualPath, &klDataPath, &outputPath, &bench)) {
        PrintTripleUsage(argv[0]);
        return 1;
    }

    if (n < 3 || n > 8) {
        fprintf(stderr, "n must be between 3 and 8.\n");
        return 1;
    }

    char defaultDual[128];
    char defaultKlData[32];
    char defaultOut[128];
    snprintf(defaultDual, sizeof(defaultDual), "dual_kl_bulk_n%d.bin", n);
    snprintf(defaultKlData, sizeof(defaultKlData), "S%d.txt", n);
    snprintf(defaultOut, sizeof(defaultOut), "kl_triple_matches_nonzero_n%d.csv", n);

    if (!dualPath) {
        dualPath = defaultDual;
    }
    if (!klDataPath) {
        klDataPath = defaultKlData;
    }
    if (!outputPath) {
        outputPath = defaultOut;
    }

    FastContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = n;
    ctx.count = fac(n);
    ctx.maxLength = n * (n - 1) / 2;

    if (bench) fprintf(stderr, "Phase: Building context tables...\n");
    phaseStart = time(NULL);

    ctx.lengths = (int *)calloc((size_t)ctx.count, sizeof(int));
    if (!ctx.lengths) {
        fprintf(stderr, "Out of memory for length table.\n");
        return 1;
    }

    for (int i = 0; i < ctx.count; i++) {
        ctx.lengths[i] = IndexToLengthSafeLocal(n, i);
    }

    if (!BuildLengthOrder(&ctx) || !BuildGeneratorTables(&ctx) || !BuildExpressionTables(&ctx)) {
        fprintf(stderr, "Failed to build context tables.\n");
        FreeContext(&ctx);
        return 1;
    }
    phaseEnd = time(NULL);
    if (bench) fprintf(stderr, "Phase done: Building context tables (%.0f s)\n", difftime(phaseEnd, phaseStart));

    if (n >= 4) {
        if (bench) {
            fprintf(stderr, "Phase: Loading KL nontrivial data from %s...\n", klDataPath);
            phaseStart = time(NULL);
        }
        if (!LoadNonTrivialData(&ctx, klDataPath)) {
            FreeContext(&ctx);
            return 1;
        }
        if (bench) {
            phaseEnd = time(NULL);
            fprintf(stderr, "Phase done: Loading KL data (%.0f s)\n", difftime(phaseEnd, phaseStart));
        }
    }

    if (!MaybeBuildKLBasisCache(&ctx)) {
        fprintf(stderr, "KL cache build failed; continuing without cache.\n");
    }

    if (bench) fprintf(stderr, "Phase: Loading dual cache from %s...\n", dualPath);
    phaseStart = time(NULL);

    DualCache dualCache;
    memset(&dualCache, 0, sizeof(dualCache));
    int dualN = 0;
    if (!LoadDualCache(dualPath, ctx.count, &dualCache, &dualN)) {
        FreeContext(&ctx);
        return 1;
    }
    if (dualN != n) {
        fprintf(stderr, "Dual cache n=%d does not match requested n=%d.\n", dualN, n);
        DualCacheFree(&dualCache);
        FreeContext(&ctx);
        return 1;
    }
    phaseEnd = time(NULL);
    if (bench) fprintf(stderr, "Phase done: Loading dual cache (%.0f s)\n", difftime(phaseEnd, phaseStart));

    TableauEntry *entries = NULL;
    int64_t pairCount = 0;
    if (bench) fprintf(stderr, "Phase: Building same-left-cell buckets...\n");
    phaseStart = time(NULL);
    if (!BuildSameCellEntries(n, ctx.count, &entries, &pairCount)) {
        fprintf(stderr, "Failed to build same-left-cell pair buckets.\n");
        DualCacheFree(&dualCache);
        FreeContext(&ctx);
        return 1;
    }
    phaseEnd = time(NULL);
    if (bench) fprintf(stderr, "Phase done: Building buckets (%.0f s)\n", difftime(phaseEnd, phaseStart));

    FILE *out = fopen(outputPath, "w");
    if (!out) {
        fprintf(stderr, "Could not open output file %s\n", outputPath);
        free(entries);
        DualCacheFree(&dualCache);
        FreeContext(&ctx);
        return 1;
    }

    fprintf(out, "w,x,y,w_perm,x_perm,y_perm\n");

    char *permTable = BuildPermStringTable(n, ctx.count);
    if (!permTable) {
        fprintf(stderr, "Failed to build permutation string table.\n");
        fclose(out);
        free(entries);
        DualCacheFree(&dualCache);
        FreeContext(&ctx);
        return 1;
    }

    SparseHecke kxTemp, tmp1, tmp2, product;
    DenseAccum accum;
    LeftMulCache leftCache;
    if (!SparseInit(&kxTemp, ctx.count) || !SparseInit(&tmp1, ctx.count) ||
        !SparseInit(&tmp2, ctx.count) || !SparseInit(&product, ctx.count) || !DenseAccumInit(&accum, ctx.count) ||
        !LeftMulCacheInit(&leftCache, ctx.count)) {
        fprintf(stderr, "Failed to allocate multiplication scratch space.\n");
        free(permTable);
        fclose(out);
        free(entries);
        DualCacheFree(&dualCache);
        FreeContext(&ctx);
        return 1;
    }

    CompareBenchStats benchStats;
    memset(&benchStats, 0, sizeof(benchStats));

    /* count total involutions present in the dual cache for progress reporting */
    for (int ii = 0; ii < ctx.count; ii++) {
        if (dualCache.present[ii]) totalInvolutions++;
    }

    /* compute total tasks = sum over buckets (if bucketSize>1) of totalInvolutions */
    for (int ii = 0; ii < ctx.count; ) {
        int jj = ii + 1;
        while (jj < ctx.count && memcmp(entries[ii].key.bytes, entries[jj].key.bytes, sizeof(entries[ii].key.bytes)) == 0) jj++;
        int bsize = jj - ii;
        if (bsize > 1) {
            totalBuckets++;
            totalTasks += (int64_t)totalInvolutions;
        }
        ii = jj;
    }

    if (bench) {
        fprintf(stderr, "Phase: Comparing buckets and involutions... total tasks=%lld (buckets=%d, involutions=%d)\n", (long long)totalTasks, totalBuckets, totalInvolutions);
        phaseStart = time(NULL);
    }

    for (int i = 0; i < ctx.count; ) {
        int j = i + 1;
        while (j < ctx.count && memcmp(entries[i].key.bytes, entries[j].key.bytes, sizeof(entries[i].key.bytes)) == 0) {
            j++;
        }

        int bucketSize = j - i;
        if (bucketSize > 1) {
            for (int w = 0; w < ctx.count; w++) {
                if (!dualCache.present[w]) {
                    continue;
                }

                processedTasks++;
                if ((processedTasks % 5) == 0) {
                    fprintf(stderr, "Progress: processed %lld/%lld tasks (last w=%d)\n", (long long)processedTasks, (long long)totalTasks, w);
                }

                if (!CompareBucketProducts(out, &ctx, w, &dualCache.records[w], &entries[i], bucketSize, permTable, n + 1, &accum, &leftCache, &tmp1, &tmp2, &kxTemp, &product, bench ? &benchStats : NULL)) {
                    fprintf(stderr, "Failed while comparing bucket starting at %d for w=%d.\n", i, w);
                    SparseFree(&kxTemp);
                    SparseFree(&tmp1);
                    SparseFree(&tmp2);
                    SparseFree(&product);
                    DenseAccumFree(&accum);
                    LeftMulCacheFree(&leftCache);
                    free(permTable);
                    fclose(out);
                    free(entries);
                    DualCacheFree(&dualCache);
                    FreeContext(&ctx);
                    return 1;
                }
            }
        }

        i = j;
    }

    fclose(out);
    if (bench) {
        phaseEnd = time(NULL);
        fprintf(stderr, "Phase done: Comparing buckets (%.0f s). Total tasks processed: %lld\n", difftime(phaseEnd, phaseStart), (long long)processedTasks);
    }
    if (bench) {
        fprintf(stderr, "Bench summary for compare phase:\n");
        fprintf(stderr, "  BuildKLElem: count=%ld time=%.6f s\n", benchStats.buildKLElemCount, benchStats.buildKLElemTime);
        fprintf(stderr, "  Multiply:    count=%ld time=%.6f s\n", benchStats.multiplyCount, benchStats.multiplyTime);
        fprintf(stderr, "  ZeroSkip:    count=%ld (%.1f%% of products skipped CompactFromSparse)\n", 
                benchStats.zeroProductSkipCount, 
                benchStats.multiplyCount > 0 ? (100.0 * benchStats.zeroProductSkipCount / benchStats.multiplyCount) : 0.0);
        fprintf(stderr, "  CompactFrom: count=%ld time=%.6f s\n", benchStats.compactFromSparseCount, benchStats.compactFromSparseTime);
        fprintf(stderr, "  CompactEq:   count=%ld time=%.6f s\n", benchStats.compactEqualsCount, benchStats.compactEqualsTime);
        fprintf(stderr, "  WriteTriple: count=%ld time=%.6f s\n", benchStats.writeTripleCount, benchStats.writeTripleTime);
    }
    printf("Processed %lld same-left-cell pairs and wrote matching nonzero triples to %s\n", (long long)pairCount, outputPath);

    free(permTable);
    free(entries);
    SparseFree(&kxTemp);
    SparseFree(&tmp1);
    SparseFree(&tmp2);
    SparseFree(&product);
    DenseAccumFree(&accum);
    LeftMulCacheFree(&leftCache);
    DualCacheFree(&dualCache);
    FreeContext(&ctx);
    return 0;
}
