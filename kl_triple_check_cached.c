#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

#include <omp.h>

typedef struct {
    unsigned char bytes[72];
} TableauKey;

typedef struct {
    int index;
    TableauKey key;
} TableauEntry;

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
    uint64_t build_count;
    uint64_t multiply_count;
    uint64_t multiply_miss_count;
    uint64_t compact_from_count;
    uint64_t compact_eq_count;
    uint64_t write_count;
    double build_time;
    double multiply_time;
    double leftmul_time;
    double mulPolyWall;
    double mulAccumWall;
    double mulCommitWall;
    double compact_from_time;
    double compact_eq_time;
    double write_time;
} CompareBenchStats;

typedef struct {
    SparseHecke *records;
    unsigned char *present;
    int count;
} DualCache;

static int IndexToLengthSafeLocal(int n, int index) {
    return IndexToLength(n, index);
}

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
    const int ia = ta->index;
    const int ib = tb->index;
    return ia - ib;
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
        fprintf(stderr, "Dual KL cache count=%u does not match n=%d.\n", header.count, count);
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

static void PrintOptimizedUsage(const char *prog) {
    printf("Usage:\n");
    printf("  %s n [--dual-bin PATH] [--kl-data PATH] [--out PATH] [--bench] [--materialize-top K]\n", prog);
    printf("\n");
    printf("Optimized checker with cached kH_x and OpenMP parallelism across involutions.\n");
    printf("Default dual cache: dual_kl_bulk_n<n>.bin\n");
    printf("Default KL data: S<n>.txt\n");
    printf("Default output: kl_triple_matches_cached_n<n>.csv\n");
}

static int ParseArgsCached(
    int argc,
    char *argv[],
    int *nOut,
    const char **dualPathOut,
    const char **klDataOut,
    const char **outputPathOut,
    int *benchOut,
    int *materializeTopOut
) {
    if (argc < 2) {
        return 0;
    }

    *nOut = atoi(argv[1]);
    *dualPathOut = NULL;
    *klDataOut = NULL;
    *outputPathOut = NULL;
    if (benchOut) *benchOut = 0;
    if (materializeTopOut) *materializeTopOut = 0;

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
            *klDataOut = argv[++i];
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

        if (strcmp(argv[i], "--materialize-top") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--materialize-top requires a number K.\n");
                return 0;
            }
            if (materializeTopOut) *materializeTopOut = atoi(argv[++i]);
            continue;
        }

        fprintf(stderr, "Unknown option: %s\n", argv[i]);
        return 0;
    }

    return 1;
}

static int WriteTripleSafe(FILE *out, int w, int x, int y, const char *permTable, int stride) {
    const char *pw = permTable + (size_t)w * (size_t)stride;
    const char *px = permTable + (size_t)x * (size_t)stride;
    const char *py = permTable + (size_t)y * (size_t)stride;
    return fprintf(out, "%d,%d,%d,%s,%s,%s\n", w, x, y, pw, px, py) > 0;
}

static int CollectInvolutionIndices(const DualCache *dualCache, int count, int **indicesOut, int *indicesCountOut) {
    int *indices = (int *)calloc((size_t)count, sizeof(int));
    if (!indices) {
        return 0;
    }

    int indicesCount = 0;
    for (int i = 0; i < count; i++) {
        if (dualCache->present[i]) {
            indices[indicesCount++] = i;
        }
    }

    *indicesOut = indices;
    *indicesCountOut = indicesCount;
    return 1;
}

static void FreeSparseArray(SparseHecke *array, int count) {
    if (!array) {
        return;
    }

    for (int i = 0; i < count; i++) {
        SparseFree(&array[i]);
    }
    free(array);
}

static int MultiplyUsingLeftCache(
    const FastContext *ctx,
    int rightId,
    const SparseHecke *left,
    const SparseHecke *right,
    SparseHecke *product,
    SparseHecke *tmp1,
    SparseHecke *tmp2,
    DenseAccum *accum,
    LeftMulCache *leftCache,
    CompareBenchStats *stats,
    int *leftFreq,
    SparseHecke *materialized, /* flattened [hotCount * bucketSize] or NULL */
    const int *hotMap, /* maps leftIndex -> hotIdx or -1 */
    int bucketSize,
    int rightPos,
    int histogramOnly
) {
    /* Emulate MultiplyHeckeSparse but use LeftMulGetCached to reuse left*right rows */
    DenseAccumBegin(accum);

    for (int lp = 0; lp < left->supportSize; lp++) {
        int i = left->support[lp];
        if (leftFreq) leftFreq[i]++;
        if (histogramOnly) continue;
        const Poly57 *scalar = &left->coeff[i];
        int cacheMiss = 0;
        double lt0 = 0.0;
        if (stats) lt0 = omp_get_wtime();
        /* Check materialized hot rows first */
        const SparseHecke *lm = NULL;
        if (hotMap && hotMap[i] >= 0 && materialized) {
            int h = hotMap[i];
            lm = &materialized[(size_t)h * (size_t)bucketSize + (size_t)rightPos];
            cacheMiss = 0;
        } else {
            lm = LeftMulGetCached(ctx, rightId, i, right, tmp1, tmp2, leftCache, &cacheMiss);
        }
        if (stats) {
            stats->multiply_count++;
            if (cacheMiss) stats->multiply_miss_count++;
            stats->leftmul_time += omp_get_wtime() - lt0;
        }

        for (int rp = 0; rp < lm->supportSize; rp++) {
            int j = lm->support[rp];
            Poly57 scaled;
            double t0 = 0.0;
            if (stats) t0 = omp_get_wtime();
            PolyMultiply(scalar, &lm->coeff[j], &scaled);
            if (stats) stats->mulPolyWall += omp_get_wtime() - t0;

            if (stats) t0 = omp_get_wtime();
            DenseAccumAdd(accum, j, &scaled);
            if (stats) stats->mulAccumWall += omp_get_wtime() - t0;
        }
    }

    if (stats) {
        double t0 = omp_get_wtime();
        DenseAccumCommitToSparse(accum, product);
        stats->mulCommitWall += omp_get_wtime() - t0;
    } else {
        DenseAccumCommitToSparse(accum, product);
    }
    return 1;
}

static int CompareBucketForInvolution(
    FILE *out,
    const FastContext *ctx,
    int w,
    const SparseHecke *dw,
    const TableauEntry *bucketEntries,
    int bucketSize,
    const SparseHecke *cachedKx,
    const char *permTable,
    int permStride,
    SparseHecke *tmp1,
    SparseHecke *tmp2,
    SparseHecke *product,
    DenseAccum *accum,
    CompareBenchStats *stats,
    LeftMulCache *leftCache,
    int *leftFreq,
    SparseHecke *materialized,
    const int *hotMap,
    int hotCount,
    int histogramOnly
) {
    CompactHecke *products = (CompactHecke *)calloc((size_t)bucketSize, sizeof(CompactHecke));
    if (!products) {
        return 0;
    }

        for (int i = 0; i < bucketSize; i++) {
        SparseClear(product);
        if (leftCache) {
            if (!MultiplyUsingLeftCache(ctx, bucketEntries[i].index, dw, &cachedKx[i], product, tmp1, tmp2, accum, leftCache, stats, leftFreq, materialized, hotMap, bucketSize, i, histogramOnly)) {
                for (int k = 0; k < i; k++) CompactFree(&products[k]);
                free(products);
                return 0;
            }
            if (histogramOnly) continue;
            if (stats) {
                double t1 = omp_get_wtime();
                if (!CompactFromSparse(product, &products[i])) {
                    for (int k = 0; k < i; k++) CompactFree(&products[k]);
                    free(products);
                    return 0;
                }
                stats->compact_from_count++;
                stats->compact_from_time += omp_get_wtime() - t1;
            } else {
                if (!CompactFromSparse(product, &products[i])) {
                    for (int k = 0; k < i; k++) CompactFree(&products[k]);
                    free(products);
                    return 0;
                }
            }
        } else {
            /* fallback to existing MultiplyHeckeSparse */
            if (histogramOnly) continue;
            if (stats) {
                double t0 = omp_get_wtime();
                MultiplyHeckeSparse(ctx, bucketEntries[i].index, dw, &cachedKx[i], product, tmp1, tmp2, accum, NULL, 0);
                stats->multiply_count++;
                stats->multiply_time += omp_get_wtime() - t0;

                double t1 = omp_get_wtime();
                if (!CompactFromSparse(product, &products[i])) {
                    for (int k = 0; k < i; k++) CompactFree(&products[k]);
                    free(products);
                    return 0;
                }
                stats->compact_from_count++;
                stats->compact_from_time += omp_get_wtime() - t1;
            } else {
                MultiplyHeckeSparse(ctx, bucketEntries[i].index, dw, &cachedKx[i], product, tmp1, tmp2, accum, NULL, 0);
                if (!CompactFromSparse(product, &products[i])) {
                    for (int k = 0; k < i; k++) CompactFree(&products[k]);
                    free(products);
                    return 0;
                }
            }
        }
    }

    int ok = 1;
    for (int i = 0; i < bucketSize; i++) {
        for (int j = i + 1; j < bucketSize; j++) {
            if (stats) {
                double t0 = omp_get_wtime();
                int eq = CompactEquals(&products[i], &products[j]);
                stats->compact_eq_count++;
                stats->compact_eq_time += omp_get_wtime() - t0;
                if (!eq) continue;
            } else {
                if (!CompactEquals(&products[i], &products[j])) continue;
            }

            int x = bucketEntries[i].index;
            int y = bucketEntries[j].index;
            if (x > y) {
                int tmp = x;
                x = y;
                y = tmp;
            }

            #pragma omp critical(triple_output)
            {
                if (stats) {
                    double t1 = omp_get_wtime();
                    if (!WriteTripleSafe(out, w, x, y, permTable, permStride)) {
                        ok = 0;
                    }
                    stats->write_count++;
                    stats->write_time += omp_get_wtime() - t1;
                } else {
                    if (!WriteTripleSafe(out, w, x, y, permTable, permStride)) {
                        ok = 0;
                    }
                }
            }
        }
    }

    for (int k = 0; k < bucketSize; k++) {
        CompactFree(&products[k]);
    }
    free(products);
    return ok;
}

static int ProcessBucketParallel(
    FILE *out,
    const FastContext *ctx,
    const DualCache *dualCache,
    const int *involutionIndices,
    int involutionCount,
    const TableauEntry *bucketEntries,
    int bucketSize,
    const char *permTable,
    int permStride,
    int bench,
    CompareBenchStats *benchStats,
    int64_t *processedTasks,
    int64_t totalTasks,
    int *globalLeftFreq,
    const int *hotIndices,
    int hotCount,
    int histogramOnly
) {
    SparseHecke *cachedKx = (SparseHecke *)calloc((size_t)bucketSize, sizeof(SparseHecke));
    if (!cachedKx) {
        return 0;
    }

    for (int i = 0; i < bucketSize; i++) {
        if (!SparseInit(&cachedKx[i], ctx->count)) {
            FreeSparseArray(cachedKx, i);
            return 0;
        }
        if (bench && benchStats) {
            double t0 = omp_get_wtime();
            BuildKLElementSparse(ctx, bucketEntries[i].index, &cachedKx[i]);
            benchStats->build_count++;
            benchStats->build_time += omp_get_wtime() - t0;
        } else {
            BuildKLElementSparse(ctx, bucketEntries[i].index, &cachedKx[i]);
        }
    }

    /* If requested, materialize hot-left rows for this bucket */
    SparseHecke *materialized = NULL;
    int *hotMap = NULL;
    if (!histogramOnly && hotIndices && hotCount > 0) {
        materialized = (SparseHecke *)calloc((size_t)hotCount * (size_t)bucketSize, sizeof(SparseHecke));
        hotMap = (int *)malloc((size_t)ctx->count * sizeof(int));
        if (!materialized || !hotMap) {
            if (materialized) free(materialized);
            if (hotMap) free(hotMap);
            materialized = NULL; hotMap = NULL;
            /* continue without materialization */
        } else {
            for (int ii = 0; ii < ctx->count; ii++) hotMap[ii] = -1;
            for (int h = 0; h < hotCount; h++) {
                int li = hotIndices[h];
                if (li >= 0 && li < ctx->count) hotMap[li] = h;
            }
            /* initialize and compute materialized entries */
            /* prepare a scratch for LeftMultiplyByIndex */
            SparseHecke scratch;
            if (!SparseInit(&scratch, ctx->count)) {
                free(materialized);
                free(hotMap);
                materialized = NULL; hotMap = NULL;
            } else {
                for (int h = 0; h < hotCount; h++) {
                    for (int j = 0; j < bucketSize; j++) {
                        SparseHecke *dest = &materialized[(size_t)h * (size_t)bucketSize + (size_t)j];
                        if (!SparseInit(dest, ctx->count)) {
                        /* fall back: free allocated materialized */
                            for (int hh = 0; hh < h; hh++) {
                                for (int jj = 0; jj < bucketSize; jj++) SparseFree(&materialized[(size_t)hh * (size_t)bucketSize + (size_t)jj]);
                            }
                            SparseFree(&scratch);
                            free(materialized);
                            free(hotMap);
                            materialized = NULL; hotMap = NULL;
                            break;
                    }
                        LeftMultiplyByIndex(ctx, hotIndices[h], &cachedKx[j], dest, &scratch);
                }
                    if (!materialized) break;
                }
                if (materialized) SparseFree(&scratch);
            }
        }
    }

    int ok = 1;

    #pragma omp parallel
    {
        SparseHecke tmp1;
        SparseHecke tmp2;
        SparseHecke product;
        DenseAccum accum;
        LeftMulCache leftCacheThread;
        int *leftFreq = NULL;
        int localOk = 1;
        CompareBenchStats threadStats;
        memset(&threadStats, 0, sizeof(threadStats));

        if (!SparseInit(&tmp1, ctx->count) || !SparseInit(&tmp2, ctx->count) ||
            !SparseInit(&product, ctx->count) || !DenseAccumInit(&accum, ctx->count)) {
            localOk = 0;
        }
        memset(&leftCacheThread, 0, sizeof(leftCacheThread));
        if (!LeftMulCacheInit(&leftCacheThread, ctx->count)) {
            localOk = 0;
        }
        if (bench && globalLeftFreq) {
            leftFreq = (int *)calloc((size_t)ctx->count, sizeof(int));
            if (!leftFreq) leftFreq = NULL;
        }

        #pragma omp for schedule(dynamic)
        for (int wi = 0; wi < involutionCount; wi++) {
            if (!localOk || !ok) {
                continue;
            }

            int w = involutionIndices[wi];
            if (!CompareBucketForInvolution(out, ctx, w, &dualCache->records[w], bucketEntries, bucketSize, cachedKx, permTable, permStride, &tmp1, &tmp2, &product, &accum, (bench && benchStats) ? &threadStats : NULL, &leftCacheThread, leftFreq, materialized, hotMap, hotCount, histogramOnly)) {
                ok = 0;
            }

            int64_t done = 0;
            /* Count each (w,x) product as a task: add bucketSize per processed involution */
            #pragma omp atomic capture
            done = (*processedTasks) += bucketSize;

            if (bench && (done % 5) == 0) {
                #pragma omp critical(progress_output)
                fprintf(stderr, "Progress: processed %lld/%lld tasks (last w=%d)\n", (long long)done, (long long)totalTasks, w);
            }
        }

        /* Aggregate thread-local stats into shared benchStats */
        if (bench && benchStats) {
            #pragma omp critical(bench_update)
            {
                benchStats->build_count += threadStats.build_count;
                benchStats->multiply_count += threadStats.multiply_count;
                benchStats->multiply_miss_count += threadStats.multiply_miss_count;
                benchStats->compact_from_count += threadStats.compact_from_count;
                benchStats->compact_eq_count += threadStats.compact_eq_count;
                benchStats->write_count += threadStats.write_count;
                benchStats->build_time += threadStats.build_time;
                benchStats->multiply_time += threadStats.multiply_time;
                benchStats->leftmul_time += threadStats.leftmul_time;
                benchStats->mulPolyWall += threadStats.mulPolyWall;
                benchStats->mulAccumWall += threadStats.mulAccumWall;
                benchStats->mulCommitWall += threadStats.mulCommitWall;
                benchStats->compact_from_time += threadStats.compact_from_time;
                benchStats->compact_eq_time += threadStats.compact_eq_time;
                benchStats->write_time += threadStats.write_time;
            }
        }

        LeftMulCacheFree(&leftCacheThread);
        if (leftFreq) {
            /* aggregate per-thread counts into global histogram */
            #pragma omp critical(hist_update)
            {
                for (int ii = 0; ii < ctx->count; ii++) {
                    if (leftFreq[ii]) globalLeftFreq[ii] += leftFreq[ii];
                }
            }
            free(leftFreq);
        }
        SparseFree(&tmp1);
        SparseFree(&tmp2);
        SparseFree(&product);
        DenseAccumFree(&accum);
    }
        FreeSparseArray(cachedKx, bucketSize);

        if (materialized) {
            for (int h = 0; h < hotCount; h++) {
                for (int j = 0; j < bucketSize; j++) {
                    SparseFree(&materialized[(size_t)h * (size_t)bucketSize + (size_t)j]);
                }
            }
            free(materialized);
        }
        if (hotMap) free(hotMap);
    return ok;
}

static int RunCachedChecker(int argc, char *argv[]) {
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
    int totalCachedElements = 0;

    int materializeTop = 0;
    if (!ParseArgsCached(argc, argv, &n, &dualPath, &klDataPath, &outputPath, &bench, &materializeTop)) {
        PrintOptimizedUsage(argv[0]);
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
    snprintf(defaultOut, sizeof(defaultOut), "kl_triple_matches_cached_n%d.csv", n);

    if (!dualPath) dualPath = defaultDual;
    if (!klDataPath) klDataPath = defaultKlData;
    if (!outputPath) outputPath = defaultOut;

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

    if (bench) {
        fprintf(stderr, "Phase: Loading dual cache from %s...\n", dualPath);
        phaseStart = time(NULL);
    }

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

    if (bench) {
        phaseEnd = time(NULL);
        fprintf(stderr, "Phase done: Loading dual cache (%.0f s)\n", difftime(phaseEnd, phaseStart));
    }

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

    int *involutionIndices = NULL;
    int involutionCount = 0;
    if (!CollectInvolutionIndices(&dualCache, ctx.count, &involutionIndices, &involutionCount)) {
        fprintf(stderr, "Failed to collect involution indices.\n");
        free(entries);
        DualCacheFree(&dualCache);
        FreeContext(&ctx);
        return 1;
    }

    FILE *out = fopen(outputPath, "w");
    if (!out) {
        fprintf(stderr, "Could not open output file %s\n", outputPath);
        free(involutionIndices);
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
        free(involutionIndices);
        free(entries);
        DualCacheFree(&dualCache);
        FreeContext(&ctx);
        return 1;
    }

    for (int ii = 0; ii < ctx.count; ii++) {
        if (dualCache.present[ii]) {
            totalInvolutions++;
        }
    }

    for (int ii = 0; ii < ctx.count; ) {
        int jj = ii + 1;
        while (jj < ctx.count && memcmp(entries[ii].key.bytes, entries[jj].key.bytes, sizeof(entries[ii].key.bytes)) == 0) {
            jj++;
        }
        int bucketSize = jj - ii;
        if (bucketSize > 1) {
            totalBuckets++;
            totalCachedElements += bucketSize;
            totalTasks += (int64_t)bucketSize * (int64_t)involutionCount;
        }
        ii = jj;
    }

    if (bench) {
        fprintf(stderr, "Phase: Comparing buckets and involutions... total tasks=%lld (buckets=%d, involutions=%d, cached x=%d)\n", (long long)totalTasks, totalBuckets, involutionCount, totalCachedElements);
        phaseStart = time(NULL);
    }

    int ok = 1;
    CompareBenchStats compareBench;
    memset(&compareBench, 0, sizeof(compareBench));
    int *globalLeftFreq = NULL;
    if (bench) {
        globalLeftFreq = (int *)calloc((size_t)ctx.count, sizeof(int));
    }

    if (materializeTop > 0 && bench && globalLeftFreq) {
        /* First pass: histogram-only to gather left-index frequencies */
        for (int i = 0; i < ctx.count; ) {
            int j = i + 1;
            while (j < ctx.count && memcmp(entries[i].key.bytes, entries[j].key.bytes, sizeof(entries[i].key.bytes)) == 0) {
                j++;
            }

            int bucketSize = j - i;
            if (bucketSize > 1) {
                if (!ProcessBucketParallel(out, &ctx, &dualCache, involutionIndices, involutionCount, &entries[i], bucketSize, permTable, n + 1, 1, NULL, &processedTasks, totalTasks, globalLeftFreq, NULL, 0, 1)) {
                    ok = 0;
                    break;
                }
            }

            i = j;
        }

        if (ok) {
            /* choose top-K hot indices */
            int K = materializeTop;
            if (K > ctx.count) K = ctx.count;
            int *idxs = (int *)calloc((size_t)ctx.count, sizeof(int));
            int *topIndices = NULL;
            if (idxs) {
                for (int ii = 0; ii < ctx.count; ii++) idxs[ii] = ii;
                for (int a = 0; a < K; a++) {
                    int best = a;
                    for (int b = a+1; b < ctx.count; b++) {
                        if (globalLeftFreq[idxs[b]] > globalLeftFreq[idxs[best]]) best = b;
                    }
                    int tmp = idxs[a]; idxs[a] = idxs[best]; idxs[best] = tmp;
                }
                topIndices = (int *)malloc((size_t)K * sizeof(int));
                for (int t = 0; t < K; t++) topIndices[t] = idxs[t];
                free(idxs);
            }

            /* Reset progress counters and run a second pass that materializes top rows */
            processedTasks = 0;
            for (int i = 0; i < ctx.count; ) {
                int j = i + 1;
                while (j < ctx.count && memcmp(entries[i].key.bytes, entries[j].key.bytes, sizeof(entries[i].key.bytes)) == 0) {
                    j++;
                }

                int bucketSize = j - i;
                if (bucketSize > 1) {
                    if (!ProcessBucketParallel(out, &ctx, &dualCache, involutionIndices, involutionCount, &entries[i], bucketSize, permTable, n + 1, bench, (bench ? &compareBench : NULL), &processedTasks, totalTasks, NULL, topIndices, K, 0)) {
                        ok = 0;
                        break;
                    }
                }

                i = j;
            }

            if (topIndices) free(topIndices);
        }
    } else {
        for (int i = 0; i < ctx.count; ) {
            int j = i + 1;
            while (j < ctx.count && memcmp(entries[i].key.bytes, entries[j].key.bytes, sizeof(entries[i].key.bytes)) == 0) {
                j++;
            }

            int bucketSize = j - i;
                if (bucketSize > 1) {
                if (!ProcessBucketParallel(out, &ctx, &dualCache, involutionIndices, involutionCount, &entries[i], bucketSize, permTable, n + 1, bench, (bench ? &compareBench : NULL), &processedTasks, totalTasks, globalLeftFreq, NULL, 0, 0)) {
                    ok = 0;
                    break;
                }
            }

            i = j;
        }
    }

    fclose(out);
    if (bench) {
        phaseEnd = time(NULL);
        fprintf(stderr, "Phase done: Comparing buckets (%.0f s). Total tasks processed: %lld\n", difftime(phaseEnd, phaseStart), (long long)processedTasks);
    }

    if (bench && globalLeftFreq) {
        /* print top-20 left indices by frequency */
        int topN = 20;
        int *idxs = (int *)calloc((size_t)ctx.count, sizeof(int));
        if (idxs) {
            for (int ii = 0; ii < ctx.count; ii++) idxs[ii] = ii;
            /* simple partial sort by frequency descending */
            for (int a = 0; a < topN && a < ctx.count; a++) {
                int best = a;
                for (int b = a+1; b < ctx.count; b++) {
                    if (globalLeftFreq[idxs[b]] > globalLeftFreq[idxs[best]]) best = b;
                }
                int tmp = idxs[a]; idxs[a] = idxs[best]; idxs[best] = tmp;
            }

            fprintf(stderr, "\nTop-%d left-index frequencies:\n", topN);
            for (int k = 0; k < topN && k < ctx.count; k++) {
                int ii = idxs[k];
                if (globalLeftFreq[ii] == 0) break;
                fprintf(stderr, "  leftIndex=%d freq=%d\n", ii, globalLeftFreq[ii]);
            }
            free(idxs);
        }
    }

    if (bench) {
        fprintf(stderr, "\nCompare phase bench summary:\n");
        fprintf(stderr, "  BuildKLElem: count=%llu time=%.6f s\n", (unsigned long long)compareBench.build_count, compareBench.build_time);
        fprintf(stderr, "  LeftMulGetCached: count=%llu (misses=%llu) time=%.6f s\n", (unsigned long long)compareBench.multiply_count, (unsigned long long)compareBench.multiply_miss_count, compareBench.leftmul_time);
        fprintf(stderr, "  Multiply internals (wall): poly=%.6f s, accum=%.6f s, commit=%.6f s\n", compareBench.mulPolyWall, compareBench.mulAccumWall, compareBench.mulCommitWall);
        fprintf(stderr, "  CompactFrom: count=%llu time=%.6f s\n", (unsigned long long)compareBench.compact_from_count, compareBench.compact_from_time);
        fprintf(stderr, "  CompactEq:   count=%llu time=%.6f s\n", (unsigned long long)compareBench.compact_eq_count, compareBench.compact_eq_time);
        fprintf(stderr, "  WriteTriple: count=%llu time=%.6f s\n", (unsigned long long)compareBench.write_count, compareBench.write_time);
    }

    printf("Processed %lld same-left-cell pairs and wrote matching triples to %s\n", (long long)pairCount, outputPath);

    free(permTable);
    free(involutionIndices);
    free(entries);
    DualCacheFree(&dualCache);
    FreeContext(&ctx);
    return ok ? 0 : 1;
}

int main(int argc, char *argv[]) {
    return RunCachedChecker(argc, argv);
}
