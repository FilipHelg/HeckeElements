#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Reuse the optimized parallel implementation internals in this translation unit. */
#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

static void EstimatePrintUsage(const char *prog) {
    printf("Usage:\n");
    printf("  %s n [--sample-rows K] [--threads N] [--row-schedule block|weighted] [--baseline-seconds T] [--baseline-n M]\n", prog);
    printf("\n");
    printf("n must be in {4,5,6,7,8,9}.\n");
    printf("Samples K x-rows (full y-loop for each sampled row), then extrapolates full S<n> loop time.\n");
    printf("--row-schedule weighted is default and samples by row-weight quantiles.\n");
    printf("If --baseline-seconds and --baseline-n are provided, also fit a in T ~ (n!)^a and predict S<n+1>, S<n+2>.\n");
}

static int BuildUniformSampleRows(int count, int sampleCount, int *rows) {
    if (!rows || sampleCount <= 0 || count <= 0) {
        return 0;
    }

    if (sampleCount >= count) {
        for (int i = 0; i < count; i++) {
            rows[i] = i;
        }
        return count;
    }

    int usedCount = 0;
    char *used = (char *)calloc((size_t)count, sizeof(char));
    if (!used) {
        return 0;
    }

    for (int i = 0; i < sampleCount; i++) {
        int x = (sampleCount == 1) ? (count / 2) : (int)(((int64_t)i * (int64_t)(count - 1)) / (int64_t)(sampleCount - 1));
        if (!used[x]) {
            used[x] = 1;
            rows[usedCount++] = x;
        }
    }

    for (int x = 0; x < count && usedCount < sampleCount; x++) {
        if (!used[x]) {
            used[x] = 1;
            rows[usedCount++] = x;
        }
    }

    free(used);
    return usedCount;
}

static int FindNearestUnused(const char *used, int count, int preferred) {
    if (preferred < 0) preferred = 0;
    if (preferred >= count) preferred = count - 1;
    if (!used[preferred]) {
        return preferred;
    }

    for (int d = 1; d < count; d++) {
        int left = preferred - d;
        int right = preferred + d;
        if (left >= 0 && !used[left]) {
            return left;
        }
        if (right < count && !used[right]) {
            return right;
        }
    }

    return -1;
}

static int BuildWeightedSampleRows(const int *weights, int count, int sampleCount, int *rows) {
    if (!weights || !rows || sampleCount <= 0 || count <= 0) {
        return 0;
    }

    if (sampleCount >= count) {
        for (int i = 0; i < count; i++) {
            rows[i] = i;
        }
        return count;
    }

    int64_t totalWeight = 0;
    for (int i = 0; i < count; i++) {
        int w = weights[i] > 0 ? weights[i] : 1;
        totalWeight += (int64_t)w;
    }

    if (totalWeight <= 0) {
        return BuildUniformSampleRows(count, sampleCount, rows);
    }

    int usedCount = 0;
    char *used = (char *)calloc((size_t)count, sizeof(char));
    if (!used) {
        return 0;
    }

    int x = 0;
    int64_t prefix = (count > 0) ? (int64_t)(weights[0] > 0 ? weights[0] : 1) : 0;

    for (int i = 0; i < sampleCount; i++) {
        int64_t target = ((int64_t)(2 * i + 1) * totalWeight) / (int64_t)(2 * sampleCount);
        while (x + 1 < count && prefix < target) {
            x++;
            prefix += (int64_t)(weights[x] > 0 ? weights[x] : 1);
        }

        int pick = FindNearestUnused(used, count, x);
        if (pick >= 0) {
            used[pick] = 1;
            rows[usedCount++] = pick;
        }
    }

    for (int i = 0; i < count && usedCount < sampleCount; i++) {
        if (!used[i]) {
            used[i] = 1;
            rows[usedCount++] = i;
        }
    }

    free(used);
    return usedCount;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        EstimatePrintUsage(argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    if (n < 4 || n > 9) {
        fprintf(stderr, "n must be between 4 and 9.\n");
        return 1;
    }

    int sampleRowsRequested = 64;
    int requestedThreads = 0;
    RowScheduleMode rowScheduleMode = ROW_SCHEDULE_WEIGHTED;
    double baselineSeconds = 0.0;
    int baselineN = 0;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--sample-rows") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--sample-rows requires a positive integer.\n");
                return 1;
            }
            sampleRowsRequested = atoi(argv[++i]);
            if (sampleRowsRequested <= 0) {
                fprintf(stderr, "--sample-rows requires a positive integer.\n");
                return 1;
            }
            continue;
        }

        if (strcmp(argv[i], "--threads") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--threads requires a positive integer.\n");
                return 1;
            }
            requestedThreads = atoi(argv[++i]);
            if (requestedThreads <= 0) {
                fprintf(stderr, "--threads requires a positive integer.\n");
                return 1;
            }
            continue;
        }

        if (strcmp(argv[i], "--row-schedule") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--row-schedule requires one of: block, weighted\n");
                return 1;
            }
            const char *mode = argv[++i];
            if (strcmp(mode, "block") == 0) {
                rowScheduleMode = ROW_SCHEDULE_BLOCK;
            } else if (strcmp(mode, "weighted") == 0) {
                rowScheduleMode = ROW_SCHEDULE_WEIGHTED;
            } else {
                fprintf(stderr, "Unknown --row-schedule mode '%s' (use block or weighted).\n", mode);
                return 1;
            }
            continue;
        }

        if (strcmp(argv[i], "--baseline-seconds") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--baseline-seconds requires a positive number.\n");
                return 1;
            }
            baselineSeconds = atof(argv[++i]);
            if (baselineSeconds <= 0.0) {
                fprintf(stderr, "--baseline-seconds requires a positive number.\n");
                return 1;
            }
            continue;
        }

        if (strcmp(argv[i], "--baseline-n") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--baseline-n requires an integer in [4,9].\n");
                return 1;
            }
            baselineN = atoi(argv[++i]);
            if (baselineN < 4 || baselineN > 9) {
                fprintf(stderr, "--baseline-n requires an integer in [4,9].\n");
                return 1;
            }
            continue;
        }

        fprintf(stderr, "Unknown argument: %s\n", argv[i]);
        EstimatePrintUsage(argv[0]);
        return 1;
    }

    if ((baselineSeconds > 0.0) != (baselineN > 0)) {
        fprintf(stderr, "Use --baseline-seconds and --baseline-n together.\n");
        return 1;
    }

    FastContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = n;
    ctx.count = fac(n);
    ctx.maxLength = n * (n - 1) / 2;
    ctx.useCache = 0;

    ctx.lengths = (int *)calloc((size_t)ctx.count, sizeof(int));
    if (!ctx.lengths) {
        fprintf(stderr, "Out of memory for length table.\n");
        return 1;
    }
    for (int i = 0; i < ctx.count; i++) {
        ctx.lengths[i] = IndexToLength(n, i);
    }

#ifdef _OPENMP
    int availableThreads = omp_get_max_threads();
    int desiredThreads = requestedThreads > 0 ? requestedThreads : availableThreads;
    int effectiveThreads = desiredThreads < 1 ? 1 : desiredThreads;
    if (effectiveThreads > ctx.count) {
        effectiveThreads = ctx.count;
    }
    omp_set_dynamic(0);
    omp_set_num_threads(effectiveThreads);
#else
    int availableThreads = 1;
    int effectiveThreads = 1;
    if (requestedThreads > 1) {
        fprintf(stderr, "Warning: binary built without OpenMP support; using 1 thread.\n");
    }
#endif

    double setupStart = WallNowSeconds();

    if (!BuildLengthOrder(&ctx) || !BuildGeneratorTables(&ctx) || !BuildExpressionTables(&ctx)) {
        fprintf(stderr, "Failed to build precomputed tables.\n");
        FreeContext(&ctx);
        return 1;
    }

    char dataFile[32];
    snprintf(dataFile, sizeof(dataFile), "S%d.txt", n);
    if (!LoadNonTrivialData(&ctx, dataFile)) {
        FreeContext(&ctx);
        return 1;
    }

    if (!MaybeBuildKLBasisCache(&ctx)) {
        fprintf(stderr, "Failed while preparing KL cache.\n");
        FreeContext(&ctx);
        return 1;
    }

    int sampleCount = sampleRowsRequested;
    if (sampleCount > ctx.count) {
        sampleCount = ctx.count;
    }

    int *weights = (int *)calloc((size_t)ctx.count, sizeof(int));
    int *sampleRows = (int *)calloc((size_t)sampleCount, sizeof(int));
    if (!weights || !sampleRows) {
        fprintf(stderr, "Out of memory preparing sample rows.\n");
        free(weights);
        free(sampleRows);
        FreeContext(&ctx);
        return 1;
    }

    if (!BuildRowWeights(&ctx, weights)) {
        fprintf(stderr, "Failed to build row weights.\n");
        free(weights);
        free(sampleRows);
        FreeContext(&ctx);
        return 1;
    }

    int actualSampleCount = 0;
    if (rowScheduleMode == ROW_SCHEDULE_WEIGHTED) {
        actualSampleCount = BuildWeightedSampleRows(weights, ctx.count, sampleCount, sampleRows);
    } else {
        actualSampleCount = BuildUniformSampleRows(ctx.count, sampleCount, sampleRows);
    }

    if (actualSampleCount <= 0) {
        fprintf(stderr, "Could not create sample rows.\n");
        free(weights);
        free(sampleRows);
        FreeContext(&ctx);
        return 1;
    }

    double setupEnd = WallNowSeconds();

    int workerError = 0;
    int doneRows = 0;
    int64_t donePairs = 0;
    int64_t fullPairs = (int64_t)ctx.count * (int64_t)ctx.count;

    double sampleStart = WallNowSeconds();

    #pragma omp parallel
    {
        SparseHecke Cx, Cy, product, tmp1, tmp2, residual, basisTmp;
        DenseAccum accum;
        LeftMulCache leftMulCache;
        Poly57 *coeffs = NULL;
        int localInitOk = 1;

        memset(&Cx, 0, sizeof(Cx));
        memset(&Cy, 0, sizeof(Cy));
        memset(&product, 0, sizeof(product));
        memset(&tmp1, 0, sizeof(tmp1));
        memset(&tmp2, 0, sizeof(tmp2));
        memset(&residual, 0, sizeof(residual));
        memset(&basisTmp, 0, sizeof(basisTmp));
        memset(&accum, 0, sizeof(accum));
        memset(&leftMulCache, 0, sizeof(leftMulCache));

        if (!SparseInit(&Cx, ctx.count) || !SparseInit(&Cy, ctx.count) || !SparseInit(&product, ctx.count) ||
            !SparseInit(&tmp1, ctx.count) || !SparseInit(&tmp2, ctx.count) || !SparseInit(&residual, ctx.count) ||
            !SparseInit(&basisTmp, ctx.count) || !DenseAccumInit(&accum, ctx.count) ||
            !LeftMulCacheInit(&leftMulCache, ctx.count)) {
            localInitOk = 0;
        }

        coeffs = (Poly57 *)calloc((size_t)ctx.count, sizeof(Poly57));
        if (!coeffs) {
            localInitOk = 0;
        }

        if (!localInitOk) {
            #pragma omp critical
            {
                workerError = 1;
            }
        }

        #pragma omp for schedule(static)
        for (int si = 0; si < actualSampleCount; si++) {
            if (workerError) {
                continue;
            }

            int x = sampleRows[si];

            const SparseHecke *CxPtr = NULL;
            SparseHecke CxView;
            if (ctx.useCache) {
                const Poly57 *rowX = ctx.klCache + (size_t)x * (size_t)ctx.count;
                int sx = ctx.klSupportOffsets[x];
                int ex = ctx.klSupportOffsets[x + 1];
                CxView.count = ctx.count;
                CxView.coeff = (Poly57 *)rowX;
                CxView.support = ctx.klSupportIndices + sx;
                CxView.posMap = NULL;
                CxView.supportSize = ex - sx;
                CxPtr = &CxView;
            } else {
                BuildKLElementSparse(&ctx, x, &Cx);
                CxPtr = &Cx;
            }

            int64_t localPairs = 0;
            for (int y = 0; y < ctx.count; y++) {
                const SparseHecke *CyPtr = NULL;
                SparseHecke CyView;
                if (ctx.useCache) {
                    const Poly57 *rowY = ctx.klCache + (size_t)y * (size_t)ctx.count;
                    int sy = ctx.klSupportOffsets[y];
                    int ey = ctx.klSupportOffsets[y + 1];
                    CyView.count = ctx.count;
                    CyView.coeff = (Poly57 *)rowY;
                    CyView.support = ctx.klSupportIndices + sy;
                    CyView.posMap = NULL;
                    CyView.supportSize = ey - sy;
                    CyPtr = &CyView;
                } else {
                    BuildKLElementSparse(&ctx, y, &Cy);
                    CyPtr = &Cy;
                }

                MultiplyHeckeSparse(
                    &ctx,
                    y,
                    CxPtr,
                    CyPtr,
                    &product,
                    &tmp1,
                    &tmp2,
                    &accum,
                    &leftMulCache
                );

                DecomposeToKLBasis(&ctx, &product, coeffs, &residual, &basisTmp);
                localPairs++;
            }

            #pragma omp atomic
            donePairs += localPairs;

            int sampledNow;
            #pragma omp atomic capture
            sampledNow = ++doneRows;

            if (sampledNow % 8 == 0 || sampledNow == actualSampleCount) {
                #pragma omp critical
                {
                    double elapsed = WallNowSeconds() - sampleStart;
                    fprintf(stderr, "Sample progress: %d/%d rows, elapsed %.1fs\n", sampledNow, actualSampleCount, elapsed);
                }
            }
        }

        free(coeffs);
        LeftMulCacheFree(&leftMulCache);
        DenseAccumFree(&accum);
        SparseFree(&Cx);
        SparseFree(&Cy);
        SparseFree(&product);
        SparseFree(&tmp1);
        SparseFree(&tmp2);
        SparseFree(&residual);
        SparseFree(&basisTmp);
    }

    double sampleWall = WallNowSeconds() - sampleStart;

    if (workerError || donePairs <= 0) {
        fprintf(stderr, "Sampling failed (worker initialization/execution error).\n");
        free(weights);
        free(sampleRows);
        FreeContext(&ctx);
        return 1;
    }

    double setupWall = setupEnd - setupStart;
    double estLoopWall = sampleWall * ((double)fullPairs / (double)donePairs);
    double estTotalWall = setupWall + estLoopWall;

    printf("\n=== Runtime Estimate for S%d ===\n", n);
    printf("threads               : %d (available %d)\n", effectiveThreads, availableThreads);
    printf("row-schedule          : %s\n", rowScheduleMode == ROW_SCHEDULE_WEIGHTED ? "weighted" : "block");
    printf("sample rows           : %d/%d\n", actualSampleCount, ctx.count);
    printf("sampled pairs         : %lld / %lld (%.2f%%)\n",
           (long long)donePairs,
           (long long)fullPairs,
           100.0 * (double)donePairs / (double)fullPairs);
    printf("setup wall            : %.3fs\n", setupWall);
    printf("sample loop wall      : %.3fs\n", sampleWall);
    printf("estimated full loop   : %.3fs\n", estLoopWall);
    printf("estimated total       : %.3fs (%.2f min, %.2f h)\n",
           estTotalWall,
           estTotalWall / 60.0,
           estTotalWall / 3600.0);

    if (baselineN > 0 && baselineSeconds > 0.0 && baselineN != n) {
        double fn = (double)fac(n);
        double fb = (double)fac(baselineN);
        if (fn > 0.0 && fb > 0.0) {
            double a = log(estTotalWall / baselineSeconds) / log(fn / fb);
            printf("\nFitted exponent a from baseline S%d -> S%d: %.4f\n", baselineN, n, a);

            if (n + 1 <= 9) {
                double f1 = (double)fac(n + 1);
                double t1 = estTotalWall * pow(f1 / fn, a);
                printf("predicted S%d total   : %.3fs (%.2f min, %.2f h, %.2f d)\n",
                       n + 1, t1, t1 / 60.0, t1 / 3600.0, t1 / 86400.0);
            }
            if (n + 2 <= 9) {
                double f2 = (double)fac(n + 2);
                double t2 = estTotalWall * pow(f2 / fn, a);
                printf("predicted S%d total   : %.3fs (%.2f min, %.2f h, %.2f d)\n",
                       n + 2, t2, t2 / 60.0, t2 / 3600.0, t2 / 86400.0);
            }
        }
    }

    free(weights);
    free(sampleRows);
    FreeContext(&ctx);
    return 0;
}
