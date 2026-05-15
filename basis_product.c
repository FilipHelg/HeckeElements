#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

#define KL_INVOLUTION_CHECK_SKIP_INCLUDE
#define main kl_involution_check_embedded_main
#include "kl_involution_check.c"
#undef main

#include <stdio.h>
#include <string.h>

/* forward declaration for exact division helper used in dual decomposition */
static int PolyDividesExact(const Poly57 *divisor, const Poly57 *numerator, Poly57 *quotient);

static void BasisProductPrintUsage(const char *prog) {
    fprintf(stderr,
        "Usage: %s n x y mode [--kl-data PATH] [--dual-bin PATH]\n"
        "  mode: t (standard H basis), k (KL basis), d (dual KL basis), kd (KL left, dual right), dk (dual left, KL right)\n"
        "  optional: --gap-words to print standard basis terms as T(reduced-word)\n"
        "  optional: --gap-factors to print the two input factors in GAP-style normalized form\n"
        "  optional: --gap-output to print the product in GAP-style normalized form\n"
        "Examples:\n"
        "  %s 5 7 9 t\n"
        "  %s 5 7 9 k\n"
        "  %s 5 7 9 d --dual-bin dual_kl_bulk_n5.bin\n"
        "  %s 5 7 9 kd --dual-bin dual_kl_bulk_n5.bin\n"
        "  %s 5 7 9 dk --dual-bin dual_kl_bulk_n5.bin\n",
        prog, prog, prog, prog, prog, prog);
}

static int ParseArgs(
    int argc,
    char *argv[],
    int *nOut,
    int *xOut,
    int *yOut,
    const char **modeOut,
    const char **klDataOut,
    const char **dualBinOut,
    int *gapWordsOut,
    int *printFactorsOut,
    int *gapFactorsOut,
    int *gapOutputOut
) {
    if (argc < 5) {
        return 0;
    }

    *nOut = atoi(argv[1]);
    *xOut = atoi(argv[2]);
    *yOut = atoi(argv[3]);
    *modeOut = argv[4];
    *klDataOut = NULL;
    *dualBinOut = NULL;
    *gapWordsOut = 0;
    *printFactorsOut = 0;
    *gapFactorsOut = 0;
    *gapOutputOut = 0;

    for (int i = 5; i < argc; i++) {
        if (strcmp(argv[i], "--kl-data") == 0) {
            if (i + 1 >= argc) return 0;
            *klDataOut = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "--dual-bin") == 0) {
            if (i + 1 >= argc) return 0;
            *dualBinOut = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "--gap-words") == 0) {
            *gapWordsOut = 1;
            continue;
        }
        if (strcmp(argv[i], "--print-factors") == 0) {
            *printFactorsOut = 1;
            continue;
        }
        if (strcmp(argv[i], "--gap-factors") == 0) {
            *printFactorsOut = 1;
            *gapFactorsOut = 1;
            continue;
        }
        if (strcmp(argv[i], "--gap-output") == 0) {
            *gapOutputOut = 1;
            continue;
        }
        return 0;
    }

    return 1;
}

static void BuildStandardBasisElement(int index, SparseHecke *out) {
    Poly57 one;
    PolyFromMonomial(&one, 0, 1);
    SparseClear(out);
    SparseAddPoly(out, index, &one);
}

static int CompareInt(const void *a, const void *b) {
    const int ia = *(const int *)a;
    const int ib = *(const int *)b;
    return ia - ib;
}

static void PrintReducedWordContents(int n, int index) {
    char expr[28] = {0};
    int len = ReducedExpression(n, index, expr);
    for (int i = 0; i < len; i++) {
        if (i > 0) {
            printf(",");
        }
        printf("%d", (int)expr[i]);
    }
}

static void PrintReducedWord(int n, int index) {
    printf("T(");
    PrintReducedWordContents(n, index);
    printf(")");
}

static void ShiftPolyByExponent(const Poly57 *src, int shift, Poly57 *dst) {
    PolyZero(dst);
    if (PolyIsZero(src)) {
        return;
    }

    int seen = 0;
    int minPos = 57;
    int maxPos = -1;

    for (int i = src->minPos; i <= src->maxPos; i++) {
        int c = src->coeff[i];
        if (c == 0) {
            continue;
        }
        int j = i + shift;
        if (j < 0 || j >= 57) {
            continue;
        }
        dst->coeff[j] += c;
        if (!seen) {
            minPos = j;
            maxPos = j;
            seen = 1;
        } else {
            if (j < minPos) minPos = j;
            if (j > maxPos) maxPos = j;
        }
    }

    if (!seen) {
        PolyZero(dst);
        return;
    }

    dst->minPos = minPos;
    dst->maxPos = maxPos;
    PolyNormalize(dst);
}

static void PrintSparseAsStandard(const FastContext *ctx, const SparseHecke *h, const char *permTable, int gapWords) {
    if (h->supportSize <= 0) {
        printf("0\n");
        return;
    }

    int *indices = (int *)calloc((size_t)h->supportSize, sizeof(int));
    if (!indices) {
        printf("<out of memory while printing>\n");
        return;
    }

    for (int i = 0; i < h->supportSize; i++) {
        indices[i] = h->support[i];
    }
    qsort(indices, (size_t)h->supportSize, sizeof(int), CompareInt);

    for (int i = 0; i < h->supportSize; i++) {
        int idx = indices[i];
        const Poly57 *coeff = &h->coeff[idx];
        if (i > 0) {
            printf(" + ");
        }
        if (coeff->minPos < coeff->maxPos) {
            printf("(");
            FprintPoly(stdout, coeff);
            printf(")");
        } else {
            FprintPoly(stdout, coeff);
        }
        if (gapWords) {
            printf("*");
            PrintReducedWord(ctx->n, idx);
        } else {
            printf("*H_{%s}", permTable + (size_t)idx * (size_t)(ctx->n + 1));
        }
    }
    printf("\n");

    free(indices);
}

static void PrintSparseAsGapNormalized(const FastContext *ctx, const SparseHecke *h) {
    if (h->supportSize <= 0) {
        printf("0\n");
        return;
    }

    int *indices = (int *)calloc((size_t)h->supportSize, sizeof(int));
    if (!indices) {
        printf("<out of memory while printing>\n");
        return;
    }

    for (int i = 0; i < h->supportSize; i++) {
        indices[i] = h->support[i];
    }
    qsort(indices, (size_t)h->supportSize, sizeof(int), CompareInt);

    for (int i = 0; i < h->supportSize; i++) {
        int idx = indices[i];
        Poly57 normalized;
        ShiftPolyByExponent(&h->coeff[idx], ctx->lengths[idx], &normalized);

        if (i > 0) {
            printf(" + ");
        }
        if (normalized.minPos < normalized.maxPos) {
            printf("(");
            FprintPoly(stdout, &normalized);
            printf(")");
        } else {
            FprintPoly(stdout, &normalized);
        }
        printf("*T(");
        PrintReducedWordContents(ctx->n, idx);
        printf(")");
    }
    printf("\n");

    free(indices);
}

static void PrintLabeledGapFactor(const char *label, const FastContext *ctx, const SparseHecke *value) {
    printf("%s\n", label);
    PrintSparseAsGapNormalized(ctx, value);
}

static void PrintLabeledFactor(const char *label, const FastContext *ctx, const SparseHecke *value, const char *permTable, int gapWords) {
    printf("%s\n", label);
    PrintSparseAsStandard(ctx, value, permTable, gapWords);
}

static void PrintCoeffArrayAsKLBasis(const FastContext *ctx, const Poly57 *coeffs, const char *permTable) {
    int first = 1;
    for (int i = 0; i < ctx->count; i++) {
        if (PolyIsZero(&coeffs[i])) {
            continue;
        }
        if (!first) {
            printf(" + ");
        }
        if (coeffs[i].minPos < coeffs[i].maxPos) {
            printf("(");
            FprintPoly(stdout, &coeffs[i]);
            printf(")");
        } else {
            FprintPoly(stdout, &coeffs[i]);
        }
        printf("*C_{%s}", permTable + (size_t)i * (size_t)(ctx->n + 1));
        first = 0;
    }
    if (first) {
        printf("0");
    }
    printf("\n");
}

static int PolyIsUnitMonomial(const Poly57 *p, int *expOut, int *signOut) {
    if (PolyIsZero(p)) {
        return 0;
    }

    int found = 0;
    int exp = 0;
    int sign = 0;
    for (int i = p->minPos; i <= p->maxPos; i++) {
        int c = p->coeff[i];
        if (c == 0) {
            continue;
        }
        if (!(c == 1 || c == -1)) {
            return 0;
        }
        found = 1;
        exp = i - 28;
        sign = c;
    }

    if (!found) {
        return 0;
    }

    if (expOut) *expOut = exp;
    if (signOut) *signOut = sign;
    return 1;
}

static int DivideByUnitMonomial(const Poly57 *numerator, int exp, int sign, Poly57 *out) {
    Poly57 inv;
    PolyZero(&inv);

    int invPos = 28 - exp;
    if (invPos < 0 || invPos >= 57) {
        return 0;
    }
    inv.coeff[invPos] = sign;
    inv.minPos = invPos;
    inv.maxPos = invPos;

    PolyMultiply(numerator, &inv, out);
    return 1;
}

static int DecomposeToDualBasis(
    const FastContext *ctx,
    const DualCache *dualCache,
    const SparseHecke *element,
    Poly57 *coeffs,
    SparseHecke *residual
) {
    ClearCoeffArray(coeffs, ctx->count);
    SparseAssign(residual, element);

    for (int ord = 0; ord < ctx->count; ord++) {
        int z = ctx->indicesByLengthDesc[ord];
        if (!dualCache->present[z]) {
            continue;
        }
        if (residual->posMap[z] < 0) {
            continue;
        }

        Poly57 top;
        PolyCopy(&top, &residual->coeff[z]);
        if (PolyIsZero(&top)) {
            continue;
        }

        const SparseHecke *dz = &dualCache->records[z];
        int diagPos = dz->posMap[z];
        if (diagPos < 0) {
            return 0;
        }

        /* Require exact divisibility in Z[v,v^-1]: find q with dz->coeff[z] * q = top. */
        Poly57 q;
        PolyZero(&q);
        if (!PolyDividesExact(&dz->coeff[z], &top, &q)) {
            return 0;
        }
        coeffs[z] = q;

        Poly57 minusQ;
        PolyNegate(&minusQ, &q);
        for (int p = 0; p < dz->supportSize; p++) {
            int y = dz->support[p];
            Poly57 corr;
            PolyMultiply(&minusQ, &dz->coeff[y], &corr);
            SparseAddPoly(residual, y, &corr);
        }
    }

    return 1;
}

static void PrintCoeffArrayAsDualBasis(const FastContext *ctx, const Poly57 *coeffs, const char *permTable) {
    int first = 1;
    for (int i = 0; i < ctx->count; i++) {
        if (PolyIsZero(&coeffs[i])) {
            continue;
        }
        if (!first) {
            printf(" + ");
        }
        if (coeffs[i].minPos < coeffs[i].maxPos) {
            printf("(");
            FprintPoly(stdout, &coeffs[i]);
            printf(")");
        } else {
            FprintPoly(stdout, &coeffs[i]);
        }
        printf("*d_{%s}", permTable + (size_t)i * (size_t)(ctx->n + 1));
        first = 0;
    }
    if (first) {
        printf("0");
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    int n = 0;
    int x = 0;
    int y = 0;
    const char *mode = NULL;
    const char *klDataPath = NULL;
    const char *dualBinPath = NULL;
    int gapWords = 0;
    int printFactors = 0;
    int gapFactors = 0;
    int gapOutput = 0;

    if (!ParseArgs(argc, argv, &n, &x, &y, &mode, &klDataPath, &dualBinPath, &gapWords, &printFactors, &gapFactors, &gapOutput)) {
        BasisProductPrintUsage(argv[0]);
        return 1;
    }

    if (n < 3 || n > 8) {
        fprintf(stderr, "n must be between 3 and 8.\n");
        return 1;
    }

    if (!(strcmp(mode, "t") == 0 || strcmp(mode, "k") == 0 || strcmp(mode, "d") == 0 || strcmp(mode, "kd") == 0 || strcmp(mode, "dk") == 0)) {
        fprintf(stderr, "mode must be one of: t, k, d, kd, dk\n");
        return 1;
    }

    int count = fac(n);
    if (x < 0 || x >= count || y < 0 || y >= count) {
        fprintf(stderr, "x and y must be in range [0, %d].\n", count - 1);
        return 1;
    }

    char defaultKlData[32];
    char defaultDualBin[128];
    if (!klDataPath) {
        snprintf(defaultKlData, sizeof(defaultKlData), "S%d.txt", n);
        klDataPath = defaultKlData;
    }
    if (!dualBinPath) {
        snprintf(defaultDualBin, sizeof(defaultDualBin), "dual_kl_bulk_n%d.bin", n);
        dualBinPath = defaultDualBin;
    }

    FastContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = n;
    ctx.count = count;
    ctx.maxLength = n * (n - 1) / 2;

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

    if (n >= 4) {
        if (!LoadNonTrivialData(&ctx, klDataPath)) {
            fprintf(stderr, "Failed to load KL data from %s.\n", klDataPath);
            FreeContext(&ctx);
            return 1;
        }
    }

    if (!MaybeBuildKLBasisCache(&ctx)) {
        fprintf(stderr, "Warning: KL basis cache could not be built. Continuing.\n");
    }

    SparseHecke left, right, product, tmp1, tmp2, residual, basisTmp;
    DenseAccum accum;
    if (!SparseInit(&left, ctx.count) || !SparseInit(&right, ctx.count) || !SparseInit(&product, ctx.count) ||
        !SparseInit(&tmp1, ctx.count) || !SparseInit(&tmp2, ctx.count) || !SparseInit(&residual, ctx.count) ||
        !SparseInit(&basisTmp, ctx.count) || !DenseAccumInit(&accum, ctx.count)) {
        fprintf(stderr, "Failed to allocate multiplication workspace.\n");
        FreeContext(&ctx);
        return 1;
    }

    DualCache dualCache;
    memset(&dualCache, 0, sizeof(dualCache));

    if (strcmp(mode, "t") == 0) {
        BuildStandardBasisElement(x, &left);
        BuildStandardBasisElement(y, &right);
    } else if (strcmp(mode, "k") == 0) {
        BuildKLElementSparse(&ctx, x, &left);
        BuildKLElementSparse(&ctx, y, &right);
    } else {
        /* mode == "d", "kd", or "dk" -> need dual cache */
        int loadedN = 0;
        if (!LoadDualCache(dualBinPath, ctx.count, &dualCache, &loadedN)) {
            fprintf(stderr, "Failed to load dual cache from %s.\n", dualBinPath);
            SparseFree(&left);
            SparseFree(&right);
            SparseFree(&product);
            SparseFree(&tmp1);
            SparseFree(&tmp2);
            SparseFree(&residual);
            SparseFree(&basisTmp);
            DenseAccumFree(&accum);
            FreeContext(&ctx);
            return 1;
        }
        if (loadedN != n) {
            fprintf(stderr, "Dual cache n=%d does not match requested n=%d.\n", loadedN, n);
            DualCacheFree(&dualCache);
            SparseFree(&left);
            SparseFree(&right);
            SparseFree(&product);
            SparseFree(&tmp1);
            SparseFree(&tmp2);
            SparseFree(&residual);
            SparseFree(&basisTmp);
            DenseAccumFree(&accum);
            FreeContext(&ctx);
            return 1;
        }
        if (strcmp(mode, "d") == 0) {
            if (!dualCache.present[x] || !dualCache.present[y]) {
                fprintf(stderr, "For mode d, x and y must be indices present in dual cache (typically involutions).\n");
                DualCacheFree(&dualCache);
                SparseFree(&left);
                SparseFree(&right);
                SparseFree(&product);
                SparseFree(&tmp1);
                SparseFree(&tmp2);
                SparseFree(&residual);
                SparseFree(&basisTmp);
                DenseAccumFree(&accum);
                FreeContext(&ctx);
                return 1;
            }
            SparseAssign(&left, &dualCache.records[x]);
            SparseAssign(&right, &dualCache.records[y]);
        } else if (strcmp(mode, "kd") == 0) {
            /* mode == kd: left = KL, right = dual */
            BuildKLElementSparse(&ctx, x, &left);
            if (!dualCache.present[y]) {
                fprintf(stderr, "For mode kd, right index y must be present in dual cache.\n");
                DualCacheFree(&dualCache);
                SparseFree(&left);
                SparseFree(&right);
                SparseFree(&product);
                SparseFree(&tmp1);
                SparseFree(&tmp2);
                SparseFree(&residual);
                SparseFree(&basisTmp);
                DenseAccumFree(&accum);
                FreeContext(&ctx);
                return 1;
            }
            SparseAssign(&right, &dualCache.records[y]);
        } else {
            /* mode == dk: left = dual, right = KL */
            if (!dualCache.present[x]) {
                fprintf(stderr, "For mode dk, left index x must be present in dual cache.\n");
                DualCacheFree(&dualCache);
                SparseFree(&left);
                SparseFree(&right);
                SparseFree(&product);
                SparseFree(&tmp1);
                SparseFree(&tmp2);
                SparseFree(&residual);
                SparseFree(&basisTmp);
                DenseAccumFree(&accum);
                FreeContext(&ctx);
                return 1;
            }
            SparseAssign(&left, &dualCache.records[x]);
            BuildKLElementSparse(&ctx, y, &right);
        }
    }

    SparseClear(&product);
    MultiplyHeckeSparse(&ctx, y, &left, &right, &product, &tmp1, &tmp2, &accum, NULL, 0);

    char *permTable = BuildPermStringTable(n, ctx.count);
    if (!permTable) {
        fprintf(stderr, "Failed to build permutation table.\n");
        if (strcmp(mode, "d") == 0 || strcmp(mode, "kd") == 0 || strcmp(mode, "dk") == 0) DualCacheFree(&dualCache);
        SparseFree(&left);
        SparseFree(&right);
        SparseFree(&product);
        SparseFree(&tmp1);
        SparseFree(&tmp2);
        SparseFree(&residual);
        SparseFree(&basisTmp);
        DenseAccumFree(&accum);
        FreeContext(&ctx);
        return 1;
    }

    printf("n=%d, x=%d, y=%d, mode=%s\n", n, x, y, mode);
        printf("Reduced expression for x=%d: ", x);
        PrintReducedWord(n, x);
        printf("\n");
        printf("Reduced expression for y=%d: ", y);
        PrintReducedWord(n, y);
        printf("\n");

    if (printFactors) {
        if (gapFactors) {
            PrintLabeledGapFactor("Left factor in GAP-style normalized form:", &ctx, &left);
            PrintLabeledGapFactor("Right factor in GAP-style normalized form:", &ctx, &right);
        } else {
            PrintLabeledFactor("Left factor in standard basis:", &ctx, &left, permTable, gapWords);
            PrintLabeledFactor("Right factor in standard basis:", &ctx, &right, permTable, gapWords);
        }
        }

    if (strcmp(mode, "t") == 0) {
        printf("H_x * H_y in standard basis:\n");
        if (gapOutput) {
            PrintSparseAsGapNormalized(&ctx, &product);
        } else {
            PrintSparseAsStandard(&ctx, &product, permTable, gapWords);
        }
    } else if (strcmp(mode, "k") == 0) {
        Poly57 *coeffs = (Poly57 *)calloc((size_t)ctx.count, sizeof(Poly57));
        if (!coeffs) {
            fprintf(stderr, "Out of memory for KL decomposition.\n");
            free(permTable);
            SparseFree(&left);
            SparseFree(&right);
            SparseFree(&product);
            SparseFree(&tmp1);
            SparseFree(&tmp2);
            SparseFree(&residual);
            SparseFree(&basisTmp);
            DenseAccumFree(&accum);
            FreeContext(&ctx);
            return 1;
        }

        DecomposeToKLBasis(&ctx, &product, coeffs, &residual, &basisTmp);
        printf("C_x * C_y in KL basis:\n");
        if (gapOutput) {
            PrintSparseAsGapNormalized(&ctx, &product);
        } else {
            PrintCoeffArrayAsKLBasis(&ctx, coeffs, permTable);
        }
        free(coeffs);
    } else if (strcmp(mode, "d") == 0) {
        Poly57 *dualCoeffs = (Poly57 *)calloc((size_t)ctx.count, sizeof(Poly57));
        if (!dualCoeffs) {
            fprintf(stderr, "Out of memory for dual decomposition.\n");
            free(permTable);
            DualCacheFree(&dualCache);
            SparseFree(&left);
            SparseFree(&right);
            SparseFree(&product);
            SparseFree(&tmp1);
            SparseFree(&tmp2);
            SparseFree(&residual);
            SparseFree(&basisTmp);
            DenseAccumFree(&accum);
            FreeContext(&ctx);
            return 1;
        }

        int decomposeOk = DecomposeToDualBasis(&ctx, &dualCache, &product, dualCoeffs, &residual);
        if (!decomposeOk) {
            fprintf(stderr, "Error: dual-basis decomposition failed (non-unit diagonal pivot in dual cache).\n");
            PrintSparseAsStandard(&ctx, &product, permTable, gapWords);
            free(dualCoeffs);
            DualCacheFree(&dualCache);
            SparseFree(&left);
            SparseFree(&right);
            SparseFree(&product);
            SparseFree(&tmp1);
            SparseFree(&tmp2);
            SparseFree(&residual);
            SparseFree(&basisTmp);
            DenseAccumFree(&accum);
            FreeContext(&ctx);
            return 2;
        }

        if (residual.supportSize > 0) {
            fprintf(stderr, "Error: dual-basis decomposition produced non-zero residual (unexpected).\n");
            PrintCoeffArrayAsDualBasis(&ctx, dualCoeffs, permTable);
            if (gapOutput) {
                PrintSparseAsGapNormalized(&ctx, &residual);
            } else {
                PrintSparseAsStandard(&ctx, &residual, permTable, gapWords);
            }
            free(dualCoeffs);
            DualCacheFree(&dualCache);
            SparseFree(&left);
            SparseFree(&right);
            SparseFree(&product);
            SparseFree(&tmp1);
            SparseFree(&tmp2);
            SparseFree(&residual);
            SparseFree(&basisTmp);
            DenseAccumFree(&accum);
            FreeContext(&ctx);
            return 3;
        }

        if (strcmp(mode, "d") == 0) {
            printf("d_x * d_y in dual basis:\n");
        }
        PrintCoeffArrayAsDualBasis(&ctx, dualCoeffs, permTable);

        free(dualCoeffs);
    } else if (strcmp(mode, "kd") == 0) {
        printf("C_x * d_y in standard basis:\n");
        if (gapOutput) {
            PrintSparseAsGapNormalized(&ctx, &product);
        } else {
            PrintSparseAsStandard(&ctx, &product, permTable, gapWords);
        }
    } else {
        printf("d_x * C_y in standard basis:\n");
        if (gapOutput) {
            PrintSparseAsGapNormalized(&ctx, &product);
        } else {
            PrintSparseAsStandard(&ctx, &product, permTable, gapWords);
        }
    }

    free(permTable);
    if (strcmp(mode, "d") == 0 || strcmp(mode, "kd") == 0 || strcmp(mode, "dk") == 0) {
        DualCacheFree(&dualCache);
    }
    SparseFree(&left);
    SparseFree(&right);
    SparseFree(&product);
    SparseFree(&tmp1);
    SparseFree(&tmp2);
    SparseFree(&residual);
    SparseFree(&basisTmp);
    DenseAccumFree(&accum);
    FreeContext(&ctx);
    return 0;
}

/* Check if `divisor * quotient == numerator` exactly over Z[v,v^-1]. If so set `quotient` and return 1. */
static int PolyDividesExact(const Poly57 *divisor, const Poly57 *numerator, Poly57 *quotient) {
    if (PolyIsZero(divisor)) return 0;

    Poly57 rem;
    PolyCopy(&rem, numerator);
    PolyZero(quotient);

    int degDiv = divisor->maxPos;
    int degRem = rem.maxPos;

    if (PolyIsZero(&rem)) {
        PolyZero(quotient);
        return 1;
    }

    while (!PolyIsZero(&rem) && rem.maxPos >= degDiv) {
        int degR = rem.maxPos;
        int ca = divisor->coeff[degDiv];
        int cr = rem.coeff[degR];

        if (ca == 0) return 0;
        if (cr % ca != 0) return 0; /* requires exact integer division */

        int t = cr / ca;
        int shift = degR - degDiv; /* index shift in array coordinates */

        /* place t into quotient at position `shift` (exponent = shift - 28 offset) */
        int qpos = shift + 28;
        if (qpos < 0 || qpos >= 57) return 0;
        quotient->coeff[qpos] += t;
        if (quotient->minPos > qpos) quotient->minPos = qpos;
        if (quotient->maxPos < qpos) quotient->maxPos = qpos;

        /* subtract divisor * t shifted by `shift` from remainder */
        for (int i = divisor->minPos; i <= divisor->maxPos; i++) {
            int dst = i + shift;
            if (dst < 0 || dst >= 57) return 0;
            rem.coeff[dst] -= divisor->coeff[i] * t;
        }

        PolyNormalize(&rem);
    }

    if (!PolyIsZero(&rem)) return 0;

    if (quotient->minPos > quotient->maxPos) {
        PolyZero(quotient);
    } else {
        PolyNormalize(quotient);
    }

    return 1;
}

/* forward declaration for exact division helper used earlier */
