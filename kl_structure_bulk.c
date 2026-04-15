#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "lehmer.h"
#include "laurent.h"

typedef struct {
    int64_t key;
    Laurent_t poly;
} NonTrivialEntry;

typedef struct {
    int n;
    int count;
    int maxLength;
    int useCache;

    int *lengths;
    int *indicesByLengthDesc;
    int *exprLen;
    int *exprGenerators;
    char *permStrings;

    NonTrivialEntry *entries;
    size_t entryCount;

    Laurent_t *klBasisCache;
} BulkContext;

static int g_sortCount = 0;
static int *g_sortLengths = NULL;

static int64_t PairKey(int count, int y, int w) {
    return (int64_t)y * (int64_t)count + (int64_t)w;
}

static Laurent_t NegateLaurent(Laurent_t p) {
    Laurent_t result = ZeroInitializeLaurent();
    for (int i = 0; i < 57; i++) {
        result.coeff[i] = -p.coeff[i];
    }
    return result;
}

static int ExponentToPos(int exponent) {
    int pos = 28 + exponent;
    if (pos < 0 || pos >= 57) {
        return -1;
    }
    return pos;
}

static void AddAtExponent(Laurent_t *poly, int exponent, int coeff) {
    int pos = ExponentToPos(exponent);
    if (pos >= 0) {
        poly->coeff[pos] += coeff;
    }
}

static int ParsePermutationString(int n, const char *token, int *outIndex) {
    size_t tokenLen = strlen(token);
    if ((int)tokenLen != n || n > 8) {
        return 0;
    }

    int seen[8] = {0};
    char perm[8] = {0};
    for (int i = 0; i < n; i++) {
        if (!isdigit((unsigned char)token[i])) {
            return 0;
        }

        int digit = token[i] - '0';
        if (digit < 0 || digit >= n || seen[digit]) {
            return 0;
        }

        seen[digit] = 1;
        perm[i] = (char)digit;
    }

    *outIndex = PermToIndex(n, perm);
    return 1;
}

static void BuildPermutationStringTable(BulkContext *ctx) {
    int n = ctx->n;
    int count = ctx->count;
    ctx->permStrings = (char *)calloc((size_t)count * (size_t)(n + 1), sizeof(char));
    if (!ctx->permStrings) {
        return;
    }

    for (int i = 0; i < count; i++) {
        char perm[8] = {0};
        IndexToPerm(n, i, perm);
        char *target = ctx->permStrings + (size_t)i * (size_t)(n + 1);
        for (int j = 0; j < n; j++) {
            target[j] = (char)('0' + perm[j]);
        }
        target[n] = '\0';
    }
}

static const char *PermString(const BulkContext *ctx, int index) {
    return ctx->permStrings + (size_t)index * (size_t)(ctx->n + 1);
}

static int CompareEntryByKey(const void *a, const void *b) {
    const NonTrivialEntry *ea = (const NonTrivialEntry *)a;
    const NonTrivialEntry *eb = (const NonTrivialEntry *)b;
    if (ea->key < eb->key) {
        return -1;
    }
    if (ea->key > eb->key) {
        return 1;
    }
    return 0;
}

static const Laurent_t *FindNonTrivial(const BulkContext *ctx, int y, int w) {
    if (ctx->entryCount == 0) {
        return NULL;
    }

    int64_t key = PairKey(ctx->count, y, w);
    NonTrivialEntry needle;
    needle.key = key;
    NonTrivialEntry *found = (NonTrivialEntry *)bsearch(
        &needle,
        ctx->entries,
        ctx->entryCount,
        sizeof(NonTrivialEntry),
        CompareEntryByKey
    );

    if (!found) {
        return NULL;
    }
    return &found->poly;
}

static Laurent_t KLCoefficient(const BulkContext *ctx, int y, int w) {
    Laurent_t result = ZeroInitializeLaurent();

    if (y == w) {
        result.coeff[28] = 1;
        return result;
    }

    if (!BruhatSmaller2(ctx->n, y, w)) {
        return result;
    }

    const Laurent_t *nonTrivial = FindNonTrivial(ctx, y, w);
    if (nonTrivial) {
        return *nonTrivial;
    }

    int diff = ctx->lengths[w] - ctx->lengths[y];
    AddAtExponent(&result, diff, 1);
    return result;
}

static int LoadNonTrivialData(BulkContext *ctx, const char *filename) {
    FILE *fptr = fopen(filename, "r");
    if (!fptr) {
        fprintf(stderr, "Could not open %s\n", filename);
        return 0;
    }

    size_t cap = 1024;
    size_t used = 0;
    NonTrivialEntry *entries = (NonTrivialEntry *)calloc(cap, sizeof(NonTrivialEntry));
    if (!entries) {
        fclose(fptr);
        fprintf(stderr, "Out of memory while loading %s\n", filename);
        return 0;
    }

    char line[4096];
    while (fgets(line, sizeof(line), fptr)) {
        char yPerm[32] = {0};
        char wPerm[32] = {0};
        char coeffToken[3072] = {0};

        if (sscanf(line, "%31s %31s %3071s", yPerm, wPerm, coeffToken) != 3) {
            continue;
        }

        int y = -1;
        int w = -1;
        if (!ParsePermutationString(ctx->n, yPerm, &y) || !ParsePermutationString(ctx->n, wPerm, &w)) {
            continue;
        }

        Laurent_t poly = ZeroInitializeLaurent();
        int diff = ctx->lengths[w] - ctx->lengths[y];

        char coeffBuffer[3072] = {0};
        snprintf(coeffBuffer, sizeof(coeffBuffer), "%s", coeffToken);
        char *part = strtok(coeffBuffer, ",");
        int k = 0;
        while (part != NULL) {
            int coeff = atoi(part);
            int exponent = diff - 2 * k;
            AddAtExponent(&poly, exponent, coeff);
            part = strtok(NULL, ",");
            k++;
        }

        if (used == cap) {
            size_t newCap = cap * 2;
            NonTrivialEntry *newEntries = (NonTrivialEntry *)realloc(entries, newCap * sizeof(NonTrivialEntry));
            if (!newEntries) {
                free(entries);
                fclose(fptr);
                fprintf(stderr, "Out of memory while growing nontrivial entry storage.\n");
                return 0;
            }
            entries = newEntries;
            cap = newCap;
        }

        entries[used].key = PairKey(ctx->count, y, w);
        entries[used].poly = poly;
        used++;
    }

    fclose(fptr);

    qsort(entries, used, sizeof(NonTrivialEntry), CompareEntryByKey);
    ctx->entries = entries;
    ctx->entryCount = used;
    return 1;
}

static void BuildReducedExpressionCache(BulkContext *ctx) {
    int count = ctx->count;
    int maxLen = ctx->maxLength;
    ctx->exprLen = (int *)calloc((size_t)count, sizeof(int));
    ctx->exprGenerators = (int *)calloc((size_t)count * (size_t)maxLen, sizeof(int));
    if (!ctx->exprLen || !ctx->exprGenerators) {
        return;
    }

    for (int i = 0; i < count; i++) {
        char expression[40] = {0};
        int len = ReducedExpression(ctx->n, i, expression);
        ctx->exprLen[i] = len;

        for (int k = 0; k < len; k++) {
            int generator = (int)expression[k];
            ctx->exprGenerators[(size_t)i * (size_t)maxLen + (size_t)k] = fac(ctx->n - generator);
        }
    }
}

static int CompareIndicesByLengthDesc(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    int la = g_sortLengths[ia];
    int lb = g_sortLengths[ib];
    if (la != lb) {
        return lb - la;
    }
    return ia - ib;
}

static void BuildDescendingLengthOrder(BulkContext *ctx) {
    ctx->indicesByLengthDesc = (int *)calloc((size_t)ctx->count, sizeof(int));
    if (!ctx->indicesByLengthDesc) {
        return;
    }

    for (int i = 0; i < ctx->count; i++) {
        ctx->indicesByLengthDesc[i] = i;
    }

    g_sortCount = ctx->count;
    g_sortLengths = ctx->lengths;
    qsort(ctx->indicesByLengthDesc, (size_t)g_sortCount, sizeof(int), CompareIndicesByLengthDesc);
    g_sortCount = 0;
    g_sortLengths = NULL;
}

static void FillKLElement(const BulkContext *ctx, int w, Laurent_t element[]) {
    for (int y = 0; y < ctx->count; y++) {
        element[y] = KLCoefficient(ctx, y, w);
    }
}

static Laurent_t *CopyHeckeElement(int count, const Laurent_t source[]) {
    Laurent_t *copy = (Laurent_t *)calloc((size_t)count, sizeof(Laurent_t));
    if (!copy) {
        return NULL;
    }
    for (int i = 0; i < count; i++) {
        copy[i] = source[i];
    }
    return copy;
}

static Laurent_t *MultiplySimpleHeckeLeft(const BulkContext *ctx, int s, const Laurent_t H[]) {
    int count = ctx->count;
    Laurent_t *product = (Laurent_t *)calloc((size_t)count, sizeof(Laurent_t));
    if (!product) {
        return NULL;
    }

    Laurent_t lowerTerm = ZeroInitializeLaurent();
    lowerTerm.coeff[27] = 1;
    lowerTerm.coeff[29] = -1;

    for (int i = 0; i < count; i++) {
        if (!HasNonZero(&H[i])) {
            continue;
        }

        int composition = MultiplyIndex(ctx->n, s, i);
        product[composition] = SumLaurent(product[composition], H[i]);

        if (ctx->lengths[composition] != ctx->lengths[i] + 1) {
            Laurent_t corr = MultiplyLaurent(lowerTerm, H[i]);
            product[i] = SumLaurent(product[i], corr);
        }
    }

    return product;
}

static Laurent_t *LeftMultiplyByIndex(const BulkContext *ctx, int index, const Laurent_t H[]) {
    Laurent_t *term = CopyHeckeElement(ctx->count, H);
    if (!term) {
        return NULL;
    }

    int len = ctx->exprLen[index];
    for (int k = 0; k < len; k++) {
        int s = ctx->exprGenerators[(size_t)index * (size_t)ctx->maxLength + (size_t)k];
        Laurent_t *next = MultiplySimpleHeckeLeft(ctx, s, term);
        free(term);
        if (!next) {
            return NULL;
        }
        term = next;
    }

    return term;
}

static Laurent_t *MultiplyHeckeElements(const BulkContext *ctx, const Laurent_t H1[], const Laurent_t H2[]) {
    int count = ctx->count;
    Laurent_t *product = (Laurent_t *)calloc((size_t)count, sizeof(Laurent_t));
    if (!product) {
        return NULL;
    }

    for (int i = 0; i < count; i++) {
        if (!HasNonZero(&H1[i])) {
            continue;
        }

        Laurent_t *left = LeftMultiplyByIndex(ctx, i, H2);
        if (!left) {
            free(product);
            return NULL;
        }

        Laurent_t scalar = H1[i];
        for (int j = 0; j < count; j++) {
            if (!HasNonZero(&left[j])) {
                continue;
            }
            Laurent_t scaled = MultiplyLaurent(scalar, left[j]);
            product[j] = SumLaurent(product[j], scaled);
        }

        free(left);
    }

    return product;
}

static const Laurent_t *GetKLBasisElement(const BulkContext *ctx, int w, Laurent_t *scratch) {
    if (ctx->useCache) {
        return ctx->klBasisCache + (size_t)w * (size_t)ctx->count;
    }

    FillKLElement(ctx, w, scratch);
    return scratch;
}

static Laurent_t *DecomposeToKLBasis(const BulkContext *ctx, const Laurent_t heckeElement[]) {
    int count = ctx->count;
    Laurent_t *residual = CopyHeckeElement(count, heckeElement);
    Laurent_t *coeffs = (Laurent_t *)calloc((size_t)count, sizeof(Laurent_t));
    Laurent_t *scratch = NULL;

    if (!residual || !coeffs) {
        free(residual);
        free(coeffs);
        return NULL;
    }

    if (!ctx->useCache) {
        scratch = (Laurent_t *)calloc((size_t)count, sizeof(Laurent_t));
        if (!scratch) {
            free(residual);
            free(coeffs);
            return NULL;
        }
    }

    for (int idx = 0; idx < count; idx++) {
        int z = ctx->indicesByLengthDesc[idx];
        Laurent_t top = residual[z];
        if (!HasNonZero(&top)) {
            continue;
        }

        coeffs[z] = top;
        Laurent_t minusTop = NegateLaurent(top);
        const Laurent_t *basisZ = GetKLBasisElement(ctx, z, scratch);

        for (int y = 0; y < count; y++) {
            if (!HasNonZero(&basisZ[y])) {
                continue;
            }
            Laurent_t corr = MultiplyLaurent(minusTop, basisZ[y]);
            residual[y] = SumLaurent(residual[y], corr);
        }
    }

    free(scratch);
    free(residual);
    return coeffs;
}

static void FprintLaurent(FILE *f, Laurent_t poly) {
    int first = 1;
    if (!HasNonZero(&poly)) {
        fprintf(f, "0");
        return;
    }

    for (int i = 0; i < 57; i++) {
        if (poly.coeff[i] == 0) {
            continue;
        }
        if (!first) {
            fprintf(f, " + ");
        }
        fprintf(f, "%d*v^%d", poly.coeff[i], i - 28);
        first = 0;
    }
}

static void WriteTextProduct(FILE *out, const BulkContext *ctx, int x, int y, const Laurent_t coeffs[]) {
    fprintf(out, "C_{%s} * C_{%s} = ", PermString(ctx, x), PermString(ctx, y));

    int first = 1;
    for (int z = 0; z < ctx->count; z++) {
        if (!HasNonZero(&coeffs[z])) {
            continue;
        }

        if (!first) {
            fprintf(out, " + ");
        }

        fprintf(out, "(");
        FprintLaurent(out, coeffs[z]);
        fprintf(out, ")C_{%s}", PermString(ctx, z));
        first = 0;
    }

    if (first) {
        fprintf(out, "0");
    }
    fprintf(out, "\n");
}

static void WriteCsvTerms(FILE *out, const BulkContext *ctx, int x, int y, const Laurent_t coeffs[]) {
    int wrote = 0;
    for (int z = 0; z < ctx->count; z++) {
        if (!HasNonZero(&coeffs[z])) {
            continue;
        }

        fprintf(out, "\"%s\",\"%s\",\"%s\",\"", PermString(ctx, x), PermString(ctx, y), PermString(ctx, z));
        FprintLaurent(out, coeffs[z]);
        fprintf(out, "\"\n");
        wrote = 1;
    }

    if (!wrote) {
        fprintf(out, "\"%s\",\"%s\",\"\",\"0\"\n", PermString(ctx, x), PermString(ctx, y));
    }
}

static int MaybeBuildBasisCache(BulkContext *ctx) {
    const size_t cacheLimitBytes = (size_t)512 * (size_t)1024 * (size_t)1024;
    size_t needed = (size_t)ctx->count * (size_t)ctx->count * sizeof(Laurent_t);
    if (needed > cacheLimitBytes) {
        ctx->useCache = 0;
        return 1;
    }

    ctx->klBasisCache = (Laurent_t *)calloc((size_t)ctx->count * (size_t)ctx->count, sizeof(Laurent_t));
    if (!ctx->klBasisCache) {
        ctx->useCache = 0;
        return 1;
    }

    ctx->useCache = 1;
    for (int w = 0; w < ctx->count; w++) {
        Laurent_t *target = ctx->klBasisCache + (size_t)w * (size_t)ctx->count;
        FillKLElement(ctx, w, target);
    }
    return 1;
}

static void FreeContext(BulkContext *ctx) {
    free(ctx->lengths);
    free(ctx->indicesByLengthDesc);
    free(ctx->exprLen);
    free(ctx->exprGenerators);
    free(ctx->permStrings);
    free(ctx->entries);
    free(ctx->klBasisCache);
}

static void PrintUsage(const char *prog) {
    printf("Usage:\n");
    printf("  %s n [output_txt] [output_csv]\n", prog);
    printf("\n");
    printf("n must be in {4,5,6,7,8,9}.\n");
    printf("If outputs are omitted, defaults are S<n>_structure_constants_bulk.txt/.csv.\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2 || argc > 4) {
        PrintUsage(argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    if (n < 4 || n > 9) {
        fprintf(stderr, "n must be between 4 and 9.\n");
        return 1;
    }

    char dataFile[32];
    snprintf(dataFile, sizeof(dataFile), "S%d.txt", n);

    char defaultTxt[128];
    char defaultCsv[128];
    snprintf(defaultTxt, sizeof(defaultTxt), "S%d_structure_constants_bulk.txt", n);
    snprintf(defaultCsv, sizeof(defaultCsv), "S%d_structure_constants_bulk.csv", n);

    const char *txtPath = (argc >= 3) ? argv[2] : defaultTxt;
    const char *csvPath = (argc >= 4) ? argv[3] : defaultCsv;

    BulkContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = n;
    ctx.count = fac(n);
    ctx.maxLength = n * (n - 1) / 2;
    ctx.useCache = 0;

    ctx.lengths = (int *)calloc((size_t)ctx.count, sizeof(int));
    if (!ctx.lengths) {
        fprintf(stderr, "Out of memory allocating lengths.\n");
        FreeContext(&ctx);
        return 1;
    }

    for (int i = 0; i < ctx.count; i++) {
        ctx.lengths[i] = IndexToLength(n, i);
    }

    BuildDescendingLengthOrder(&ctx);
    BuildReducedExpressionCache(&ctx);
    BuildPermutationStringTable(&ctx);
    if (!ctx.indicesByLengthDesc || !ctx.exprLen || !ctx.exprGenerators || !ctx.permStrings) {
        fprintf(stderr, "Out of memory while building caches.\n");
        FreeContext(&ctx);
        return 1;
    }

    if (!LoadNonTrivialData(&ctx, dataFile)) {
        FreeContext(&ctx);
        return 1;
    }

    if (!MaybeBuildBasisCache(&ctx)) {
        FreeContext(&ctx);
        return 1;
    }

    FILE *txtOut = fopen(txtPath, "w");
    FILE *csvOut = fopen(csvPath, "w");
    if (!txtOut || !csvOut) {
        fprintf(stderr, "Could not open output files.\n");
        if (txtOut) fclose(txtOut);
        if (csvOut) fclose(csvOut);
        FreeContext(&ctx);
        return 1;
    }

    fprintf(txtOut, "# All KL-basis products in S%d\n", n);
    fprintf(txtOut, "# Format: C_{x} * C_{y} = sum_z h_{x,y}^z(v) C_{z}\n\n");
    fprintf(csvOut, "\"x\",\"y\",\"z\",\"h\"\n");

    clock_t start = clock();
    Laurent_t *scratchX = NULL;
    Laurent_t *scratchY = NULL;
    if (!ctx.useCache) {
        scratchX = (Laurent_t *)calloc((size_t)ctx.count, sizeof(Laurent_t));
        scratchY = (Laurent_t *)calloc((size_t)ctx.count, sizeof(Laurent_t));
        if (!scratchX || !scratchY) {
            fprintf(stderr, "Out of memory in non-cached mode.\n");
            free(scratchX);
            free(scratchY);
            fclose(txtOut);
            fclose(csvOut);
            FreeContext(&ctx);
            return 1;
        }
    }

    int64_t totalPairs = (int64_t)ctx.count * (int64_t)ctx.count;
    int64_t donePairs = 0;

    for (int x = 0; x < ctx.count; x++) {
        const Laurent_t *Cx = NULL;
        if (ctx.useCache) {
            Cx = ctx.klBasisCache + (size_t)x * (size_t)ctx.count;
        } else {
            FillKLElement(&ctx, x, scratchX);
            Cx = scratchX;
        }

        for (int y = 0; y < ctx.count; y++) {
            const Laurent_t *Cy = NULL;
            if (ctx.useCache) {
                Cy = ctx.klBasisCache + (size_t)y * (size_t)ctx.count;
            } else {
                FillKLElement(&ctx, y, scratchY);
                Cy = scratchY;
            }

            Laurent_t *productH = MultiplyHeckeElements(&ctx, Cx, Cy);
            Laurent_t *coeffs = DecomposeToKLBasis(&ctx, productH);
            if (!productH || !coeffs) {
                fprintf(stderr, "Out of memory during multiplication/decomposition.\n");
                free(productH);
                free(coeffs);
                free(scratchX);
                free(scratchY);
                fclose(txtOut);
                fclose(csvOut);
                FreeContext(&ctx);
                return 1;
            }

            WriteTextProduct(txtOut, &ctx, x, y, coeffs);
            WriteCsvTerms(csvOut, &ctx, x, y, coeffs);

            free(productH);
            free(coeffs);
            donePairs++;
        }

        if ((x + 1) % 5 == 0 || x + 1 == ctx.count) {
            double elapsedSec = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
            fprintf(stderr, "Progress: %d/%d rows of left factors, %.2f%% pairs, elapsed %.1fs\n",
                x + 1,
                ctx.count,
                100.0 * (double)donePairs / (double)totalPairs,
                elapsedSec
            );
        }
    }

    double elapsedSec = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
    fprintf(stderr, "Done. Wrote %s and %s in %.1fs\n", txtPath, csvPath, elapsedSec);

    free(scratchX);
    free(scratchY);
    fclose(txtOut);
    fclose(csvOut);
    FreeContext(&ctx);
    return 0;
}