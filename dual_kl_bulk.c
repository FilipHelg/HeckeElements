#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    int coeff[57];
    int minPos;
    int maxPos;
} Poly57;

typedef struct {
    int64_t key;
    Poly57 poly;
} NonTrivialEntry;

typedef struct {
    int n;
    int count;
    int maxLength;

    int *lengths;
    int *order;
    int *inverseIndex;

    int involutionCount;
    int *involutions;
    int *involutionWeights;

    NonTrivialEntry *entries;
    size_t entryCount;
} BulkContext;

typedef struct {
    FILE *file;
    char path[512];
} ThreadOutput;

typedef struct {
    char magic[8];
    uint32_t version;
    uint32_t n;
    uint32_t count;
    uint32_t involutionCount;
    uint32_t recordFormat;
    uint32_t reserved[3];
} DualBulkHeader;

static int *g_sortLengths = NULL;

static double WallNowSeconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static double CpuNowSeconds(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

static int FacLocal(int n) {
    int value = 1;
    for (int i = 2; i <= n; i++) {
        value *= i;
    }
    return value;
}

static void PolyZero(Poly57 *p) {
    for (int i = 0; i < 57; i++) {
        p->coeff[i] = 0;
    }
    p->minPos = 57;
    p->maxPos = -1;
}

static int PolyIsZero(const Poly57 *p) {
    return p->maxPos < p->minPos;
}

static void PolyNormalize(Poly57 *p) {
    if (PolyIsZero(p)) {
        p->minPos = 57;
        p->maxPos = -1;
        return;
    }

    int left = p->minPos;
    int right = p->maxPos;
    while (left <= right && p->coeff[left] == 0) {
        left++;
    }
    while (right >= left && p->coeff[right] == 0) {
        right--;
    }

    if (left > right) {
        p->minPos = 57;
        p->maxPos = -1;
    } else {
        p->minPos = left;
        p->maxPos = right;
    }
}

static void PolyCopy(Poly57 *dst, const Poly57 *src) {
    *dst = *src;
}

static void PolyNegate(Poly57 *dst, const Poly57 *src) {
    if (PolyIsZero(src)) {
        PolyZero(dst);
        return;
    }

    *dst = *src;
    for (int i = src->minPos; i <= src->maxPos; i++) {
        dst->coeff[i] = -dst->coeff[i];
    }
}

static void PolyAddInplace(Poly57 *dst, const Poly57 *src) {
    if (PolyIsZero(src)) {
        return;
    }

    if (PolyIsZero(dst)) {
        *dst = *src;
        return;
    }

    for (int i = src->minPos; i <= src->maxPos; i++) {
        if (src->coeff[i] != 0) {
            dst->coeff[i] += src->coeff[i];
        }
    }

    if (src->minPos < dst->minPos) {
        dst->minPos = src->minPos;
    }
    if (src->maxPos > dst->maxPos) {
        dst->maxPos = src->maxPos;
    }
    PolyNormalize(dst);
}

static void PolyAdd(Poly57 *out, const Poly57 *a, const Poly57 *b) {
    PolyCopy(out, a);
    PolyAddInplace(out, b);
}

static void PolyMultiply(const Poly57 *a, const Poly57 *b, Poly57 *out) {
    PolyZero(out);
    if (PolyIsZero(a) || PolyIsZero(b)) {
        return;
    }

    int seen = 0;
    int minPos = 57;
    int maxPos = -1;

    for (int i = a->minPos; i <= a->maxPos; i++) {
        int ai = a->coeff[i];
        if (ai == 0) {
            continue;
        }
        for (int j = b->minPos; j <= b->maxPos; j++) {
            int bj = b->coeff[j];
            if (bj == 0) {
                continue;
            }

            int pos = i + j - 28;
            if (pos < 0 || pos >= 57) {
                continue;
            }

            out->coeff[pos] += ai * bj;
            if (!seen) {
                minPos = pos;
                maxPos = pos;
                seen = 1;
            } else {
                if (pos < minPos) minPos = pos;
                if (pos > maxPos) maxPos = pos;
            }
        }
    }

    if (!seen) {
        PolyZero(out);
        return;
    }

    out->minPos = minPos;
    out->maxPos = maxPos;
    PolyNormalize(out);
}

static int PolyTermCount(const Poly57 *p) {
    if (PolyIsZero(p)) {
        return 0;
    }

    int count = 0;
    for (int i = p->minPos; i <= p->maxPos; i++) {
        if (p->coeff[i] != 0) {
            count++;
        }
    }
    return count;
}

static void IndexToPermSafe(int n, int index, int perm[]) {
    int code[8] = {0};
    int remaining = index;

    for (int i = 0; i < n - 1; i++) {
        int fact = 1;
        for (int j = 1; j <= n - 1 - i; j++) {
            fact *= j;
        }
        code[i] = remaining / fact;
        remaining %= fact;
    }

    int available[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    int availableCount = n;
    for (int i = 0; i < n; i++) {
        int j = code[i];
        if (j < 0 || j >= availableCount) {
            perm[i] = 0;
            continue;
        }
        perm[i] = available[j];
        for (int k = j; k < availableCount - 1; k++) {
            available[k] = available[k + 1];
        }
        availableCount--;
    }
}

static int PermToIndexSafe(int n, const int perm[]) {
    int code[8] = {0};
    for (int i = 0; i < n; i++) {
        int sum = 0;
        for (int j = i + 1; j < n; j++) {
            if (perm[j] < perm[i]) {
                sum++;
            }
        }
        code[i] = sum;
    }

    int index = 0;
    for (int i = 0; i < n; i++) {
        int fact = 1;
        for (int j = 1; j <= n - 1 - i; j++) {
            fact *= j;
        }
        index += code[i] * fact;
    }
    return index;
}

static int IndexToLengthSafe(int n, int index) {
    int code[8] = {0};
    int remaining = index;
    int length = 0;

    for (int i = 0; i < n - 1; i++) {
        int fact = 1;
        for (int j = 1; j <= n - 1 - i; j++) {
            fact *= j;
        }
        code[i] = remaining / fact;
        remaining %= fact;
        length += code[i];
    }

    return length;
}

static int CheckBruhatSmallerSafe(int n, int w1, int w2) {
    int permx[8] = {0};
    int permy[8] = {0};
    IndexToPermSafe(n, w1, permx);
    IndexToPermSafe(n, w2, permy);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int rx = 0;
            int ry = 0;
            for (int k = 0; k <= i; k++) {
                if (permx[k] >= j) {
                    rx++;
                }
                if (permy[k] >= j) {
                    ry++;
                }
            }
            if (rx > ry) {
                return 0;
            }
        }
    }

    return 1;
}

static int CompareEntries(const void *a, const void *b) {
    const NonTrivialEntry *ea = (const NonTrivialEntry *)a;
    const NonTrivialEntry *eb = (const NonTrivialEntry *)b;
    if (ea->key < eb->key) return -1;
    if (ea->key > eb->key) return 1;
    return 0;
}

static int ParsePermutationString(int n, const char *token, int *outIndex) {
    if ((int)strlen(token) != n || n > 8) {
        return 0;
    }

    int seen[8] = {0};
    int perm[8] = {0};
    for (int i = 0; i < n; i++) {
        if (!isdigit((unsigned char)token[i])) {
            return 0;
        }

        int d = token[i] - '0';
        if (d < 0 || d >= n || seen[d]) {
            return 0;
        }

        seen[d] = 1;
        perm[i] = d;
    }

    *outIndex = PermToIndexSafe(n, perm);
    return 1;
}

static int64_t PairKey(int count, int y, int w) {
    return (int64_t)y * (int64_t)count + (int64_t)w;
}

static const Poly57 *FindNonTrivial(const BulkContext *ctx, int y, int w) {
    if (ctx->entryCount == 0) {
        return NULL;
    }

    NonTrivialEntry needle;
    needle.key = PairKey(ctx->count, y, w);

    NonTrivialEntry *found = (NonTrivialEntry *)bsearch(
        &needle,
        ctx->entries,
        ctx->entryCount,
        sizeof(NonTrivialEntry),
        CompareEntries
    );

    if (!found) {
        return NULL;
    }
    return &found->poly;
}

static void AddAtExponent(Poly57 *poly, int exponent, int coeff) {
    int pos = 28 + exponent;
    if (pos < 0 || pos >= 57 || coeff == 0) {
        return;
    }
    poly->coeff[pos] += coeff;
    if (poly->minPos > poly->maxPos) {
        poly->minPos = pos;
        poly->maxPos = pos;
    } else {
        if (pos < poly->minPos) poly->minPos = pos;
        if (pos > poly->maxPos) poly->maxPos = pos;
    }
}

static Poly57 KLCoefficient(
    const BulkContext *ctx,
    int y,
    int w
) {
    Poly57 p;
    PolyZero(&p);

    if (y == w) {
        p.coeff[28] = 1;
        p.minPos = 28;
        p.maxPos = 28;
        return p;
    }

    if (!CheckBruhatSmallerSafe(ctx->n, y, w)) {
        return p;
    }

    const Poly57 *nt = FindNonTrivial(ctx, y, w);
    if (nt) {
        return *nt;
    }

    int exponent = ctx->lengths[w] - ctx->lengths[y];
    int pos = 28 + exponent;
    if (pos >= 0 && pos < 57) {
        p.coeff[pos] = 1;
        p.minPos = pos;
        p.maxPos = pos;
    }
    return p;
}

static int LoadNonTrivialData(
    BulkContext *ctx,
    const char *filename
) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Could not open %s\n", filename);
        return 0;
    }

    size_t cap = 1024;
    size_t used = 0;
    NonTrivialEntry *entries = (NonTrivialEntry *)calloc(cap, sizeof(NonTrivialEntry));
    if (!entries) {
        fclose(f);
        fprintf(stderr, "Out of memory while loading nontrivial KL data.\n");
        return 0;
    }

    char line[16384];
    while (fgets(line, sizeof(line), f)) {
        char yPerm[32] = {0};
        char wPerm[32] = {0};
        char coeffToken[16000] = {0};

        if (sscanf(line, "%31s %31s %15999s", yPerm, wPerm, coeffToken) != 3) {
            continue;
        }

        int y = -1;
        int w = -1;
        if (!ParsePermutationString(ctx->n, yPerm, &y) || !ParsePermutationString(ctx->n, wPerm, &w)) {
            continue;
        }

        Poly57 p;
        PolyZero(&p);
        int diff = ctx->lengths[w] - ctx->lengths[y];

        char coeffBuf[16000] = {0};
        snprintf(coeffBuf, sizeof(coeffBuf), "%s", coeffToken);

        char *part = strtok(coeffBuf, ",");
        int k = 0;
        while (part != NULL) {
            int coeff = atoi(part);
            int exponent = diff - 2 * k;
            AddAtExponent(&p, exponent, coeff);
            part = strtok(NULL, ",");
            k++;
        }

        if (used == cap) {
            size_t newCap = cap * 2;
            NonTrivialEntry *newEntries = (NonTrivialEntry *)realloc(entries, newCap * sizeof(NonTrivialEntry));
            if (!newEntries) {
                free(entries);
                fclose(f);
                fprintf(stderr, "Out of memory while growing nontrivial KL table.\n");
                return 0;
            }
            entries = newEntries;
            cap = newCap;
        }

        entries[used].key = PairKey(ctx->count, y, w);
        entries[used].poly = p;
        used++;
    }

    fclose(f);

    qsort(entries, used, sizeof(NonTrivialEntry), CompareEntries);
    ctx->entries = entries;
    ctx->entryCount = used;
    return 1;
}

static int CompareByLengthThenIndex(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    int la = g_sortLengths[ia];
    int lb = g_sortLengths[ib];
    if (la != lb) {
        return la - lb;
    }
    return ia - ib;
}

static int BuildLengthOrder(BulkContext *ctx) {
    ctx->order = (int *)calloc((size_t)ctx->count, sizeof(int));
    if (!ctx->order) {
        return 0;
    }

    for (int i = 0; i < ctx->count; i++) {
        ctx->order[i] = i;
    }

    g_sortLengths = ctx->lengths;
    qsort(ctx->order, (size_t)ctx->count, sizeof(int), CompareByLengthThenIndex);
    g_sortLengths = NULL;
    return 1;
}

static int BuildInverseTable(BulkContext *ctx) {
    ctx->inverseIndex = (int *)calloc((size_t)ctx->count, sizeof(int));
    if (!ctx->inverseIndex) {
        return 0;
    }

    for (int w = 0; w < ctx->count; w++) {
        int perm[8] = {0};
        int inv[8] = {0};
        IndexToPermSafe(ctx->n, w, perm);
        for (int i = 0; i < ctx->n; i++) {
            inv[perm[i]] = i;
        }
        ctx->inverseIndex[w] = PermToIndexSafe(ctx->n, inv);
    }

    return 1;
}

static int BuildInvolutionList(BulkContext *ctx) {
    ctx->involutions = (int *)calloc((size_t)ctx->count, sizeof(int));
    ctx->involutionWeights = (int *)calloc((size_t)ctx->count, sizeof(int));
    if (!ctx->involutions || !ctx->involutionWeights) {
        return 0;
    }

    int used = 0;
    for (int w = 0; w < ctx->count; w++) {
        if (ctx->inverseIndex[w] == w) {
            ctx->involutions[used] = w;
            ctx->involutionWeights[used] = ctx->lengths[w] + 1;
            used++;
        }
    }

    ctx->involutionCount = used;
    return 1;
}

static void BuildThreadRanges(
    int count,
    int threadCount,
    const int *weights,
    int *starts,
    int *ends
) {
    if (threadCount <= 0) {
        return;
    }

    int64_t totalWeight = 0;
    for (int i = 0; i < count; i++) {
        totalWeight += (int64_t)(weights[i] > 0 ? weights[i] : 1);
    }

    int current = 0;
    int64_t prefix = 0;
    for (int t = 0; t < threadCount; t++) {
        starts[t] = current;

        if (t == threadCount - 1) {
            ends[t] = count;
            continue;
        }

        int64_t targetPrefix = ((int64_t)(t + 1) * totalWeight) / (int64_t)threadCount;
        while (current < count && prefix < targetPrefix) {
            prefix += (int64_t)(weights[current] > 0 ? weights[current] : 1);
            current++;
        }

        if (current < starts[t]) {
            current = starts[t];
        }
        ends[t] = current;
    }
}

static int MakeDirectoryIfNeeded(const char *path) {
#ifdef _WIN32
    if (_mkdir(path) == 0 || errno == EEXIST) {
        return 1;
    }
#else
    if (mkdir(path, 0777) == 0 || errno == EEXIST) {
        return 1;
    }
#endif
    return 0;
}

static void BuildThreadTempPath(char *out, size_t outSize, const char *dir, int tid) {
    snprintf(out, outSize, "%s/thread_%03d.bin", dir, tid);
}

static int OpenThreadOutput(ThreadOutput *out, const char *dir, int tid) {
    BuildThreadTempPath(out->path, sizeof(out->path), dir, tid);
    out->file = fopen(out->path, "wb");
    if (!out->file) {
        return 0;
    }
    setvbuf(out->file, NULL, _IOFBF, 1 << 20);
    return 1;
}

static void CloseThreadOutput(ThreadOutput *out) {
    if (out->file) {
        fclose(out->file);
        out->file = NULL;
    }
}

static int MergeBinaryFile(FILE *dst, const char *path) {
    FILE *src = fopen(path, "rb");
    if (!src) {
        return 0;
    }

    char buffer[1 << 15];
    size_t got;
    while ((got = fread(buffer, 1, sizeof(buffer), src)) > 0) {
        if (fwrite(buffer, 1, got, dst) != got) {
            fclose(src);
            return 0;
        }
    }

    if (ferror(src)) {
        fclose(src);
        return 0;
    }

    fclose(src);
    return 1;
}

static void RemoveTempFiles(const char *dir, int threadCount) {
    for (int t = 0; t < threadCount; t++) {
        char path[512];
        BuildThreadTempPath(path, sizeof(path), dir, t);
        remove(path);
    }
}

static void ClearCoeffArray(Poly57 *coeffs, int count) {
    for (int i = 0; i < count; i++) {
        PolyZero(&coeffs[i]);
    }
}

static void SolveDualElement(
    const BulkContext *ctx,
    int x,
    Poly57 *coeffs,
    int *active,
    int *activeCount
) {
    ClearCoeffArray(coeffs, ctx->count);
    *activeCount = 0;

    for (int oi = 0; oi < ctx->count; oi++) {
        int y = ctx->order[oi];

        Poly57 sum;
        PolyZero(&sum);

        for (int ai = 0; ai < *activeCount; ai++) {
            int z = active[ai];
            if (ctx->lengths[z] > ctx->lengths[y]) {
                continue;
            }

            Poly57 pyz = KLCoefficient(ctx, z, y);
            if (PolyIsZero(&pyz) || PolyIsZero(&coeffs[z])) {
                continue;
            }

            Poly57 term;
            PolyMultiply(&pyz, &coeffs[z], &term);
            PolyAddInplace(&sum, &term);
        }

        Poly57 rhs;
        PolyZero(&rhs);
        if (y == x) {
            rhs.coeff[28] = 1;
            rhs.minPos = 28;
            rhs.maxPos = 28;
        }

        Poly57 negSum;
        PolyNegate(&negSum, &sum);
        PolyAdd(&coeffs[y], &rhs, &negSum);

        if (!PolyIsZero(&coeffs[y])) {
            active[(*activeCount)++] = y;
        }
    }
}

static int WriteU32(FILE *out, uint32_t value) {
    return fwrite(&value, sizeof(value), 1, out) == 1;
}

static int WriteI32(FILE *out, int32_t value) {
    return fwrite(&value, sizeof(value), 1, out) == 1;
}

static int WriteI8(FILE *out, int8_t value) {
    return fwrite(&value, sizeof(value), 1, out) == 1;
}

static int WritePolyBinary(FILE *out, const Poly57 *p) {
    uint8_t termCount = (uint8_t)PolyTermCount(p);
    if (fwrite(&termCount, sizeof(termCount), 1, out) != 1) {
        return 0;
    }

    if (termCount == 0) {
        return 1;
    }

    for (int i = p->minPos; i <= p->maxPos; i++) {
        int coeff = p->coeff[i];
        if (coeff == 0) {
            continue;
        }

        int8_t exponent = (int8_t)(i - 28);
        if (!WriteI8(out, exponent) || !WriteI32(out, (int32_t)coeff)) {
            return 0;
        }
    }

    return 1;
}

static int WriteDualRecord(FILE *out, int x, const Poly57 *coeffs, int count) {
    uint32_t supportCount = 0;
    for (int z = 0; z < count; z++) {
        if (!PolyIsZero(&coeffs[z])) {
            supportCount++;
        }
    }

    if (!WriteU32(out, (uint32_t)x) || !WriteU32(out, supportCount)) {
        return 0;
    }

    for (int z = 0; z < count; z++) {
        if (PolyIsZero(&coeffs[z])) {
            continue;
        }

        if (!WriteU32(out, (uint32_t)z) || !WritePolyBinary(out, &coeffs[z])) {
            return 0;
        }
    }

    return 1;
}

static int ParseThreadsAndOptions(
    int argc,
    char *argv[],
    const char **outputPath,
    int *requestedThreads,
    int *useWeightedSchedule,
    int *benchEnabled
) {
    const char *positional[1] = {0};
    int positionalCount = 0;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--bench") == 0) {
            *benchEnabled = 1;
            continue;
        }

        if (strcmp(argv[i], "--threads") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--threads requires a positive integer.\n");
                return 0;
            }
            *requestedThreads = atoi(argv[++i]);
            if (*requestedThreads <= 0) {
                fprintf(stderr, "--threads requires a positive integer.\n");
                return 0;
            }
            continue;
        }

        if (strcmp(argv[i], "--row-schedule") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--row-schedule requires one of: block, weighted\n");
                return 0;
            }

            const char *mode = argv[++i];
            if (strcmp(mode, "block") == 0) {
                *useWeightedSchedule = 0;
            } else if (strcmp(mode, "weighted") == 0) {
                *useWeightedSchedule = 1;
            } else {
                fprintf(stderr, "Unknown --row-schedule mode '%s' (use block or weighted).\n", mode);
                return 0;
            }
            continue;
        }

        if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return 0;
        }

        if (positionalCount < 1) {
            positional[positionalCount++] = argv[i];
        } else {
            fprintf(stderr, "Too many positional arguments.\n");
            return 0;
        }
    }

    if (positionalCount == 1) {
        *outputPath = positional[0];
    }

    return 1;
}

static void PrintUsage(const char *prog) {
    printf("Usage:\n");
    printf("  %s n [output_bin] [--threads N] [--row-schedule block|weighted] [--bench]\n", prog);
    printf("\n");
    printf("n must be between 3 and 8.\n");
    printf("Default output is dual_kl_bulk_n<n>.bin.\n");
    printf("The program computes dual KL-basis elements for involutions only and writes a sparse binary cache.\n");
}

static void FreeContext(BulkContext *ctx) {
    free(ctx->lengths);
    free(ctx->order);
    free(ctx->inverseIndex);
    free(ctx->involutions);
    free(ctx->involutionWeights);
    free(ctx->entries);
}

int main(int argc, char *argv[]) {
    if (argc < 2 || argc > 8) {
        PrintUsage(argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    if (n < 3 || n > 8) {
        fprintf(stderr, "n must be between 3 and 8 for this bulk generator.\n");
        return 1;
    }

    char defaultOutput[128];
    snprintf(defaultOutput, sizeof(defaultOutput), "dual_kl_bulk_n%d.bin", n);
    const char *outputPath = defaultOutput;
    int requestedThreads = 0;
    int useWeightedSchedule = 1;
    int benchEnabled = 0;

    if (!ParseThreadsAndOptions(argc, argv, &outputPath, &requestedThreads, &useWeightedSchedule, &benchEnabled)) {
        return 1;
    }

    BulkContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = n;
    ctx.count = FacLocal(n);
    ctx.maxLength = n * (n - 1) / 2;

    ctx.lengths = (int *)calloc((size_t)ctx.count, sizeof(int));
    if (!ctx.lengths) {
        fprintf(stderr, "Out of memory for length table.\n");
        FreeContext(&ctx);
        return 1;
    }

    for (int i = 0; i < ctx.count; i++) {
        ctx.lengths[i] = IndexToLengthSafe(n, i);
    }

    if (!BuildLengthOrder(&ctx) || !BuildInverseTable(&ctx) || !BuildInvolutionList(&ctx)) {
        fprintf(stderr, "Failed to build permutation tables.\n");
        FreeContext(&ctx);
        return 1;
    }

    if (n >= 4) {
        char dataFile[32];
        snprintf(dataFile, sizeof(dataFile), "S%d.txt", n);
        if (!LoadNonTrivialData(&ctx, dataFile)) {
            FreeContext(&ctx);
            return 1;
        }
    }

    int threadCount = 1;
#ifdef _OPENMP
    int availableThreads = omp_get_max_threads();
    int desiredThreads = requestedThreads > 0 ? requestedThreads : availableThreads;
    threadCount = desiredThreads < 1 ? 1 : desiredThreads;
    if (threadCount > ctx.involutionCount) {
        threadCount = ctx.involutionCount;
    }
    omp_set_dynamic(0);
    omp_set_num_threads(threadCount);
#else
    int availableThreads = 1;
    if (requestedThreads > 0 && requestedThreads != 1) {
        fprintf(stderr, "Warning: binary built without OpenMP support; using 1 thread.\n");
    }
#endif

    if (threadCount < 1) {
        threadCount = 1;
    }

    int *targetStarts = (int *)calloc((size_t)threadCount, sizeof(int));
    int *targetEnds = (int *)calloc((size_t)threadCount, sizeof(int));
    if (!targetStarts || !targetEnds) {
        fprintf(stderr, "Out of memory while preparing thread ranges.\n");
        free(targetStarts);
        free(targetEnds);
        FreeContext(&ctx);
        return 1;
    }

    BuildThreadRanges(ctx.involutionCount, threadCount, ctx.involutionWeights, targetStarts, targetEnds);

    char tempDir[256];
    snprintf(tempDir, sizeof(tempDir), "tmp_dual_bulk_n%d_%lld", n, (long long)time(NULL));
    if (!MakeDirectoryIfNeeded(tempDir)) {
        fprintf(stderr, "Could not create temporary directory %s\n", tempDir);
        free(targetStarts);
        free(targetEnds);
        FreeContext(&ctx);
        return 1;
    }

    ThreadOutput *threadOutputs = (ThreadOutput *)calloc((size_t)threadCount, sizeof(ThreadOutput));
    if (!threadOutputs) {
        fprintf(stderr, "Out of memory allocating thread output states.\n");
        RemoveTempFiles(tempDir, threadCount);
        free(targetStarts);
        free(targetEnds);
        FreeContext(&ctx);
        return 1;
    }

    int workerError = 0;
    double startWall = WallNowSeconds();
    double startCpu = CpuNowSeconds();
    int globalCompleted = 0;
    const int progressInterval = 5;

    #pragma omp parallel
    {
        int tid = 0;
    #ifdef _OPENMP
        tid = omp_get_thread_num();
    #endif

        ThreadOutput *out = &threadOutputs[tid];
        Poly57 *coeffs = (Poly57 *)calloc((size_t)ctx.count, sizeof(Poly57));
        int *active = (int *)calloc((size_t)ctx.count, sizeof(int));
        int activeCount = 0;
        int localError = 0;

        if (!coeffs || !active || !OpenThreadOutput(out, tempDir, tid)) {
            localError = 1;
        }

        if (localError) {
            #pragma omp critical
            {
                workerError = 1;
            }
        } else {
            int start = targetStarts[tid];
            int end = targetEnds[tid];
            for (int i = start; i < end; i++) {
                int x = ctx.involutions[i];
                SolveDualElement(&ctx, x, coeffs, active, &activeCount);
                if (!WriteDualRecord(out->file, x, coeffs, ctx.count)) {
                    #pragma omp critical
                    {
                        workerError = 1;
                    }
                    localError = 1;
                    break;
                }
                else {
                    int done = 0;
                    #pragma omp atomic capture
                    done = ++globalCompleted;
                    if (done % progressInterval == 0) {
                        #pragma omp critical
                        {
                            fprintf(stderr, "Progress: %d/%d involutions completed\n", done, ctx.involutionCount);
                        }
                    }
                }
            }
        }

        if (out->file) {
            fclose(out->file);
            out->file = NULL;
        }
        free(coeffs);
        free(active);

        if (localError) {
            #pragma omp critical
            {
                workerError = 1;
            }
        }
    }

    if (workerError) {
        fprintf(stderr, "Bulk computation failed before merge.\n");
        RemoveTempFiles(tempDir, threadCount);
        for (int t = 0; t < threadCount; t++) {
            CloseThreadOutput(&threadOutputs[t]);
        }
        free(threadOutputs);
        free(targetStarts);
        free(targetEnds);
        FreeContext(&ctx);
        return 1;
    }

    FILE *finalOut = fopen(outputPath, "wb");
    if (!finalOut) {
        fprintf(stderr, "Could not open output file %s\n", outputPath);
        RemoveTempFiles(tempDir, threadCount);
        free(threadOutputs);
        free(targetStarts);
        free(targetEnds);
        FreeContext(&ctx);
        return 1;
    }
    setvbuf(finalOut, NULL, _IOFBF, 1 << 20);

    DualBulkHeader header;
    memset(&header, 0, sizeof(header));
    memcpy(header.magic, "DKBULK1", 7);
    header.version = 1;
    header.n = (uint32_t)n;
    header.count = (uint32_t)ctx.count;
    header.involutionCount = (uint32_t)ctx.involutionCount;
    header.recordFormat = 1;

    if (fwrite(&header, sizeof(header), 1, finalOut) != 1) {
        fprintf(stderr, "Failed while writing output header.\n");
        fclose(finalOut);
        RemoveTempFiles(tempDir, threadCount);
        free(threadOutputs);
        free(targetStarts);
        free(targetEnds);
        FreeContext(&ctx);
        return 1;
    }

    for (int t = 0; t < threadCount; t++) {
        if (!MergeBinaryFile(finalOut, threadOutputs[t].path)) {
            fprintf(stderr, "Failed while merging %s.\n", threadOutputs[t].path);
            fclose(finalOut);
            RemoveTempFiles(tempDir, threadCount);
            free(threadOutputs);
            free(targetStarts);
            free(targetEnds);
            FreeContext(&ctx);
            return 1;
        }
    }

    fclose(finalOut);
    RemoveTempFiles(tempDir, threadCount);

#ifdef _WIN32
    _rmdir(tempDir);
#else
    rmdir(tempDir);
#endif

    double elapsedWall = WallNowSeconds() - startWall;
    double elapsedCpu = CpuNowSeconds() - startCpu;

    if (benchEnabled) {
        fprintf(stderr, "Bulk dual KL generation complete for S%d.\n", n);
        fprintf(stderr, "Targets (involutions): %d\n", ctx.involutionCount);
        fprintf(stderr, "Output: %s\n", outputPath);
        fprintf(stderr, "Threads: %d (available: %d)\n", threadCount,
    #ifdef _OPENMP
            availableThreads
    #else
            1
    #endif
        );
        fprintf(stderr, "Elapsed wall time: %.3fs\n", elapsedWall);
        fprintf(stderr, "Elapsed CPU  time: %.3fs\n", elapsedCpu);
    }

    free(threadOutputs);
    free(targetStarts);
    free(targetEnds);
    FreeContext(&ctx);
    return 0;
}
