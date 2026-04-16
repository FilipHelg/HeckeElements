#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "lehmer.h"

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
    int count;
    Poly57 *coeff;
    int *support;
    int *posMap;
    int supportSize;
} SparseHecke;

typedef struct {
    int count;
    Poly57 *coeff;
    int *support;
    unsigned int *stamp;
    int supportSize;
    unsigned int epoch;
} DenseAccum;

typedef struct {
    int slotCount;
    int *leftIndex;
    int *rightId;
    unsigned int *age;
    SparseHecke *value;
    unsigned int tick;
} LeftMulCache;

typedef struct {
    int n;
    int count;
    int maxLength;

    int *lengths;
    int *indicesByLengthDesc;

    int genCount;
    int genLehmer[8];
    int *leftAction;
    unsigned char *lengthIncrease;

    int *exprLen;
    unsigned char *exprG;

    NonTrivialEntry *entries;
    size_t entryCount;

    int useCache;
    Poly57 *klCache;
    int *klSupportOffsets;
    int *klSupportIndices;

    Poly57 lowerTerm;
} FastContext;

static int *g_sortLengths = NULL;

static void LeftMultiplyByIndex(
    const FastContext *ctx,
    int index,
    const SparseHecke *base,
    SparseHecke *result,
    SparseHecke *scratch
);

static void SparseAssign(SparseHecke *dst, const SparseHecke *src);

typedef struct {
    double buildXWall;
    double buildYWall;
    double multiplyWall;
    double decomposeWall;
    double writeTxtWall;
    double writeCsvWall;

    double buildXCpu;
    double buildYCpu;
    double multiplyCpu;
    double decomposeCpu;
    double writeTxtCpu;
    double writeCsvCpu;
} BenchStats;

static double WallNowSeconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static double CpuNowSeconds(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

static void PrintBenchReport(
    const FastContext *ctx,
    const BenchStats *stats,
    double loopWall,
    double loopCpu,
    int64_t donePairs
) {
    if (donePairs <= 0) {
        return;
    }

    double pairCount = (double)donePairs;
    double sumPhaseWall =
        stats->buildXWall + stats->buildYWall + stats->multiplyWall +
        stats->decomposeWall + stats->writeTxtWall + stats->writeCsvWall;
    double sumPhaseCpu =
        stats->buildXCpu + stats->buildYCpu + stats->multiplyCpu +
        stats->decomposeCpu + stats->writeTxtCpu + stats->writeCsvCpu;

    fprintf(stderr, "\n==== Benchmark report (S%d, %lld pairs) ====\n", ctx->n, (long long)donePairs);
    fprintf(stderr, "Total loop wall time: %.3fs\n", loopWall);
    fprintf(stderr, "Total loop CPU  time: %.3fs\n", loopCpu);

    if (loopWall > 0.0) {
        fprintf(stderr, "Estimated avg active cores: %.2f\n", loopCpu / loopWall);
    }

    fprintf(stderr, "\nPhase breakdown (wall):\n");
    fprintf(stderr, "  Build Cx     : %8.3fs (%.2f%%), %.3fus/pair\n",
        stats->buildXWall,
        (loopWall > 0.0) ? 100.0 * stats->buildXWall / loopWall : 0.0,
        1e6 * stats->buildXWall / pairCount);
    fprintf(stderr, "  Build Cy     : %8.3fs (%.2f%%), %.3fus/pair\n",
        stats->buildYWall,
        (loopWall > 0.0) ? 100.0 * stats->buildYWall / loopWall : 0.0,
        1e6 * stats->buildYWall / pairCount);
    fprintf(stderr, "  Multiply     : %8.3fs (%.2f%%), %.3fus/pair\n",
        stats->multiplyWall,
        (loopWall > 0.0) ? 100.0 * stats->multiplyWall / loopWall : 0.0,
        1e6 * stats->multiplyWall / pairCount);
    fprintf(stderr, "  Decompose    : %8.3fs (%.2f%%), %.3fus/pair\n",
        stats->decomposeWall,
        (loopWall > 0.0) ? 100.0 * stats->decomposeWall / loopWall : 0.0,
        1e6 * stats->decomposeWall / pairCount);
    fprintf(stderr, "  Write TXT    : %8.3fs (%.2f%%), %.3fus/pair\n",
        stats->writeTxtWall,
        (loopWall > 0.0) ? 100.0 * stats->writeTxtWall / loopWall : 0.0,
        1e6 * stats->writeTxtWall / pairCount);
    fprintf(stderr, "  Write CSV    : %8.3fs (%.2f%%), %.3fus/pair\n",
        stats->writeCsvWall,
        (loopWall > 0.0) ? 100.0 * stats->writeCsvWall / loopWall : 0.0,
        1e6 * stats->writeCsvWall / pairCount);

    fprintf(stderr, "\nPhase breakdown (CPU):\n");
    fprintf(stderr, "  Build Cx     : %8.3fs\n", stats->buildXCpu);
    fprintf(stderr, "  Build Cy     : %8.3fs\n", stats->buildYCpu);
    fprintf(stderr, "  Multiply     : %8.3fs\n", stats->multiplyCpu);
    fprintf(stderr, "  Decompose    : %8.3fs\n", stats->decomposeCpu);
    fprintf(stderr, "  Write TXT    : %8.3fs\n", stats->writeTxtCpu);
    fprintf(stderr, "  Write CSV    : %8.3fs\n", stats->writeCsvCpu);

    if (sumPhaseWall > 0.0) {
        fprintf(stderr, "\nInstrumented phase coverage: wall %.2f%%, CPU %.2f%%\n",
            100.0 * sumPhaseWall / loopWall,
            (loopCpu > 0.0) ? 100.0 * sumPhaseCpu / loopCpu : 0.0);
    }
}

static inline void PolyZero(Poly57 *p) {
    for (int i = 0; i < 57; i++) {
        p->coeff[i] = 0;
    }
    p->minPos = 57;
    p->maxPos = -1;
}

static inline int PolyIsZero(const Poly57 *p) {
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

static void PolyFromMonomial(Poly57 *p, int exponent, int coeff) {
    PolyZero(p);
    if (coeff == 0) {
        return;
    }

    int pos = 28 + exponent;
    if (pos < 0 || pos >= 57) {
        return;
    }

    p->coeff[pos] = coeff;
    p->minPos = pos;
    p->maxPos = pos;
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
        dst->coeff[i] = -src->coeff[i];
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

    int newMin = dst->minPos < src->minPos ? dst->minPos : src->minPos;
    int newMax = dst->maxPos > src->maxPos ? dst->maxPos : src->maxPos;

    for (int i = src->minPos; i <= src->maxPos; i++) {
        if (src->coeff[i] != 0) {
            dst->coeff[i] += src->coeff[i];
        }
    }

    dst->minPos = newMin;
    dst->maxPos = newMax;
    PolyNormalize(dst);
}

/* Specialized for multiplication by (v^-1 - v), i.e. shift right minus shift left in coeff-array coordinates. */
static void PolyMultiplyLowerTerm(const Poly57 *src, Poly57 *out) {
    PolyZero(out);
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

        int pDown = i - 1;
        if (pDown >= 0) {
            out->coeff[pDown] += c;
            if (!seen) {
                minPos = pDown;
                maxPos = pDown;
                seen = 1;
            } else {
                if (pDown < minPos) minPos = pDown;
                if (pDown > maxPos) maxPos = pDown;
            }
        }

        int pUp = i + 1;
        if (pUp < 57) {
            out->coeff[pUp] -= c;
            if (!seen) {
                minPos = pUp;
                maxPos = pUp;
                seen = 1;
            } else {
                if (pUp < minPos) minPos = pUp;
                if (pUp > maxPos) maxPos = pUp;
            }
        }
    }

    if (!seen) {
        return;
    }

    out->minPos = minPos;
    out->maxPos = maxPos;
    PolyNormalize(out);
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
        int ca = a->coeff[i];
        if (ca == 0) {
            continue;
        }
        for (int j = b->minPos; j <= b->maxPos; j++) {
            int cb = b->coeff[j];
            if (cb == 0) {
                continue;
            }

            int pos = i + j - 28;
            if (pos < 0 || pos >= 57) {
                continue;
            }

            out->coeff[pos] += ca * cb;
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

static void FprintPoly(FILE *f, const Poly57 *p) {
    if (PolyIsZero(p)) {
        fprintf(f, "0");
        return;
    }

    int first = 1;
    for (int i = p->minPos; i <= p->maxPos; i++) {
        if (p->coeff[i] == 0) {
            continue;
        }

        if (!first) {
            fprintf(f, " + ");
        }
        fprintf(f, "%d*v^%d", p->coeff[i], i - 28);
        first = 0;
    }
}

static int SparseInit(SparseHecke *h, int count) {
    h->count = count;
    h->coeff = (Poly57 *)calloc((size_t)count, sizeof(Poly57));
    h->support = (int *)calloc((size_t)count, sizeof(int));
    h->posMap = (int *)calloc((size_t)count, sizeof(int));
    h->supportSize = 0;

    if (!h->coeff || !h->support || !h->posMap) {
        free(h->coeff);
        free(h->support);
        free(h->posMap);
        h->coeff = NULL;
        h->support = NULL;
        h->posMap = NULL;
        return 0;
    }

    for (int i = 0; i < count; i++) {
        PolyZero(&h->coeff[i]);
        h->posMap[i] = -1;
    }
    return 1;
}

static void SparseFree(SparseHecke *h) {
    free(h->coeff);
    free(h->support);
    free(h->posMap);
    h->coeff = NULL;
    h->support = NULL;
    h->posMap = NULL;
    h->supportSize = 0;
}

static void SparseClear(SparseHecke *h) {
    for (int i = 0; i < h->supportSize; i++) {
        int idx = h->support[i];
        PolyZero(&h->coeff[idx]);
        h->posMap[idx] = -1;
    }
    h->supportSize = 0;
}

static void SparseRemoveIndex(SparseHecke *h, int idx) {
    int pos = h->posMap[idx];
    if (pos < 0) {
        return;
    }

    int lastIdx = h->support[h->supportSize - 1];
    h->support[pos] = lastIdx;
    h->posMap[lastIdx] = pos;

    h->supportSize--;
    h->posMap[idx] = -1;
    PolyZero(&h->coeff[idx]);
}

static void SparseAddPoly(SparseHecke *h, int idx, const Poly57 *p) {
    if (PolyIsZero(p)) {
        return;
    }

    int pos = h->posMap[idx];
    if (pos < 0) {
        int newPos = h->supportSize;
        h->support[newPos] = idx;
        h->posMap[idx] = newPos;
        h->supportSize++;
        PolyCopy(&h->coeff[idx], p);
        return;
    }

    PolyAddInplace(&h->coeff[idx], p);
    if (PolyIsZero(&h->coeff[idx])) {
        SparseRemoveIndex(h, idx);
    }
}

static int DenseAccumInit(DenseAccum *a, int count) {
    a->count = count;
    a->coeff = (Poly57 *)calloc((size_t)count, sizeof(Poly57));
    a->support = (int *)calloc((size_t)count, sizeof(int));
    a->stamp = (unsigned int *)calloc((size_t)count, sizeof(unsigned int));
    a->supportSize = 0;
    a->epoch = 1;

    if (!a->coeff || !a->support || !a->stamp) {
        free(a->coeff);
        free(a->support);
        free(a->stamp);
        a->coeff = NULL;
        a->support = NULL;
        a->stamp = NULL;
        return 0;
    }

    return 1;
}

static void DenseAccumFree(DenseAccum *a) {
    free(a->coeff);
    free(a->support);
    free(a->stamp);
    a->coeff = NULL;
    a->support = NULL;
    a->stamp = NULL;
    a->supportSize = 0;
    a->epoch = 1;
}

static void DenseAccumBegin(DenseAccum *a) {
    a->supportSize = 0;
    a->epoch++;
    if (a->epoch == 0) {
        memset(a->stamp, 0, (size_t)a->count * sizeof(unsigned int));
        a->epoch = 1;
    }
}

static void DenseAccumAdd(DenseAccum *a, int idx, const Poly57 *p) {
    if (PolyIsZero(p)) {
        return;
    }

    if (a->stamp[idx] != a->epoch) {
        a->stamp[idx] = a->epoch;
        a->support[a->supportSize++] = idx;
        a->coeff[idx] = *p;
        return;
    }

    PolyAddInplace(&a->coeff[idx], p);
}

static void DenseAccumCommitToSparse(DenseAccum *a, SparseHecke *out) {
    SparseClear(out);

    for (int k = 0; k < a->supportSize; k++) {
        int idx = a->support[k];
        Poly57 *p = &a->coeff[idx];
        if (PolyIsZero(p)) {
            continue;
        }

        int pos = out->supportSize;
        out->support[pos] = idx;
        out->posMap[idx] = pos;
        out->supportSize++;
        out->coeff[idx] = *p;
    }
}

static int LeftMulCacheInit(LeftMulCache *cache, int count) {
    const size_t budgetBytes = (size_t)96 * (size_t)1024 * (size_t)1024;
    size_t bytesPerSlot = (size_t)count * (sizeof(Poly57) + sizeof(int) + sizeof(int));
    if (bytesPerSlot == 0) {
        cache->slotCount = 0;
        return 1;
    }

    int slotCount = (int)(budgetBytes / bytesPerSlot);
    if (slotCount < 1) {
        cache->slotCount = 0;
        return 1;
    }
    if (slotCount > 64) {
        slotCount = 64;
    }

    cache->slotCount = slotCount;
    cache->leftIndex = (int *)calloc((size_t)slotCount, sizeof(int));
    cache->rightId = (int *)calloc((size_t)slotCount, sizeof(int));
    cache->age = (unsigned int *)calloc((size_t)slotCount, sizeof(unsigned int));
    cache->value = (SparseHecke *)calloc((size_t)slotCount, sizeof(SparseHecke));
    cache->tick = 1;

    if (!cache->leftIndex || !cache->rightId || !cache->age || !cache->value) {
        free(cache->leftIndex);
        free(cache->rightId);
        free(cache->age);
        free(cache->value);
        memset(cache, 0, sizeof(*cache));
        return 0;
    }

    for (int s = 0; s < slotCount; s++) {
        cache->leftIndex[s] = -1;
        cache->rightId[s] = -1;
        if (!SparseInit(&cache->value[s], count)) {
            for (int j = 0; j < s; j++) {
                SparseFree(&cache->value[j]);
            }
            free(cache->leftIndex);
            free(cache->rightId);
            free(cache->age);
            free(cache->value);
            memset(cache, 0, sizeof(*cache));
            return 0;
        }
    }

    return 1;
}

static void LeftMulCacheFree(LeftMulCache *cache) {
    for (int s = 0; s < cache->slotCount; s++) {
        SparseFree(&cache->value[s]);
    }
    free(cache->leftIndex);
    free(cache->rightId);
    free(cache->age);
    free(cache->value);
    memset(cache, 0, sizeof(*cache));
}

static const SparseHecke *LeftMulGetCached(
    const FastContext *ctx,
    int rightId,
    int leftIndex,
    const SparseHecke *right,
    SparseHecke *tmp1,
    SparseHecke *tmp2,
    LeftMulCache *cache
) {
    if (cache->slotCount <= 0) {
        LeftMultiplyByIndex(ctx, leftIndex, right, tmp1, tmp2);
        return tmp1;
    }

    int hit = -1;
    for (int s = 0; s < cache->slotCount; s++) {
        if (cache->leftIndex[s] == leftIndex && cache->rightId[s] == rightId) {
            hit = s;
            break;
        }
    }

    cache->tick++;
    if (cache->tick == 0) {
        cache->tick = 1;
    }

    if (hit >= 0) {
        cache->age[hit] = cache->tick;
        return &cache->value[hit];
    }

    int victim = -1;
    unsigned int oldest = 0u;
    for (int s = 0; s < cache->slotCount; s++) {
        if (cache->rightId[s] < 0) {
            victim = s;
            break;
        }
        if (victim < 0 || cache->age[s] < oldest) {
            victim = s;
            oldest = cache->age[s];
        }
    }

    LeftMultiplyByIndex(ctx, leftIndex, right, tmp1, tmp2);
    SparseAssign(&cache->value[victim], tmp1);
    cache->leftIndex[victim] = leftIndex;
    cache->rightId[victim] = rightId;
    cache->age[victim] = cache->tick;
    return &cache->value[victim];
}

static void SparseAssign(SparseHecke *dst, const SparseHecke *src) {
    SparseClear(dst);
    for (int i = 0; i < src->supportSize; i++) {
        int idx = src->support[i];
        int pos = dst->supportSize;
        dst->support[pos] = idx;
        dst->posMap[idx] = pos;
        dst->supportSize++;
        PolyCopy(&dst->coeff[idx], &src->coeff[idx]);
    }
}

static void SparseAssignFromCacheRow(
    SparseHecke *dst,
    const Poly57 *row,
    const int *support,
    int supportSize
) {
    SparseClear(dst);
    for (int i = 0; i < supportSize; i++) {
        int idx = support[i];
        int pos = dst->supportSize;
        dst->support[pos] = idx;
        dst->posMap[idx] = pos;
        dst->supportSize++;
        PolyCopy(&dst->coeff[idx], &row[idx]);
    }
}

static int64_t PairKey(int count, int y, int w) {
    return (int64_t)y * (int64_t)count + (int64_t)w;
}

static int CompareEntries(const void *a, const void *b) {
    const NonTrivialEntry *ea = (const NonTrivialEntry *)a;
    const NonTrivialEntry *eb = (const NonTrivialEntry *)b;
    if (ea->key < eb->key) return -1;
    if (ea->key > eb->key) return 1;
    return 0;
}

static const Poly57 *FindNonTrivial(const FastContext *ctx, int y, int w) {
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

static int ParsePermutationString(int n, const char *token, int *outIndex) {
    if ((int)strlen(token) != n || n > 8) {
        return 0;
    }

    int seen[8] = {0};
    char perm[8] = {0};

    for (int i = 0; i < n; i++) {
        if (!isdigit((unsigned char)token[i])) {
            return 0;
        }
        int d = token[i] - '0';
        if (d < 0 || d >= n || seen[d]) {
            return 0;
        }
        seen[d] = 1;
        perm[i] = (char)d;
    }

    *outIndex = PermToIndex(n, perm);
    return 1;
}

static Poly57 KLCoefficient(const FastContext *ctx, int y, int w) {
    Poly57 p;
    PolyZero(&p);

    if (y == w) {
        PolyFromMonomial(&p, 0, 1);
        return p;
    }

    if (!BruhatSmaller2(ctx->n, y, w)) {
        return p;
    }

    const Poly57 *nt = FindNonTrivial(ctx, y, w);
    if (nt) {
        return *nt;
    }

    int diff = ctx->lengths[w] - ctx->lengths[y];
    PolyFromMonomial(&p, diff, 1);
    return p;
}

static int LoadNonTrivialData(FastContext *ctx, const char *filename) {
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
        fprintf(stderr, "Out of memory loading nontrivial data.\n");
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
            int pos = 28 + exponent;
            if (coeff != 0 && pos >= 0 && pos < 57) {
                p.coeff[pos] += coeff;
                if (p.minPos > p.maxPos) {
                    p.minPos = pos;
                    p.maxPos = pos;
                } else {
                    if (pos < p.minPos) p.minPos = pos;
                    if (pos > p.maxPos) p.maxPos = pos;
                }
            }
            part = strtok(NULL, ",");
            k++;
        }
        PolyNormalize(&p);

        if (used == cap) {
            size_t newCap = cap * 2;
            NonTrivialEntry *newEntries = (NonTrivialEntry *)realloc(entries, newCap * sizeof(NonTrivialEntry));
            if (!newEntries) {
                free(entries);
                fclose(f);
                fprintf(stderr, "Out of memory while growing nontrivial table.\n");
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

static int CompareByLengthDesc(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    int la = g_sortLengths[ia];
    int lb = g_sortLengths[ib];
    if (la != lb) {
        return lb - la;
    }
    return ia - ib;
}

static int BuildLengthOrder(FastContext *ctx) {
    ctx->indicesByLengthDesc = (int *)calloc((size_t)ctx->count, sizeof(int));
    if (!ctx->indicesByLengthDesc) {
        return 0;
    }

    for (int i = 0; i < ctx->count; i++) {
        ctx->indicesByLengthDesc[i] = i;
    }

    g_sortLengths = ctx->lengths;
    qsort(ctx->indicesByLengthDesc, (size_t)ctx->count, sizeof(int), CompareByLengthDesc);
    g_sortLengths = NULL;
    return 1;
}

static int BuildGeneratorTables(FastContext *ctx) {
    ctx->genCount = ctx->n - 1;
    for (int i = 0; i < ctx->genCount; i++) {
        int gen = i + 1;
        ctx->genLehmer[i] = fac(ctx->n - gen);
    }

    size_t actionSize = (size_t)ctx->genCount * (size_t)ctx->count;
    ctx->leftAction = (int *)calloc(actionSize, sizeof(int));
    ctx->lengthIncrease = (unsigned char *)calloc(actionSize, sizeof(unsigned char));
    if (!ctx->leftAction || !ctx->lengthIncrease) {
        return 0;
    }

    for (int g = 0; g < ctx->genCount; g++) {
        int s = ctx->genLehmer[g];
        for (int w = 0; w < ctx->count; w++) {
            int comp = MultiplyIndex(ctx->n, s, w);
            size_t pos = (size_t)g * (size_t)ctx->count + (size_t)w;
            ctx->leftAction[pos] = comp;
            ctx->lengthIncrease[pos] = (ctx->lengths[comp] == ctx->lengths[w] + 1) ? 1 : 0;
        }
    }

    return 1;
}

static int BuildExpressionTables(FastContext *ctx) {
    ctx->exprLen = (int *)calloc((size_t)ctx->count, sizeof(int));
    ctx->exprG = (unsigned char *)calloc((size_t)ctx->count * (size_t)ctx->maxLength, sizeof(unsigned char));
    if (!ctx->exprLen || !ctx->exprG) {
        return 0;
    }

    for (int w = 0; w < ctx->count; w++) {
        char expr[40] = {0};
        int len = ReducedExpression(ctx->n, w, expr);
        ctx->exprLen[w] = len;
        for (int i = 0; i < len; i++) {
            int gen = (int)expr[i];
            int g = gen - 1;
            if (g < 0) g = 0;
            if (g >= ctx->genCount) g = ctx->genCount - 1;
            ctx->exprG[(size_t)w * (size_t)ctx->maxLength + (size_t)i] = (unsigned char)g;
        }
    }

    return 1;
}

static int MaybeBuildKLBasisCache(FastContext *ctx) {
    const size_t cacheLimitBytes = (size_t)512 * (size_t)1024 * (size_t)1024;
    size_t needed = (size_t)ctx->count * (size_t)ctx->count * sizeof(Poly57);

    if (needed > cacheLimitBytes) {
        ctx->useCache = 0;
        return 1;
    }

    ctx->klCache = (Poly57 *)calloc((size_t)ctx->count * (size_t)ctx->count, sizeof(Poly57));
    ctx->klSupportOffsets = (int *)calloc((size_t)ctx->count + 1, sizeof(int));
    ctx->klSupportIndices = (int *)calloc((size_t)ctx->count * (size_t)ctx->count, sizeof(int));
    if (!ctx->klCache || !ctx->klSupportOffsets || !ctx->klSupportIndices) {
        free(ctx->klCache);
        free(ctx->klSupportOffsets);
        free(ctx->klSupportIndices);
        ctx->klCache = NULL;
        ctx->klSupportOffsets = NULL;
        ctx->klSupportIndices = NULL;
        ctx->useCache = 0;
        return 1;
    }

    int supportUsed = 0;
    for (int w = 0; w < ctx->count; w++) {
        ctx->klSupportOffsets[w] = supportUsed;
        Poly57 *row = ctx->klCache + (size_t)w * (size_t)ctx->count;
        for (int y = 0; y < ctx->count; y++) {
            row[y] = KLCoefficient(ctx, y, w);
            if (!PolyIsZero(&row[y])) {
                ctx->klSupportIndices[supportUsed++] = y;
            }
        }
    }
    ctx->klSupportOffsets[ctx->count] = supportUsed;

    int *shrink = (int *)realloc(ctx->klSupportIndices, (size_t)supportUsed * sizeof(int));
    if (shrink) {
        ctx->klSupportIndices = shrink;
    }

    ctx->useCache = 1;
    return 1;
}

static void BuildKLElementSparse(const FastContext *ctx, int w, SparseHecke *out) {
    SparseClear(out);
    for (int y = 0; y < ctx->count; y++) {
        Poly57 p = KLCoefficient(ctx, y, w);
        if (!PolyIsZero(&p)) {
            SparseAddPoly(out, y, &p);
        }
    }
}

static void MultiplySimpleGenerator(
    const FastContext *ctx,
    int g,
    const SparseHecke *in,
    SparseHecke *out
) {
    SparseClear(out);

    const int *action = ctx->leftAction + (size_t)g * (size_t)ctx->count;
    const unsigned char *inc = ctx->lengthIncrease + (size_t)g * (size_t)ctx->count;

    for (int p = 0; p < in->supportSize; p++) {
        int idx = in->support[p];
        const Poly57 *poly = &in->coeff[idx];

        int comp = action[idx];
        SparseAddPoly(out, comp, poly);

        if (!inc[idx]) {
            Poly57 corr;
            PolyMultiplyLowerTerm(poly, &corr);
            SparseAddPoly(out, idx, &corr);
        }
    }
}

static void LeftMultiplyByIndex(
    const FastContext *ctx,
    int index,
    const SparseHecke *base,
    SparseHecke *result,
    SparseHecke *scratch
) {
    SparseAssign(result, base);

    int len = ctx->exprLen[index];
    SparseHecke *a = result;
    SparseHecke *b = scratch;

    for (int step = 0; step < len; step++) {
        int g = (int)ctx->exprG[(size_t)index * (size_t)ctx->maxLength + (size_t)step];
        MultiplySimpleGenerator(ctx, g, a, b);

        SparseHecke *tmp = a;
        a = b;
        b = tmp;
    }

    if (a != result) {
        SparseAssign(result, a);
    }
}

static void MultiplyHeckeSparse(
    const FastContext *ctx,
    int rightId,
    const SparseHecke *left,
    const SparseHecke *right,
    SparseHecke *product,
    SparseHecke *tmp1,
    SparseHecke *tmp2,
    DenseAccum *accum,
    LeftMulCache *leftMulCache
) {
    DenseAccumBegin(accum);

    for (int lp = 0; lp < left->supportSize; lp++) {
        int i = left->support[lp];
        const Poly57 *scalar = &left->coeff[i];

        const SparseHecke *lm = LeftMulGetCached(
            ctx,
            rightId,
            i,
            right,
            tmp1,
            tmp2,
            leftMulCache
        );

        for (int rp = 0; rp < lm->supportSize; rp++) {
            int j = lm->support[rp];
            Poly57 scaled;
            PolyMultiply(scalar, &lm->coeff[j], &scaled);
            DenseAccumAdd(accum, j, &scaled);
        }
    }

    DenseAccumCommitToSparse(accum, product);
}

static void ClearCoeffArray(Poly57 *coeffs, int count) {
    for (int i = 0; i < count; i++) {
        PolyZero(&coeffs[i]);
    }
}

static void DecomposeToKLBasis(
    const FastContext *ctx,
    const SparseHecke *element,
    Poly57 *coeffs,
    SparseHecke *residual,
    SparseHecke *basisTmp
) {
    ClearCoeffArray(coeffs, ctx->count);
    SparseAssign(residual, element);

    for (int ord = 0; ord < ctx->count; ord++) {
        int z = ctx->indicesByLengthDesc[ord];
        if (residual->posMap[z] < 0) {
            continue;
        }

        Poly57 top;
        PolyCopy(&top, &residual->coeff[z]);
        if (PolyIsZero(&top)) {
            continue;
        }

        coeffs[z] = top;

        Poly57 minusTop;
        PolyNegate(&minusTop, &top);

        if (ctx->useCache) {
            const Poly57 *row = ctx->klCache + (size_t)z * (size_t)ctx->count;
            int start = ctx->klSupportOffsets[z];
            int stop = ctx->klSupportOffsets[z + 1];
            for (int p = start; p < stop; p++) {
                int y = ctx->klSupportIndices[p];
                Poly57 corr;
                PolyMultiply(&minusTop, &row[y], &corr);
                SparseAddPoly(residual, y, &corr);
            }
        } else {
            BuildKLElementSparse(ctx, z, basisTmp);
            for (int p = 0; p < basisTmp->supportSize; p++) {
                int y = basisTmp->support[p];
                Poly57 corr;
                PolyMultiply(&minusTop, &basisTmp->coeff[y], &corr);
                SparseAddPoly(residual, y, &corr);
            }
        }
    }
}

static void WriteTextProduct(
    FILE *out,
    const FastContext *ctx,
    int x,
    int y,
    const char *permTable,
    const Poly57 *coeffs
) {
    const char *px = permTable + (size_t)x * (size_t)(ctx->n + 1);
    const char *py = permTable + (size_t)y * (size_t)(ctx->n + 1);

    fprintf(out, "C_{%s} * C_{%s} = ", px, py);

    int first = 1;
    for (int z = 0; z < ctx->count; z++) {
        if (PolyIsZero(&coeffs[z])) {
            continue;
        }

        if (!first) {
            fprintf(out, " + ");
        }

        fprintf(out, "(");
        FprintPoly(out, &coeffs[z]);
        fprintf(out, ")C_{%s}", permTable + (size_t)z * (size_t)(ctx->n + 1));
        first = 0;
    }

    if (first) {
        fprintf(out, "0");
    }
    fprintf(out, "\n");
}

static void WriteCsvTerms(
    FILE *out,
    const FastContext *ctx,
    int x,
    int y,
    const char *permTable,
    const Poly57 *coeffs
) {
    const char *px = permTable + (size_t)x * (size_t)(ctx->n + 1);
    const char *py = permTable + (size_t)y * (size_t)(ctx->n + 1);

    int wrote = 0;
    for (int z = 0; z < ctx->count; z++) {
        if (PolyIsZero(&coeffs[z])) {
            continue;
        }

        const char *pz = permTable + (size_t)z * (size_t)(ctx->n + 1);
        fprintf(out, "\"%s\",\"%s\",\"%s\",\"", px, py, pz);
        FprintPoly(out, &coeffs[z]);
        fprintf(out, "\"\n");
        wrote = 1;
    }

    if (!wrote) {
        fprintf(out, "\"%s\",\"%s\",\"\",\"0\"\n", px, py);
    }
}

static char *BuildPermStringTable(int n, int count) {
    char *table = (char *)calloc((size_t)count * (size_t)(n + 1), sizeof(char));
    if (!table) {
        return NULL;
    }

    for (int i = 0; i < count; i++) {
        char perm[8] = {0};
        IndexToPerm(n, i, perm);
        char *dst = table + (size_t)i * (size_t)(n + 1);
        for (int j = 0; j < n; j++) {
            dst[j] = (char)('0' + perm[j]);
        }
        dst[n] = '\0';
    }

    return table;
}

static void FreeContext(FastContext *ctx) {
    free(ctx->lengths);
    free(ctx->indicesByLengthDesc);
    free(ctx->leftAction);
    free(ctx->lengthIncrease);
    free(ctx->exprLen);
    free(ctx->exprG);
    free(ctx->entries);
    free(ctx->klCache);
    free(ctx->klSupportOffsets);
    free(ctx->klSupportIndices);
}

static void PrintUsage(const char *prog) {
    printf("Usage:\n");
    printf("  %s n [output_txt] [output_csv] [--bench]\n", prog);
    printf("\n");
    printf("n must be in {4,5,6,7,8,9}.\n");
    printf("If outputs are omitted, defaults are S<n>_structure_constants_bulk_opt.txt/.csv.\n");
    printf("Use --bench to print per-phase wall/CPU timings to stderr.\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2 || argc > 5) {
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
    snprintf(defaultTxt, sizeof(defaultTxt), "S%d_structure_constants_bulk_opt.txt", n);
    snprintf(defaultCsv, sizeof(defaultCsv), "S%d_structure_constants_bulk_opt.csv", n);

    const char *txtPath = defaultTxt;
    const char *csvPath = defaultCsv;
    int benchEnabled = 0;

    int outputArg = 0;
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--bench") == 0) {
            benchEnabled = 1;
            continue;
        }

        if (outputArg == 0) {
            txtPath = argv[i];
            outputArg = 1;
        } else if (outputArg == 1) {
            csvPath = argv[i];
            outputArg = 2;
        } else {
            fprintf(stderr, "Too many positional arguments.\n");
            PrintUsage(argv[0]);
            return 1;
        }
    }

    FastContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = n;
    ctx.count = fac(n);
    ctx.maxLength = n * (n - 1) / 2;
    ctx.useCache = 0;

    PolyZero(&ctx.lowerTerm);
    ctx.lowerTerm.coeff[27] = 1;
    ctx.lowerTerm.coeff[29] = -1;
    ctx.lowerTerm.minPos = 27;
    ctx.lowerTerm.maxPos = 29;

    ctx.lengths = (int *)calloc((size_t)ctx.count, sizeof(int));
    if (!ctx.lengths) {
        fprintf(stderr, "Out of memory for length table.\n");
        return 1;
    }
    for (int i = 0; i < ctx.count; i++) {
        ctx.lengths[i] = IndexToLength(n, i);
    }

    if (!BuildLengthOrder(&ctx) || !BuildGeneratorTables(&ctx) || !BuildExpressionTables(&ctx)) {
        fprintf(stderr, "Failed to build precomputed tables.\n");
        FreeContext(&ctx);
        return 1;
    }

    if (!LoadNonTrivialData(&ctx, dataFile)) {
        FreeContext(&ctx);
        return 1;
    }

    if (!MaybeBuildKLBasisCache(&ctx)) {
        fprintf(stderr, "Failed while preparing KL cache.\n");
        FreeContext(&ctx);
        return 1;
    }

    char *permTable = BuildPermStringTable(n, ctx.count);
    if (!permTable) {
        fprintf(stderr, "Out of memory for permutation string table.\n");
        FreeContext(&ctx);
        return 1;
    }

    FILE *txtOut = fopen(txtPath, "w");
    FILE *csvOut = fopen(csvPath, "w");
    if (!txtOut || !csvOut) {
        fprintf(stderr, "Could not open output files.\n");
        if (txtOut) fclose(txtOut);
        if (csvOut) fclose(csvOut);
        free(permTable);
        FreeContext(&ctx);
        return 1;
    }

    /* Use larger buffers so formatted output causes fewer flushes/syscalls. */
    setvbuf(txtOut, NULL, _IOFBF, (size_t)1 << 20);
    setvbuf(csvOut, NULL, _IOFBF, (size_t)1 << 20);

    fprintf(txtOut, "# All KL-basis products in S%d (optimized bulk)\n", n);
    fprintf(txtOut, "# Format: C_{x} * C_{y} = sum_z h_{x,y}^z(v) C_{z}\n\n");
    fprintf(csvOut, "\"x\",\"y\",\"z\",\"h\"\n");

    SparseHecke Cx, Cy, product, tmp1, tmp2, residual, basisTmp;
    if (!SparseInit(&Cx, ctx.count) || !SparseInit(&Cy, ctx.count) || !SparseInit(&product, ctx.count) ||
        !SparseInit(&tmp1, ctx.count) || !SparseInit(&tmp2, ctx.count) || !SparseInit(&residual, ctx.count) ||
        !SparseInit(&basisTmp, ctx.count)) {
        fprintf(stderr, "Out of memory creating sparse workspaces.\n");
        SparseFree(&Cx);
        SparseFree(&Cy);
        SparseFree(&product);
        SparseFree(&tmp1);
        SparseFree(&tmp2);
        SparseFree(&residual);
        SparseFree(&basisTmp);
        fclose(txtOut);
        fclose(csvOut);
        free(permTable);
        FreeContext(&ctx);
        return 1;
    }

    DenseAccum accum;
    memset(&accum, 0, sizeof(accum));
    if (!DenseAccumInit(&accum, ctx.count)) {
        fprintf(stderr, "Out of memory for dense accumulation workspace.\n");
        SparseFree(&Cx);
        SparseFree(&Cy);
        SparseFree(&product);
        SparseFree(&tmp1);
        SparseFree(&tmp2);
        SparseFree(&residual);
        SparseFree(&basisTmp);
        fclose(txtOut);
        fclose(csvOut);
        free(permTable);
        FreeContext(&ctx);
        return 1;
    }

    LeftMulCache leftMulCache;
    memset(&leftMulCache, 0, sizeof(leftMulCache));
    if (!LeftMulCacheInit(&leftMulCache, ctx.count)) {
        fprintf(stderr, "Out of memory for left-multiply cache.\n");
        DenseAccumFree(&accum);
        SparseFree(&Cx);
        SparseFree(&Cy);
        SparseFree(&product);
        SparseFree(&tmp1);
        SparseFree(&tmp2);
        SparseFree(&residual);
        SparseFree(&basisTmp);
        fclose(txtOut);
        fclose(csvOut);
        free(permTable);
        FreeContext(&ctx);
        return 1;
    }

    Poly57 *coeffs = (Poly57 *)calloc((size_t)ctx.count, sizeof(Poly57));
    if (!coeffs) {
        fprintf(stderr, "Out of memory for decomposition coefficients.\n");
        SparseFree(&Cx);
        SparseFree(&Cy);
        SparseFree(&product);
        SparseFree(&tmp1);
        SparseFree(&tmp2);
        SparseFree(&residual);
        SparseFree(&basisTmp);
        fclose(txtOut);
        fclose(csvOut);
        free(permTable);
        FreeContext(&ctx);
        return 1;
    }

    int64_t totalPairs = (int64_t)ctx.count * (int64_t)ctx.count;
    int64_t donePairs = 0;
    double startCpu = CpuNowSeconds();
    double startWall = WallNowSeconds();
    BenchStats bench;
    memset(&bench, 0, sizeof(bench));

    for (int x = 0; x < ctx.count; x++) {
        double t0w = 0.0, t0c = 0.0;
        if (benchEnabled) {
            t0w = WallNowSeconds();
            t0c = CpuNowSeconds();
        }

        if (ctx.useCache) {
            const Poly57 *rowX = ctx.klCache + (size_t)x * (size_t)ctx.count;
            int sx = ctx.klSupportOffsets[x];
            int ex = ctx.klSupportOffsets[x + 1];
            SparseAssignFromCacheRow(&Cx, rowX, ctx.klSupportIndices + sx, ex - sx);
        } else {
            BuildKLElementSparse(&ctx, x, &Cx);
        }

        if (benchEnabled) {
            bench.buildXWall += WallNowSeconds() - t0w;
            bench.buildXCpu += CpuNowSeconds() - t0c;
        }

        for (int y = 0; y < ctx.count; y++) {
            if (benchEnabled) {
                t0w = WallNowSeconds();
                t0c = CpuNowSeconds();
            }

            if (ctx.useCache) {
                const Poly57 *rowY = ctx.klCache + (size_t)y * (size_t)ctx.count;
                int sy = ctx.klSupportOffsets[y];
                int ey = ctx.klSupportOffsets[y + 1];
                SparseAssignFromCacheRow(&Cy, rowY, ctx.klSupportIndices + sy, ey - sy);
            } else {
                BuildKLElementSparse(&ctx, y, &Cy);
            }

            if (benchEnabled) {
                bench.buildYWall += WallNowSeconds() - t0w;
                bench.buildYCpu += CpuNowSeconds() - t0c;
                t0w = WallNowSeconds();
                t0c = CpuNowSeconds();
            }

            MultiplyHeckeSparse(
                &ctx,
                y,
                &Cx,
                &Cy,
                &product,
                &tmp1,
                &tmp2,
                &accum,
                &leftMulCache
            );

            if (benchEnabled) {
                bench.multiplyWall += WallNowSeconds() - t0w;
                bench.multiplyCpu += CpuNowSeconds() - t0c;
                t0w = WallNowSeconds();
                t0c = CpuNowSeconds();
            }

            DecomposeToKLBasis(&ctx, &product, coeffs, &residual, &basisTmp);

            if (benchEnabled) {
                bench.decomposeWall += WallNowSeconds() - t0w;
                bench.decomposeCpu += CpuNowSeconds() - t0c;
                t0w = WallNowSeconds();
                t0c = CpuNowSeconds();
            }

            WriteTextProduct(txtOut, &ctx, x, y, permTable, coeffs);

            if (benchEnabled) {
                bench.writeTxtWall += WallNowSeconds() - t0w;
                bench.writeTxtCpu += CpuNowSeconds() - t0c;
                t0w = WallNowSeconds();
                t0c = CpuNowSeconds();
            }

            WriteCsvTerms(csvOut, &ctx, x, y, permTable, coeffs);

            if (benchEnabled) {
                bench.writeCsvWall += WallNowSeconds() - t0w;
                bench.writeCsvCpu += CpuNowSeconds() - t0c;
            }

            donePairs++;
        }

        if ((x + 1) % 5 == 0 || x + 1 == ctx.count) {
            double elapsed = WallNowSeconds() - startWall;
            fprintf(stderr, "Progress: %d/%d x-rows, %.2f%% pairs, elapsed %.1fs\n",
                x + 1,
                ctx.count,
                100.0 * (double)donePairs / (double)totalPairs,
                elapsed
            );
        }
    }

    double elapsedWall = WallNowSeconds() - startWall;
    double elapsedCpu = CpuNowSeconds() - startCpu;
    fprintf(stderr, "Done. Wrote %s and %s in %.1fs wall (%.1fs CPU)\n", txtPath, csvPath, elapsedWall, elapsedCpu);

    if (benchEnabled) {
        PrintBenchReport(&ctx, &bench, elapsedWall, elapsedCpu, donePairs);
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
    fclose(txtOut);
    fclose(csvOut);
    free(permTable);
    FreeContext(&ctx);

    return 0;
}
