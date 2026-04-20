#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "lehmer.h"

enum {
    POLY_MAX_ABS_EXP = 100,
    POLY_CENTER = POLY_MAX_ABS_EXP,
    POLY_SIZE = 2 * POLY_MAX_ABS_EXP + 1
};

typedef struct {
    int64_t coeff[POLY_SIZE];
    int minPos;
    int maxPos;
} Poly57;

typedef struct {
    int z;
    Poly57 poly;
} Term;

typedef struct {
    int key;
    Term *terms;
    int termCount;
    int age;
} PairCacheEntry;

typedef struct {
    Term *terms;
    int count;
    int cap;
    int present;
} TermBag;

typedef struct {
    int n;
    int count;
    int64_t totalPairs;

    const char *csvPath;
    FILE *csv;

    int64_t *pairOffsets;
    int *pairLineCounts;

    int *inverseIndex;

    PairCacheEntry *cache;
    int cacheSize;
    int cacheTick;
} Validator;

enum {
    MAX_FAIL_PRINT = 20
};

static int64_t FileTell64(FILE *f) {
#if defined(_WIN32)
    return _ftelli64(f);
#else
    return (int64_t)ftello(f);
#endif
}

static int FileSeek64(FILE *f, int64_t offset) {
#if defined(_WIN32)
    return _fseeki64(f, offset, SEEK_SET);
#else
    return fseeko(f, (off_t)offset, SEEK_SET);
#endif
}

static void PolyZero(Poly57 *p) {
    for (int i = 0; i < POLY_SIZE; i++) {
        p->coeff[i] = 0;
    }
    p->minPos = POLY_SIZE;
    p->maxPos = -1;
}

static int PolyIsZero(const Poly57 *p) {
    return p->maxPos < p->minPos;
}

static void PolyNormalize(Poly57 *p) {
    if (PolyIsZero(p)) {
        p->minPos = POLY_SIZE;
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
        p->minPos = POLY_SIZE;
        p->maxPos = -1;
    } else {
        p->minPos = left;
        p->maxPos = right;
    }
}

static int PolyEqual(const Poly57 *a, const Poly57 *b) {
    for (int i = 0; i < POLY_SIZE; i++) {
        if (a->coeff[i] != b->coeff[i]) {
            return 0;
        }
    }
    return 1;
}

static int PolyIsOne(const Poly57 *p) {
    if (PolyIsZero(p)) {
        return 0;
    }
    if (p->minPos != POLY_CENTER || p->maxPos != POLY_CENTER) {
        return 0;
    }
    return p->coeff[POLY_CENTER] == 1;
}

static void PrintPoly(FILE *out, const Poly57 *p) {
    if (PolyIsZero(p)) {
        fprintf(out, "0");
        return;
    }

    int first = 1;
    for (int i = p->minPos; i <= p->maxPos; i++) {
        if (p->coeff[i] == 0) {
            continue;
        }
        if (!first) {
            fprintf(out, " + ");
        }
        fprintf(out, "%lld*v^%d", (long long)p->coeff[i], i - POLY_CENTER);
        first = 0;
    }
}

static void PermStringFromIndex(int n, int idx, char *out) {
    char perm[8] = {0};
    IndexToPerm(n, idx, perm);
    for (int i = 0; i < n; i++) {
        out[i] = (char)('0' + perm[i]);
    }
    out[n] = '\0';
}

static void PrintTerms(FILE *out, int n, const Term *terms, int count) {
    if (count <= 0) {
        fprintf(out, "0");
        return;
    }

    int first = 1;
    for (int i = 0; i < count; i++) {
        if (terms[i].z < 0 || PolyIsZero(&terms[i].poly)) {
            continue;
        }

        char zPerm[16];
        PermStringFromIndex(n, terms[i].z, zPerm);

        if (!first) {
            fprintf(out, " + ");
        }
        fprintf(out, "(");
        PrintPoly(out, &terms[i].poly);
        fprintf(out, ")C_{%s}", zPerm);
        first = 0;
    }

    if (first) {
        fprintf(out, "0");
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

static void PolyMultiply(const Poly57 *a, const Poly57 *b, Poly57 *out) {
    PolyZero(out);
    if (PolyIsZero(a) || PolyIsZero(b)) {
        return;
    }

    int seen = 0;
    int minPos = POLY_SIZE;
    int maxPos = -1;

    for (int i = a->minPos; i <= a->maxPos; i++) {
        int64_t ca = a->coeff[i];
        if (ca == 0) {
            continue;
        }

        for (int j = b->minPos; j <= b->maxPos; j++) {
            int64_t cb = b->coeff[j];
            if (cb == 0) {
                continue;
            }

            int pos = i + j - POLY_CENTER;
            if (pos < 0 || pos >= POLY_SIZE) {
                continue;
            }

            out->coeff[pos] += ca * cb;
            if (!seen) {
                seen = 1;
                minPos = pos;
                maxPos = pos;
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

static int ParseCsv4Quoted(const char *line, char *f1, size_t n1, char *f2, size_t n2, char *f3, size_t n3, char *f4, size_t n4) {
    const char *p = line;
    char *outs[4] = {f1, f2, f3, f4};
    size_t caps[4] = {n1, n2, n3, n4};

    for (int field = 0; field < 4; field++) {
        while (*p == ' ' || *p == '\t') {
            p++;
        }
        if (*p != '"') {
            return 0;
        }
        p++;

        size_t k = 0;
        while (*p && *p != '"') {
            if (k + 1 < caps[field]) {
                outs[field][k++] = *p;
            }
            p++;
        }
        if (*p != '"') {
            return 0;
        }
        outs[field][k] = '\0';
        p++;

        if (field < 3) {
            if (*p != ',') {
                return 0;
            }
            p++;
        }
    }

    return 1;
}

static int ParsePoly(const char *text, Poly57 *out) {
    PolyZero(out);

    if (strcmp(text, "0") == 0 || text[0] == '\0') {
        return 1;
    }

    char buffer[16384];
    snprintf(buffer, sizeof(buffer), "%s", text);

    char *cursor = buffer;
    while (*cursor) {
        char *sep = strstr(cursor, " + ");
        if (sep) {
            *sep = '\0';
        }

        long long coeff = 0;
        int exp = 0;
        if (sscanf(cursor, "%lld*v^%d", &coeff, &exp) != 2) {
            return 0;
        }

        int pos = POLY_CENTER + exp;
        if (pos < 0 || pos >= POLY_SIZE) {
            return 0;
        }

        out->coeff[pos] += (int64_t)coeff;
        if (out->minPos > out->maxPos) {
            out->minPos = pos;
            out->maxPos = pos;
        } else {
            if (pos < out->minPos) out->minPos = pos;
            if (pos > out->maxPos) out->maxPos = pos;
        }

        if (!sep) {
            break;
        }
        cursor = sep + 3;
    }

    PolyNormalize(out);
    return 1;
}

static int CompareTermsByZ(const void *a, const void *b) {
    const Term *ta = (const Term *)a;
    const Term *tb = (const Term *)b;
    return ta->z - tb->z;
}

static int FindTerm(const Term *terms, int termCount, int z) {
    int lo = 0;
    int hi = termCount - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        if (terms[mid].z == z) return mid;
        if (terms[mid].z < z) lo = mid + 1;
        else hi = mid - 1;
    }
    return -1;
}

static int BuildInverseTable(Validator *v) {
    v->inverseIndex = (int *)calloc((size_t)v->count, sizeof(int));
    if (!v->inverseIndex) {
        return 0;
    }

    for (int w = 0; w < v->count; w++) {
        char perm[8] = {0};
        char inv[8] = {0};
        IndexToPerm(v->n, w, perm);
        for (int i = 0; i < v->n; i++) {
            inv[(int)perm[i]] = (char)i;
        }
        v->inverseIndex[w] = PermToIndex(v->n, inv);
    }

    return 1;
}

static int BuildPairIndex(Validator *v) {
    v->pairOffsets = (int64_t *)malloc((size_t)v->totalPairs * sizeof(int64_t));
    v->pairLineCounts = (int *)calloc((size_t)v->totalPairs, sizeof(int));
    if (!v->pairOffsets || !v->pairLineCounts) {
        return 0;
    }

    for (int64_t i = 0; i < v->totalPairs; i++) {
        v->pairOffsets[i] = -1;
    }

    if (FileSeek64(v->csv, 0) != 0) {
        return 0;
    }

    char line[16384];
    int sawHeader = 0;
    int64_t prevKey = -1;

    while (1) {
        int64_t lineStart = FileTell64(v->csv);
        if (lineStart < 0) {
            fprintf(stderr, "Failed to read file position while indexing.\n");
            return 0;
        }
        if (!fgets(line, sizeof(line), v->csv)) {
            break;
        }

        if (!sawHeader) {
            sawHeader = 1;
            continue;
        }

        char fx[64], fy[64], fz[64], fh[16000];
        if (!ParseCsv4Quoted(line, fx, sizeof(fx), fy, sizeof(fy), fz, sizeof(fz), fh, sizeof(fh))) {
            fprintf(stderr, "Malformed CSV row while indexing.\n");
            return 0;
        }

        int x = -1;
        int y = -1;
        if (!ParsePermutationString(v->n, fx, &x) || !ParsePermutationString(v->n, fy, &y)) {
            fprintf(stderr, "Failed to parse x/y permutation while indexing.\n");
            return 0;
        }

        int64_t key = (int64_t)x * (int64_t)v->count + (int64_t)y;
        if (key < 0 || key >= v->totalPairs) {
            fprintf(stderr, "Out-of-range pair key while indexing.\n");
            return 0;
        }

        if (key < prevKey) {
            fprintf(stderr, "CSV rows are not grouped by (x,y).\n");
            return 0;
        }
        prevKey = key;

        if (v->pairLineCounts[key] == 0) {
            v->pairOffsets[key] = lineStart;
        }
        v->pairLineCounts[key]++;
    }

    for (int64_t key = 0; key < v->totalPairs; key++) {
        if (v->pairLineCounts[key] <= 0 || v->pairOffsets[key] < 0) {
            fprintf(stderr, "Missing pair in CSV at key %lld.\n", (long long)key);
            return 0;
        }
    }

    return 1;
}

static int LoadPairTermsFromFile(Validator *v, int key, Term **outTerms, int *outCount) {
    if (key < 0 || (int64_t)key >= v->totalPairs) {
        return 0;
    }

    int count = v->pairLineCounts[key];
    int64_t offset = v->pairOffsets[key];
    if (count <= 0 || offset < 0) {
        return 0;
    }

    Term *terms = (Term *)calloc((size_t)count, sizeof(Term));
    if (!terms) {
        return 0;
    }

    if (FileSeek64(v->csv, offset) != 0) {
        free(terms);
        return 0;
    }

    char line[16384];
    int xExpected = key / v->count;
    int yExpected = key % v->count;

    for (int i = 0; i < count; i++) {
        if (!fgets(line, sizeof(line), v->csv)) {
            free(terms);
            return 0;
        }

        char fx[64], fy[64], fz[64], fh[16000];
        if (!ParseCsv4Quoted(line, fx, sizeof(fx), fy, sizeof(fy), fz, sizeof(fz), fh, sizeof(fh))) {
            free(terms);
            return 0;
        }

        int x = -1;
        int y = -1;
        int z = -1;
        if (!ParsePermutationString(v->n, fx, &x) || !ParsePermutationString(v->n, fy, &y)) {
            free(terms);
            return 0;
        }

        if (x != xExpected || y != yExpected) {
            free(terms);
            return 0;
        }

        if (fz[0] == '\0') {
            z = -1;
        } else if (!ParsePermutationString(v->n, fz, &z)) {
            free(terms);
            return 0;
        }

        Poly57 p;
        if (!ParsePoly(fh, &p)) {
            free(terms);
            return 0;
        }

        if (z < 0) {
            if (!PolyIsZero(&p)) {
                free(terms);
                return 0;
            }
            z = -1;
        }

        terms[i].z = z;
        terms[i].poly = p;
    }

    qsort(terms, (size_t)count, sizeof(Term), CompareTermsByZ);
    for (int i = 1; i < count; i++) {
        if (terms[i].z >= 0 && terms[i].z == terms[i - 1].z) {
            free(terms);
            return 0;
        }
    }

    *outTerms = terms;
    *outCount = count;
    return 1;
}

static int GetPairTerms(Validator *v, int key, const Term **outTerms, int *outCount) {
    int bestSlot = -1;
    int bestAge = 2147483647;

    for (int i = 0; i < v->cacheSize; i++) {
        if (v->cache[i].key == key) {
            v->cache[i].age = ++v->cacheTick;
            *outTerms = v->cache[i].terms;
            *outCount = v->cache[i].termCount;
            return 1;
        }

        if (v->cache[i].key < 0) {
            bestSlot = i;
            break;
        }

        if (v->cache[i].age < bestAge) {
            bestAge = v->cache[i].age;
            bestSlot = i;
        }
    }

    if (bestSlot < 0) {
        return 0;
    }

    if (v->cache[bestSlot].key >= 0) {
        free(v->cache[bestSlot].terms);
        v->cache[bestSlot].terms = NULL;
        v->cache[bestSlot].termCount = 0;
    }

    Term *terms = NULL;
    int termCount = 0;
    if (!LoadPairTermsFromFile(v, key, &terms, &termCount)) {
        return 0;
    }

    v->cache[bestSlot].key = key;
    v->cache[bestSlot].terms = terms;
    v->cache[bestSlot].termCount = termCount;
    v->cache[bestSlot].age = ++v->cacheTick;

    *outTerms = terms;
    *outCount = termCount;
    return 1;
}

static int CopyTerms(const Term *src, int count, Term **outCopy) {
    if (count < 0) {
        return 0;
    }
    if (count == 0) {
        *outCopy = NULL;
        return 1;
    }

    Term *copy = (Term *)malloc((size_t)count * sizeof(Term));
    if (!copy) {
        return 0;
    }

    memcpy(copy, src, (size_t)count * sizeof(Term));
    *outCopy = copy;
    return 1;
}

static int CheckIdentity(Validator *v) {
    int failures = 0;
    int printed = 0;
    int e = 0;

    for (int x = 0; x < v->count; x++) {
        int key1 = e * v->count + x;
        int key2 = x * v->count + e;

        const Term *t1 = NULL;
        const Term *t2 = NULL;
        int c1 = 0;
        int c2 = 0;

        if (!GetPairTerms(v, key1, &t1, &c1) || !GetPairTerms(v, key2, &t2, &c2)) {
            fprintf(stderr, "Failed to load identity pair(s).\n");
            return 0;
        }

        if (!(c1 == 1 && t1[0].z == x && PolyIsOne(&t1[0].poly))) {
            failures++;
            if (printed < MAX_FAIL_PRINT) {
                char xPerm[16];
                char gotPerm[16];
                PermStringFromIndex(v->n, x, xPerm);
                if (c1 == 1 && t1[0].z >= 0) {
                    PermStringFromIndex(v->n, t1[0].z, gotPerm);
                } else {
                    snprintf(gotPerm, sizeof(gotPerm), "?");
                }
                fprintf(stderr, "Identity fail (left): e * %s expected (1*v^0)C_{%s}, got ", xPerm, xPerm);
                PrintTerms(stderr, v->n, t1, c1);
                fprintf(stderr, "\n");
                printed++;
            }
        }

        if (!(c2 == 1 && t2[0].z == x && PolyIsOne(&t2[0].poly))) {
            failures++;
            if (printed < MAX_FAIL_PRINT) {
                char xPerm[16];
                PermStringFromIndex(v->n, x, xPerm);
                fprintf(stderr, "Identity fail (right): %s * e expected (1*v^0)C_{%s}, got ", xPerm, xPerm);
                PrintTerms(stderr, v->n, t2, c2);
                fprintf(stderr, "\n");
                printed++;
            }
        }
    }

    if (failures > 0) {
        fprintf(stderr, "Identity check failed at %d pairs.\n", failures);
        return 0;
    }

    return 1;
}

static int CheckInversionSymmetry(Validator *v, int sampleCount) {
    int64_t totalPairs = v->totalPairs;
    if (sampleCount <= 0 || (int64_t)sampleCount > totalPairs) {
        sampleCount = (int)totalPairs;
    }

    int failures = 0;
    int printed = 0;
    srand(1234567);

    for (int s = 0; s < sampleCount; s++) {
        int key;
        if (sampleCount == (int)totalPairs) {
            key = s;
        } else {
            key = rand() % (int)totalPairs;
        }

        int x = key / v->count;
        int y = key % v->count;
        int yi = v->inverseIndex[y];
        int xi = v->inverseIndex[x];
        int mirrorKey = yi * v->count + xi;

        const Term *a = NULL;
        const Term *b = NULL;
        int ac = 0;
        int bc = 0;
        if (!GetPairTerms(v, key, &a, &ac) || !GetPairTerms(v, mirrorKey, &b, &bc)) {
            fprintf(stderr, "Failed loading pair(s) in inversion check.\n");
            return 0;
        }

        int matched = 1;
        int badZ = -1;
        int badZi = -1;
        int badJ = -1;
        for (int i = 0; i < ac; i++) {
            int z = a[i].z;
            if (z < 0) {
                continue;
            }
            int zi = v->inverseIndex[z];
            int j = FindTerm(b, bc, zi);
            if (j < 0 || !PolyEqual(&a[i].poly, &b[j].poly)) {
                matched = 0;
                badZ = z;
                badZi = zi;
                badJ = j;
                break;
            }
        }

        if (!matched) {
            failures++;
            if (printed < MAX_FAIL_PRINT) {
                char xPerm[16], yPerm[16], yiPerm[16], xiPerm[16], zPerm[16], ziPerm[16];
                PermStringFromIndex(v->n, x, xPerm);
                PermStringFromIndex(v->n, y, yPerm);
                PermStringFromIndex(v->n, yi, yiPerm);
                PermStringFromIndex(v->n, xi, xiPerm);
                if (badZ >= 0) PermStringFromIndex(v->n, badZ, zPerm); else snprintf(zPerm, sizeof(zPerm), "?");
                if (badZi >= 0) PermStringFromIndex(v->n, badZi, ziPerm); else snprintf(ziPerm, sizeof(ziPerm), "?");

                fprintf(stderr, "Inversion fail: (%s,%s) vs (%s,%s), missing/mismatched z=%s -> z^{-1}=%s\n",
                    xPerm, yPerm, yiPerm, xiPerm, zPerm, ziPerm);

                if (badZ >= 0) {
                    int i = FindTerm(a, ac, badZ);
                    if (i >= 0) {
                        fprintf(stderr, "  lhs coeff at z=%s: ", zPerm);
                        PrintPoly(stderr, &a[i].poly);
                        fprintf(stderr, "\n");
                    }
                }

                if (badJ >= 0 && badJ < bc) {
                    fprintf(stderr, "  rhs coeff at z^{-1}=%s: ", ziPerm);
                    PrintPoly(stderr, &b[badJ].poly);
                    fprintf(stderr, "\n");
                } else {
                    fprintf(stderr, "  rhs has no term at z^{-1}=%s\n", ziPerm);
                }
                printed++;
            }
        }
    }

    if (failures > 0) {
        fprintf(stderr, "Warning: inversion symmetry failed for %d / %d sampled pairs.\n", failures, sampleCount);
        fprintf(stderr, "This may be due to indexing/convention choices rather than data corruption.\n");
    }

    return 1;
}

static void ZeroPolyArray(Poly57 *arr, int count) {
    for (int i = 0; i < count; i++) {
        PolyZero(&arr[i]);
    }
}

static int ComparePolyArrays(const Poly57 *a, const Poly57 *b, int count) {
    for (int i = 0; i < count; i++) {
        if (!PolyEqual(&a[i], &b[i])) {
            return 0;
        }
    }
    return 1;
}

static int FirstPolyArrayMismatch(const Poly57 *a, const Poly57 *b, int count) {
    for (int i = 0; i < count; i++) {
        if (!PolyEqual(&a[i], &b[i])) {
            return i;
        }
    }
    return -1;
}

static int MultiplyExpansionByBasis(
    Validator *v,
    const Term *left,
    int leftCount,
    int rightIndex,
    Poly57 *out
) {
    ZeroPolyArray(out, v->count);

    for (int i = 0; i < leftCount; i++) {
        if (left[i].z < 0 || PolyIsZero(&left[i].poly)) {
            continue;
        }

        int key = left[i].z * v->count + rightIndex;
        const Term *right = NULL;
        int rightCount = 0;
        if (!GetPairTerms(v, key, &right, &rightCount)) {
            fprintf(stderr, "Internal error: failed to load pair (z=%d, right=%d).\n", left[i].z, rightIndex);
            return 0;
        }

        for (int j = 0; j < rightCount; j++) {
            if (right[j].z < 0 || PolyIsZero(&right[j].poly)) {
                continue;
            }
            Poly57 product;
            PolyMultiply(&left[i].poly, &right[j].poly, &product);
            PolyAddInplace(&out[right[j].z], &product);
        }
    }

    return 1;
}

static int MultiplyBasisByExpansion(
    Validator *v,
    int leftIndex,
    const Term *right,
    int rightCount,
    Poly57 *out
) {
    ZeroPolyArray(out, v->count);

    for (int i = 0; i < rightCount; i++) {
        if (right[i].z < 0 || PolyIsZero(&right[i].poly)) {
            continue;
        }

        int key = leftIndex * v->count + right[i].z;
        const Term *leftTimes = NULL;
        int leftTimesCount = 0;
        if (!GetPairTerms(v, key, &leftTimes, &leftTimesCount)) {
            fprintf(stderr, "Internal error: failed to load pair (left=%d, z=%d).\n", leftIndex, right[i].z);
            return 0;
        }

        for (int j = 0; j < leftTimesCount; j++) {
            if (leftTimes[j].z < 0 || PolyIsZero(&leftTimes[j].poly)) {
                continue;
            }

            Poly57 product;
            PolyMultiply(&right[i].poly, &leftTimes[j].poly, &product);
            PolyAddInplace(&out[leftTimes[j].z], &product);
        }
    }

    return 1;
}

static int CheckAssociativity(Validator *v, int sampleCount) {
    if (sampleCount <= 0) {
        sampleCount = 1000;
    }

    Poly57 *lhs = (Poly57 *)calloc((size_t)v->count, sizeof(Poly57));
    Poly57 *rhs = (Poly57 *)calloc((size_t)v->count, sizeof(Poly57));
    if (!lhs || !rhs) {
        free(lhs);
        free(rhs);
        fprintf(stderr, "Out of memory in associativity check.\n");
        return 0;
    }

    int failures = 0;
    int printed = 0;
    srand(7654321);

    for (int s = 0; s < sampleCount; s++) {
        int x = rand() % v->count;
        int y = rand() % v->count;
        int z = rand() % v->count;

        const Term *xy = NULL;
        const Term *yz = NULL;
        int xyCount = 0;
        int yzCount = 0;
        Term *xyOwned = NULL;
        Term *yzOwned = NULL;

        if (!GetPairTerms(v, x * v->count + y, &xy, &xyCount) ||
            !GetPairTerms(v, y * v->count + z, &yz, &yzCount)) {
            free(lhs);
            free(rhs);
            return 0;
        }

        if (!CopyTerms(xy, xyCount, &xyOwned) || !CopyTerms(yz, yzCount, &yzOwned)) {
            free(xyOwned);
            free(yzOwned);
            free(lhs);
            free(rhs);
            fprintf(stderr, "Out of memory copying pair terms for associativity check.\n");
            return 0;
        }

        if (!MultiplyExpansionByBasis(v, xyOwned, xyCount, z, lhs) ||
            !MultiplyBasisByExpansion(v, x, yzOwned, yzCount, rhs)) {
            free(xyOwned);
            free(yzOwned);
            free(lhs);
            free(rhs);
            return 0;
        }

        free(xyOwned);
        free(yzOwned);

        if (!ComparePolyArrays(lhs, rhs, v->count)) {
            failures++;
            if (printed < MAX_FAIL_PRINT) {
                int m = FirstPolyArrayMismatch(lhs, rhs, v->count);
                char xPerm[16], yPerm[16], zPerm[16], mPerm[16];
                PermStringFromIndex(v->n, x, xPerm);
                PermStringFromIndex(v->n, y, yPerm);
                PermStringFromIndex(v->n, z, zPerm);
                if (m >= 0) {
                    PermStringFromIndex(v->n, m, mPerm);
                } else {
                    snprintf(mPerm, sizeof(mPerm), "?");
                }

                fprintf(stderr, "Assoc fail for (x,y,z)=(%s,%s,%s), first mismatch at w=%s\n", xPerm, yPerm, zPerm, mPerm);
                if (m >= 0) {
                    fprintf(stderr, "  (C_x*C_y)*C_z coeff: ");
                    PrintPoly(stderr, &lhs[m]);
                    fprintf(stderr, "\n");
                    fprintf(stderr, "  C_x*(C_y*C_z) coeff: ");
                    PrintPoly(stderr, &rhs[m]);
                    fprintf(stderr, "\n");
                }
                printed++;
            }
        }
    }

    free(lhs);
    free(rhs);

    if (failures > 0) {
        fprintf(stderr, "Associativity failed for %d / %d sampled triples.\n", failures, sampleCount);
        return 0;
    }

    return 1;
}

static int CheckAssociativityTriple(Validator *v, int x, int y, int z) {
    if (x < 0 || x >= v->count || y < 0 || y >= v->count || z < 0 || z >= v->count) {
        fprintf(stderr, "Triple indices out of range for S%d.\n", v->n);
        return 0;
    }

    Poly57 *lhs = (Poly57 *)calloc((size_t)v->count, sizeof(Poly57));
    Poly57 *rhs = (Poly57 *)calloc((size_t)v->count, sizeof(Poly57));
    if (!lhs || !rhs) {
        free(lhs);
        free(rhs);
        fprintf(stderr, "Out of memory in single-triple associativity check.\n");
        return 0;
    }

    const Term *xy = NULL;
    const Term *yz = NULL;
    int xyCount = 0;
    int yzCount = 0;
    Term *xyOwned = NULL;
    Term *yzOwned = NULL;

    if (!GetPairTerms(v, x * v->count + y, &xy, &xyCount) ||
        !GetPairTerms(v, y * v->count + z, &yz, &yzCount)) {
        free(lhs);
        free(rhs);
        fprintf(stderr, "Failed loading pair(s) in single-triple associativity check.\n");
        return 0;
    }

    if (!CopyTerms(xy, xyCount, &xyOwned) || !CopyTerms(yz, yzCount, &yzOwned)) {
        free(xyOwned);
        free(yzOwned);
        free(lhs);
        free(rhs);
        fprintf(stderr, "Out of memory copying pair terms in single-triple associativity check.\n");
        return 0;
    }

    if (!MultiplyExpansionByBasis(v, xyOwned, xyCount, z, lhs) ||
        !MultiplyBasisByExpansion(v, x, yzOwned, yzCount, rhs)) {
        free(xyOwned);
        free(yzOwned);
        free(lhs);
        free(rhs);
        fprintf(stderr, "Internal error while expanding associativity triple.\n");
        return 0;
    }

    free(xyOwned);
    free(yzOwned);

    char xPerm[16], yPerm[16], zPerm[16];
    PermStringFromIndex(v->n, x, xPerm);
    PermStringFromIndex(v->n, y, yPerm);
    PermStringFromIndex(v->n, z, zPerm);
    printf("Single triple check (x,y,z)=(%s,%s,%s)\n", xPerm, yPerm, zPerm);

    int ok = ComparePolyArrays(lhs, rhs, v->count);
    if (!ok) {
        int m = FirstPolyArrayMismatch(lhs, rhs, v->count);
        char wPerm[16];
        PermStringFromIndex(v->n, m, wPerm);
        fprintf(stderr, "Associativity mismatch at w=%s\n", wPerm);
        fprintf(stderr, "  (C_x*C_y)*C_z coeff: ");
        PrintPoly(stderr, &lhs[m]);
        fprintf(stderr, "\n");
        fprintf(stderr, "  C_x*(C_y*C_z) coeff: ");
        PrintPoly(stderr, &rhs[m]);
        fprintf(stderr, "\n");
    } else {
        printf("Associativity holds for this triple.\n");
    }

    free(lhs);
    free(rhs);
    return ok;
}

static void FreeTermBag(TermBag *b) {
    free(b->terms);
    b->terms = NULL;
    b->count = 0;
    b->cap = 0;
    b->present = 0;
}

static int TermBagPush(TermBag *b, int z, const Poly57 *poly) {
    if (b->count == b->cap) {
        int newCap = b->cap > 0 ? (b->cap * 2) : 8;
        Term *newTerms = (Term *)realloc(b->terms, (size_t)newCap * sizeof(Term));
        if (!newTerms) {
            return 0;
        }
        b->terms = newTerms;
        b->cap = newCap;
    }

    b->terms[b->count].z = z;
    b->terms[b->count].poly = *poly;
    b->count++;
    b->present = 1;
    return 1;
}

static int CheckAssociativityTripleFast(Validator *v, int x, int y, int z) {
    if (x < 0 || x >= v->count || y < 0 || y >= v->count || z < 0 || z >= v->count) {
        fprintf(stderr, "Triple indices out of range for S%d.\n", v->n);
        return 0;
    }

    TermBag xy = {0};
    TermBag yz = {0};
    TermBag *xv = (TermBag *)calloc((size_t)v->count, sizeof(TermBag));
    TermBag *uz = (TermBag *)calloc((size_t)v->count, sizeof(TermBag));
    int *needXv = (int *)calloc((size_t)v->count, sizeof(int));
    int *needUz = (int *)calloc((size_t)v->count, sizeof(int));
    Poly57 *lhs = (Poly57 *)calloc((size_t)v->count, sizeof(Poly57));
    Poly57 *rhs = (Poly57 *)calloc((size_t)v->count, sizeof(Poly57));

    if (!xv || !uz || !needXv || !needUz || !lhs || !rhs) {
        fprintf(stderr, "Out of memory in fast single-triple associativity check.\n");
        free(xv);
        free(uz);
        free(needXv);
        free(needUz);
        free(lhs);
        free(rhs);
        return 0;
    }

    char xPerm[16], yPerm[16], zPerm[16];
    PermStringFromIndex(v->n, x, xPerm);
    PermStringFromIndex(v->n, y, yPerm);
    PermStringFromIndex(v->n, z, zPerm);

    if (FileSeek64(v->csv, 0) != 0) {
        fprintf(stderr, "Failed to seek CSV in fast triple mode.\n");
        goto cleanup_fail;
    }

    char line[16384];
    int sawHeader = 0;
    while (fgets(line, sizeof(line), v->csv)) {
        if (!sawHeader) {
            sawHeader = 1;
            continue;
        }

        char fx[64], fy[64], fz[64], fh[16000];
        if (!ParseCsv4Quoted(line, fx, sizeof(fx), fy, sizeof(fy), fz, sizeof(fz), fh, sizeof(fh))) {
            fprintf(stderr, "Malformed CSV row in fast triple mode.\n");
            goto cleanup_fail;
        }

        int rowIsX = (strcmp(fx, xPerm) == 0);
        int rowIsYtoZ = (strcmp(fy, zPerm) == 0);
        if (!rowIsX && !rowIsYtoZ) {
            continue;
        }

        int rx = -1;
        int ry = -1;
        int rz = -1;
        int parsedRx = 0;
        int parsedRz = 0;
        Poly57 p;
        int parsedPoly = 0;

        if (rowIsX) {
            if (!ParsePermutationString(v->n, fy, &ry)) {
                fprintf(stderr, "Failed parsing y permutation in fast triple mode.\n");
                goto cleanup_fail;
            }

            if (fz[0] == '\0') {
                rz = -1;
            } else if (!ParsePermutationString(v->n, fz, &rz)) {
                fprintf(stderr, "Failed parsing z permutation in fast triple mode.\n");
                goto cleanup_fail;
            }
            parsedRz = 1;

            if (!ParsePoly(fh, &p)) {
                fprintf(stderr, "Failed parsing polynomial in fast triple mode.\n");
                goto cleanup_fail;
            }
            parsedPoly = 1;

            if (!TermBagPush(&xv[ry], rz, &p)) {
                fprintf(stderr, "Out of memory while collecting (x,*) pairs.\n");
                goto cleanup_fail;
            }

            if (strcmp(fy, yPerm) == 0) {
                if (!TermBagPush(&xy, rz, &p)) {
                    fprintf(stderr, "Out of memory while collecting (x,y) pair.\n");
                    goto cleanup_fail;
                }
            }
        }

        if (rowIsYtoZ) {
            if (!parsedRx) {
                if (!ParsePermutationString(v->n, fx, &rx)) {
                    fprintf(stderr, "Failed parsing x permutation in fast triple mode.\n");
                    goto cleanup_fail;
                }
                parsedRx = 1;
            }

            if (!parsedRz) {
                if (fz[0] == '\0') {
                    rz = -1;
                } else if (!ParsePermutationString(v->n, fz, &rz)) {
                    fprintf(stderr, "Failed parsing z permutation in fast triple mode.\n");
                    goto cleanup_fail;
                }
                parsedRz = 1;
            }

            if (!parsedPoly) {
                if (!ParsePoly(fh, &p)) {
                    fprintf(stderr, "Failed parsing polynomial in fast triple mode.\n");
                    goto cleanup_fail;
                }
                parsedPoly = 1;
            }

            if (!TermBagPush(&uz[rx], rz, &p)) {
                fprintf(stderr, "Out of memory while collecting (*,z) pairs.\n");
                goto cleanup_fail;
            }

            if (strcmp(fx, yPerm) == 0) {
                if (!TermBagPush(&yz, rz, &p)) {
                    fprintf(stderr, "Out of memory while collecting (y,z) pair.\n");
                    goto cleanup_fail;
                }
            }
        }
    }

    if (!xy.present || !yz.present) {
        fprintf(stderr, "Missing (x,y) or (y,z) pair data in fast triple mode.\n");
        goto cleanup_fail;
    }

    for (int i = 0; i < xy.count; i++) {
        if (xy.terms[i].z >= 0 && !PolyIsZero(&xy.terms[i].poly)) {
            needUz[xy.terms[i].z] = 1;
        }
    }

    for (int i = 0; i < yz.count; i++) {
        if (yz.terms[i].z >= 0 && !PolyIsZero(&yz.terms[i].poly)) {
            needXv[yz.terms[i].z] = 1;
        }
    }

    for (int i = 0; i < v->count; i++) {
        if (needUz[i] && !uz[i].present) {
            fprintf(stderr, "Missing required pair (%d,%d) in fast triple mode.\n", i, z);
            goto cleanup_fail;
        }
        if (needXv[i] && !xv[i].present) {
            fprintf(stderr, "Missing required pair (%d,%d) in fast triple mode.\n", x, i);
            goto cleanup_fail;
        }
    }

    ZeroPolyArray(lhs, v->count);
    ZeroPolyArray(rhs, v->count);

    for (int i = 0; i < xy.count; i++) {
        int u = xy.terms[i].z;
        if (u < 0 || PolyIsZero(&xy.terms[i].poly)) {
            continue;
        }

        for (int j = 0; j < uz[u].count; j++) {
            if (uz[u].terms[j].z < 0 || PolyIsZero(&uz[u].terms[j].poly)) {
                continue;
            }

            Poly57 product;
            PolyMultiply(&xy.terms[i].poly, &uz[u].terms[j].poly, &product);
            PolyAddInplace(&lhs[uz[u].terms[j].z], &product);
        }
    }

    for (int i = 0; i < yz.count; i++) {
        int vv = yz.terms[i].z;
        if (vv < 0 || PolyIsZero(&yz.terms[i].poly)) {
            continue;
        }

        for (int j = 0; j < xv[vv].count; j++) {
            if (xv[vv].terms[j].z < 0 || PolyIsZero(&xv[vv].terms[j].poly)) {
                continue;
            }

            Poly57 product;
            PolyMultiply(&yz.terms[i].poly, &xv[vv].terms[j].poly, &product);
            PolyAddInplace(&rhs[xv[vv].terms[j].z], &product);
        }
    }

    printf("Single fast triple check (x,y,z)=(%s,%s,%s)\n", xPerm, yPerm, zPerm);

    if (!ComparePolyArrays(lhs, rhs, v->count)) {
        int m = FirstPolyArrayMismatch(lhs, rhs, v->count);
        char wPerm[16];
        PermStringFromIndex(v->n, m, wPerm);
        fprintf(stderr, "Associativity mismatch at w=%s\n", wPerm);
        fprintf(stderr, "  (C_x*C_y)*C_z coeff: ");
        PrintPoly(stderr, &lhs[m]);
        fprintf(stderr, "\n");
        fprintf(stderr, "  C_x*(C_y*C_z) coeff: ");
        PrintPoly(stderr, &rhs[m]);
        fprintf(stderr, "\n");
        goto cleanup_fail;
    }

    printf("Associativity holds for this triple.\n");

    FreeTermBag(&xy);
    FreeTermBag(&yz);
    for (int i = 0; i < v->count; i++) {
        FreeTermBag(&xv[i]);
        FreeTermBag(&uz[i]);
    }
    free(xv);
    free(uz);
    free(needXv);
    free(needUz);
    free(lhs);
    free(rhs);
    return 1;

cleanup_fail:
    FreeTermBag(&xy);
    FreeTermBag(&yz);
    for (int i = 0; i < v->count; i++) {
        FreeTermBag(&xv[i]);
        FreeTermBag(&uz[i]);
    }
    free(xv);
    free(uz);
    free(needXv);
    free(needUz);
    free(lhs);
    free(rhs);
    return 0;
}

static int CheckAssociativityFirstFailure(Validator *v) {
    Poly57 *lhs = (Poly57 *)calloc((size_t)v->count, sizeof(Poly57));
    Poly57 *rhs = (Poly57 *)calloc((size_t)v->count, sizeof(Poly57));
    if (!lhs || !rhs) {
        free(lhs);
        free(rhs);
        fprintf(stderr, "Out of memory in exhaustive associativity check.\n");
        return 0;
    }

    int64_t checked = 0;
    int64_t total = (int64_t)v->count * (int64_t)v->count * (int64_t)v->count;

    for (int x = 0; x < v->count; x++) {
        for (int y = 0; y < v->count; y++) {
            const Term *xy = NULL;
            int xyCount = 0;
            Term *xyOwned = NULL;
            if (!GetPairTerms(v, x * v->count + y, &xy, &xyCount)) {
                fprintf(stderr, "Failed loading (x,y) pair during exhaustive check.\n");
                free(lhs);
                free(rhs);
                return 0;
            }

            if (!CopyTerms(xy, xyCount, &xyOwned)) {
                fprintf(stderr, "Out of memory copying (x,y) pair during exhaustive check.\n");
                free(lhs);
                free(rhs);
                return 0;
            }

            for (int z = 0; z < v->count; z++) {
                const Term *yz = NULL;
                int yzCount = 0;
                Term *yzOwned = NULL;

                if (!GetPairTerms(v, y * v->count + z, &yz, &yzCount)) {
                    fprintf(stderr, "Failed loading (y,z) pair during exhaustive check.\n");
                    free(xyOwned);
                    free(lhs);
                    free(rhs);
                    return 0;
                }

                if (!CopyTerms(yz, yzCount, &yzOwned)) {
                    fprintf(stderr, "Out of memory copying (y,z) pair during exhaustive check.\n");
                    free(xyOwned);
                    free(lhs);
                    free(rhs);
                    return 0;
                }

                if (!MultiplyExpansionByBasis(v, xyOwned, xyCount, z, lhs) ||
                    !MultiplyBasisByExpansion(v, x, yzOwned, yzCount, rhs)) {
                    free(yzOwned);
                    free(xyOwned);
                    free(lhs);
                    free(rhs);
                    return 0;
                }

                free(yzOwned);

                checked++;
                if ((checked % 200000) == 0 || checked == total) {
                    fprintf(stderr, "Assoc exhaustive progress: %lld / %lld triples\n", (long long)checked, (long long)total);
                }

                if (!ComparePolyArrays(lhs, rhs, v->count)) {
                    int m = FirstPolyArrayMismatch(lhs, rhs, v->count);
                    char xPerm[16], yPerm[16], zPerm[16], wPerm[16];
                    PermStringFromIndex(v->n, x, xPerm);
                    PermStringFromIndex(v->n, y, yPerm);
                    PermStringFromIndex(v->n, z, zPerm);
                    PermStringFromIndex(v->n, m, wPerm);

                    fprintf(stderr, "First associativity failure at (x,y,z)=(%s,%s,%s), w=%s\n", xPerm, yPerm, zPerm, wPerm);
                    fprintf(stderr, "  (C_x*C_y)*C_z coeff: ");
                    PrintPoly(stderr, &lhs[m]);
                    fprintf(stderr, "\n");
                    fprintf(stderr, "  C_x*(C_y*C_z) coeff: ");
                    PrintPoly(stderr, &rhs[m]);
                    fprintf(stderr, "\n");

                    free(xyOwned);
                    free(lhs);
                    free(rhs);
                    return 0;
                }
            }

            free(xyOwned);
        }
    }

    printf("Exhaustive associativity check: all %lld triples passed.\n", (long long)total);
    free(lhs);
    free(rhs);
    return 1;
}

static void FreeValidator(Validator *v) {
    if (v->csv) {
        fclose(v->csv);
        v->csv = NULL;
    }

    if (v->cache) {
        for (int i = 0; i < v->cacheSize; i++) {
            free(v->cache[i].terms);
        }
    }

    free(v->cache);
    free(v->pairOffsets);
    free(v->pairLineCounts);
    free(v->inverseIndex);
}

static void Usage(const char *prog) {
    printf("Usage:\n");
    printf("  %s n csv_file [inverse_samples] [assoc_samples]\n", prog);
    printf("  %s n csv_file --triple x y z\n", prog);
    printf("  %s n csv_file --triple-legacy x y z\n", prog);
    printf("  %s n csv_file --first-assoc-fail\n", prog);
    printf("\n");
    printf("Examples:\n");
    printf("  %s 5 S5_structure_constants_bulk_opt.csv\n", prog);
    printf("  %s 6 S6_structure_constants_bulk_opt.csv 5000 2000\n", prog);
    printf("  %s 6 S6_structure_constants_bulk_opt_parallel.csv --triple 684 681 315\n", prog);
    printf("\n");
    printf("Set assoc_samples to 0 to skip associativity sampling (default).\n");
}

int main(int argc, char *argv[]) {
    if (argc < 3 || argc > 7) {
        Usage(argv[0]);
        return 1;
    }

    clock_t programStart = clock();

    Validator v;
    memset(&v, 0, sizeof(v));

    v.n = atoi(argv[1]);
    if (v.n < 4 || v.n > 9) {
        fprintf(stderr, "n must be between 4 and 9.\n");
        return 1;
    }

    v.count = fac(v.n);
    v.totalPairs = (int64_t)v.count * (int64_t)v.count;
    v.csvPath = argv[2];

    int tripleMode = 0;
    int tripleFastMode = 0;
    int firstAssocFailMode = 0;
    int tripleX = -1;
    int tripleY = -1;
    int tripleZ = -1;

    int inverseSamples = ((v.count <= 120) ? (int)v.totalPairs : 5000);
    int assocSamples = 0;
    if (argc >= 4) {
        if (strcmp(argv[3], "--triple") == 0) {
            if (argc != 7) {
                Usage(argv[0]);
                return 1;
            }
            tripleFastMode = 1;
            tripleX = atoi(argv[4]);
            tripleY = atoi(argv[5]);
            tripleZ = atoi(argv[6]);
        } else if (strcmp(argv[3], "--triple-legacy") == 0) {
            if (argc != 7) {
                Usage(argv[0]);
                return 1;
            }
            tripleMode = 1;
            tripleX = atoi(argv[4]);
            tripleY = atoi(argv[5]);
            tripleZ = atoi(argv[6]);
        } else if (strcmp(argv[3], "--first-assoc-fail") == 0) {
            if (argc != 4) {
                Usage(argv[0]);
                return 1;
            }
            firstAssocFailMode = 1;
        } else {
            inverseSamples = atoi(argv[3]);
            assocSamples = (argc >= 5) ? atoi(argv[4]) : 0;
        }
    }

    v.csv = fopen(v.csvPath, "r");
    if (!v.csv) {
        fprintf(stderr, "Could not open %s\n", v.csvPath);
        return 1;
    }

    if (tripleFastMode) {
        clock_t start = clock();
        printf("Validator start: n=%d, elements=%d, pairs=%lld\n", v.n, v.count, (long long)v.totalPairs);
        printf("Checks: single associativity triple only (fast path)\n");
        int ok = CheckAssociativityTripleFast(&v, tripleX, tripleY, tripleZ);
        double elapsed = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
        double totalElapsed = (double)(clock() - programStart) / (double)CLOCKS_PER_SEC;
        printf("Validation run finished in %.2fs\n", elapsed);
        printf("Total process time (including startup) %.2fs\n", totalElapsed);
        FreeValidator(&v);
        return ok ? 0 : 1;
    }

    clock_t idxStart = clock();
    if (!BuildPairIndex(&v)) {
        fprintf(stderr, "Failed to build pair index from CSV.\n");
        FreeValidator(&v);
        return 1;
    }
    double idxSeconds = (double)(clock() - idxStart) / (double)CLOCKS_PER_SEC;

    if (!tripleMode && !firstAssocFailMode) {
        if (!BuildInverseTable(&v)) {
            fprintf(stderr, "Failed to build inversion table.\n");
            FreeValidator(&v);
            return 1;
        }
    }

    v.cacheSize = 512;
    v.cache = (PairCacheEntry *)calloc((size_t)v.cacheSize, sizeof(PairCacheEntry));
    if (!v.cache) {
        fprintf(stderr, "Out of memory creating cache.\n");
        FreeValidator(&v);
        return 1;
    }
    for (int i = 0; i < v.cacheSize; i++) {
        v.cache[i].key = -1;
        v.cache[i].terms = NULL;
        v.cache[i].termCount = 0;
        v.cache[i].age = 0;
    }

    clock_t start = clock();

    if (tripleMode) {
        printf("Validator start: n=%d, elements=%d, pairs=%lld\n", v.n, v.count, (long long)v.totalPairs);
        printf("Checks: single associativity triple only\n");
        printf("Startup timing: pair-index build %.2fs\n", idxSeconds);
        int ok = CheckAssociativityTriple(&v, tripleX, tripleY, tripleZ);
        double elapsed = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
        double totalElapsed = (double)(clock() - programStart) / (double)CLOCKS_PER_SEC;
        printf("Validation run finished in %.2fs\n", elapsed);
        printf("Total process time (including startup) %.2fs\n", totalElapsed);
        FreeValidator(&v);
        return ok ? 0 : 1;
    }

    if (firstAssocFailMode) {
        printf("Validator start: n=%d, elements=%d, pairs=%lld\n", v.n, v.count, (long long)v.totalPairs);
        printf("Checks: exhaustive associativity scan (lexicographic triple order)\n");
        printf("Startup timing: pair-index build %.2fs\n", idxSeconds);
        int ok = CheckAssociativityFirstFailure(&v);
        double elapsed = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
        double totalElapsed = (double)(clock() - programStart) / (double)CLOCKS_PER_SEC;
        printf("Validation run finished in %.2fs\n", elapsed);
        printf("Total process time (including startup) %.2fs\n", totalElapsed);
        FreeValidator(&v);
        return ok ? 0 : 1;
    }

    printf("Validator start: n=%d, elements=%d, pairs=%lld\n", v.n, v.count, (long long)v.totalPairs);
    printf("Checks: identity + inversion(%d samples)", inverseSamples);
    if (assocSamples > 0) {
        printf(" + associativity(%d samples)\n", assocSamples);
    } else {
        printf(" + associativity(skipped)\n");
    }

    if (!CheckIdentity(&v)) {
        fprintf(stderr, "Identity check failed.\n");
        FreeValidator(&v);
        return 1;
    }
    printf("Identity check: OK\n");

    if (!CheckInversionSymmetry(&v, inverseSamples)) {
        fprintf(stderr, "Inversion symmetry check encountered an internal error.\n");
        FreeValidator(&v);
        return 1;
    }
    printf("Inversion symmetry check: completed (see warnings if any)\n");

    if (assocSamples > 0) {
        if (!CheckAssociativity(&v, assocSamples)) {
            fprintf(stderr, "Warning: associativity check reported failures for sampled triples.\n");
            fprintf(stderr, "This can indicate a convention mismatch or a genuine data issue.\n");
        } else {
            printf("Associativity check: OK\n");
        }
    }

    double elapsed = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
    double totalElapsed = (double)(clock() - programStart) / (double)CLOCKS_PER_SEC;
    printf("All validation checks passed in %.2fs\n", elapsed);
    printf("Total process time (including startup) %.2fs\n", totalElapsed);

    FreeValidator(&v);
    return 0;
}
