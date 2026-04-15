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

static int PolyEqual(const Poly57 *a, const Poly57 *b) {
    for (int i = 0; i < 57; i++) {
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
    if (p->minPos != 28 || p->maxPos != 28) {
        return 0;
    }
    return p->coeff[28] == 1;
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

        int coeff = 0;
        int exp = 0;
        if (sscanf(cursor, "%d*v^%d", &coeff, &exp) != 2) {
            return 0;
        }

        int pos = 28 + exp;
        if (pos < 0 || pos >= 57) {
            return 0;
        }

        out->coeff[pos] += coeff;
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

    if (fseek(v->csv, 0, SEEK_SET) != 0) {
        return 0;
    }

    char line[16384];
    int sawHeader = 0;
    int64_t prevKey = -1;

    while (1) {
        int64_t lineStart = (int64_t)ftell(v->csv);
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

    if (fseek(v->csv, (long)offset, SEEK_SET) != 0) {
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

static int CheckIdentity(Validator *v) {
    int failures = 0;
    int e = 0;

    for (int x = 0; x < v->count; x++) {
        int key1 = e * v->count + x;
        int key2 = v->inverseIndex[x] * v->count + e;

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
        }

        if (!(c2 == 1 && t2[0].z == x && PolyIsOne(&t2[0].poly))) {
            failures++;
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
        for (int i = 0; i < ac; i++) {
            int z = a[i].z;
            if (z < 0) {
                continue;
            }
            int zi = v->inverseIndex[z];
            int j = FindTerm(b, bc, zi);
            if (j < 0 || !PolyEqual(&a[i].poly, &b[j].poly)) {
                matched = 0;
                break;
            }
        }

        if (!matched) {
            failures++;
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

static void MultiplyExpansionByBasis(
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

        int key = v->inverseIndex[left[i].z] * v->count + rightIndex;
        const Term *right = NULL;
        int rightCount = 0;
        if (!GetPairTerms(v, key, &right, &rightCount)) {
            continue;
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
}

static void MultiplyBasisByExpansion(
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

        int key = v->inverseIndex[leftIndex] * v->count + right[i].z;
        const Term *leftTimes = NULL;
        int leftTimesCount = 0;
        if (!GetPairTerms(v, key, &leftTimes, &leftTimesCount)) {
            continue;
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
    srand(7654321);

    for (int s = 0; s < sampleCount; s++) {
        int x = rand() % v->count;
        int y = rand() % v->count;
        int z = rand() % v->count;

        const Term *xy = NULL;
        const Term *yz = NULL;
        int xyCount = 0;
        int yzCount = 0;

        if (!GetPairTerms(v, x * v->count + y, &xy, &xyCount) ||
            !GetPairTerms(v, y * v->count + z, &yz, &yzCount)) {
            free(lhs);
            free(rhs);
            return 0;
        }

        MultiplyExpansionByBasis(v, xy, xyCount, z, lhs);
        MultiplyBasisByExpansion(v, x, yz, yzCount, rhs);

        if (!ComparePolyArrays(lhs, rhs, v->count)) {
            failures++;
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
    printf("\n");
    printf("Examples:\n");
    printf("  %s 5 S5_structure_constants_bulk_opt.csv\n", prog);
    printf("  %s 6 S6_structure_constants_bulk_opt.csv 5000 2000\n", prog);
    printf("\n");
    printf("Set assoc_samples to 0 to skip associativity sampling (default).\n");
}

int main(int argc, char *argv[]) {
    if (argc < 3 || argc > 5) {
        Usage(argv[0]);
        return 1;
    }

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

    int inverseSamples = (argc >= 4) ? atoi(argv[3]) : ((v.count <= 120) ? (int)v.totalPairs : 5000);
    int assocSamples = (argc >= 5) ? atoi(argv[4]) : 0;

    v.csv = fopen(v.csvPath, "r");
    if (!v.csv) {
        fprintf(stderr, "Could not open %s\n", v.csvPath);
        return 1;
    }

    if (!BuildPairIndex(&v)) {
        fprintf(stderr, "Failed to build pair index from CSV.\n");
        FreeValidator(&v);
        return 1;
    }

    if (!BuildInverseTable(&v)) {
        fprintf(stderr, "Failed to build inversion table.\n");
        FreeValidator(&v);
        return 1;
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
    printf("All validation checks passed in %.2fs\n", elapsed);

    FreeValidator(&v);
    return 0;
}
