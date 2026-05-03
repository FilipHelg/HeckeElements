#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lehmer.h"
#include "laurent.h"

typedef struct {
    int64_t key;
    Laurent_t poly;
} NonTrivialEntry;

static int64_t PairKey(int count, int y, int w) {
    return (int64_t)y * (int64_t)count + (int64_t)w;
}

static int FacLocal(int n) {
    int value = 1;
    for (int i = 2; i <= n; i++) {
        value *= i;
    }
    return value;
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

static void AddAtExponent(Laurent_t *poly, int exponent, int coeff) {
    int pos = 28 + exponent;
    if (pos < 0 || pos >= 57 || coeff == 0) {
        return;
    }
    poly->coeff[pos] += coeff;
}

static Laurent_t MultiplyLaurentSafe(Laurent_t l1, Laurent_t l2) {
    Laurent_t product = ZeroInitializeLaurent();

    for (int i = 0; i < 57; i++) {
        if (l1.coeff[i] == 0) {
            continue;
        }
        for (int j = 0; j < 57; j++) {
            if (l2.coeff[j] == 0) {
                continue;
            }

            int index = i + j - 28;
            if (index >= 0 && index < 57) {
                product.coeff[index] += l1.coeff[i] * l2.coeff[j];
            }
        }
    }

    return product;
}

static const Laurent_t *FindNonTrivial(
    const NonTrivialEntry *entries,
    size_t entryCount,
    int count,
    int y,
    int w
) {
    if (entryCount == 0) {
        return NULL;
    }

    NonTrivialEntry needle;
    needle.key = PairKey(count, y, w);

    NonTrivialEntry *found = (NonTrivialEntry *)bsearch(
        &needle,
        entries,
        entryCount,
        sizeof(NonTrivialEntry),
        CompareEntries
    );

    if (!found) {
        return NULL;
    }
    return &found->poly;
}

static Laurent_t KLCoefficient(
    int n,
    int count,
    const int *lengths,
    const NonTrivialEntry *entries,
    size_t entryCount,
    int y,
    int w
) {
    Laurent_t p = ZeroInitializeLaurent();

    if (y == w) {
        p.coeff[28] = 1;
        return p;
    }

    if (!CheckBruhatSmallerSafe(n, y, w)) {
        return p;
    }

    const Laurent_t *nt = FindNonTrivial(entries, entryCount, count, y, w);
    if (nt) {
        return *nt;
    }

        AddAtExponent(&p, lengths[w] - lengths[y], 1);
    return p;
}

static int LoadNonTrivialData(
    int n,
    int count,
    const int *lengths,
    const char *filename,
    NonTrivialEntry **outEntries,
    size_t *outEntryCount
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
        if (!ParsePermutationString(n, yPerm, &y) || !ParsePermutationString(n, wPerm, &w)) {
            continue;
        }

        Laurent_t p = ZeroInitializeLaurent();
        int diff = lengths[w] - lengths[y];

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

        entries[used].key = PairKey(count, y, w);
        entries[used].poly = p;
        used++;
    }

    fclose(f);

    qsort(entries, used, sizeof(NonTrivialEntry), CompareEntries);
    *outEntries = entries;
    *outEntryCount = used;
    return 1;
}

static int BuildLengthAscendingOrder(int count, const int *lengths, int *order) {
    for (int i = 0; i < count; i++) {
        order[i] = i;
    }

    /* Portable fallback: insertion sort for small proof-of-concept sizes. */
    for (int i = 1; i < count; i++) {
        int key = order[i];
        int j = i - 1;
        while (j >= 0) {
            int oj = order[j];
            int less = 0;
            if (lengths[oj] > lengths[key]) {
                less = 1;
            } else if (lengths[oj] == lengths[key] && oj > key) {
                less = 1;
            }
            if (!less) {
                break;
            }
            order[j + 1] = order[j];
            j--;
        }
        order[j + 1] = key;
    }

    return 1;
}

static void PrintPermutation(int n, int index) {
    int perm[8] = {0};
    IndexToPermSafe(n, index, perm);
    for (int i = 0; i < n; i++) {
        printf("%d", perm[i]);
    }
}

static void PrintDualElement(int n, int x, int count, const Laurent_t *coeffs) {
    printf("Dual KL element for x = %d (", x);
    PrintPermutation(n, x);
    printf("):\n");

    printf("D_{x} = ");
    int first = 1;
    for (int z = 0; z < count; z++) {
        if (!HasNonZero(&coeffs[z])) {
            continue;
        }

        if (!first) {
            printf(" + ");
        }

        printf("(");
        DisplayLaurentPoly(coeffs[z]);
        printf(")H[");
        PrintPermutation(n, z);
        printf("]");
        first = 0;
    }

    if (first) {
        printf("0");
    }
    printf("\n");
}

static void PrintUsage(const char *prog) {
    printf("Usage:\n");
    printf("  %s n x\n", prog);
    printf("\n");
    printf("n determines S_n and must be between 3 and 8.\n");
    printf("x must be a Lehmer index in [0, n! - 1].\n");
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        PrintUsage(argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    if (n < 3 || n > 8) {
        fprintf(stderr, "n must be between 3 and 8 for this proof-of-concept.\n");
        return 1;
    }

    int count = FacLocal(n);

    char *endptr = NULL;
    long xLong = strtol(argv[2], &endptr, 10);
    if (endptr == argv[2] || *endptr != '\0') {
        fprintf(stderr, "x must be a Lehmer index (integer).\n");
        return 1;
    }

    if (xLong < 0 || xLong >= count) {
        fprintf(stderr, "x must satisfy 0 <= x < %d.\n", count);
        return 1;
    }
    int x = (int)xLong;

    int *lengths = (int *)calloc((size_t)count, sizeof(int));
    int *order = (int *)calloc((size_t)count, sizeof(int));
    Laurent_t *dualCoeff = (Laurent_t *)calloc((size_t)count, sizeof(Laurent_t));
    if (!lengths || !order || !dualCoeff) {
        fprintf(stderr, "Out of memory while allocating workspace.\n");
        free(lengths);
        free(order);
        free(dualCoeff);
        return 1;
    }

    for (int i = 0; i < count; i++) {
        lengths[i] = IndexToLengthSafe(n, i);
    }

    if (!BuildLengthAscendingOrder(count, lengths, order)) {
        fprintf(stderr, "Failed to build Bruhat-length order.\n");
        free(lengths);
        free(order);
        free(dualCoeff);
        return 1;
    }

    NonTrivialEntry *entries = NULL;
    size_t entryCount = 0;
    if (n >= 4) {
        char dataFile[32];
        snprintf(dataFile, sizeof(dataFile), "S%d.txt", n);
        if (!LoadNonTrivialData(n, count, lengths, dataFile, &entries, &entryCount)) {
            free(lengths);
            free(order);
            free(dualCoeff);
            return 1;
        }
    }

    /* Solve P * d = e_x by forward substitution in increasing Bruhat length order. */
    for (int oi = 0; oi < count; oi++) {
        int y = order[oi];

        Laurent_t rhs = ZeroInitializeLaurent();
        if (y == x) {
            rhs.coeff[28] = 1;
        }

        Laurent_t sum = ZeroInitializeLaurent();
        for (int oj = 0; oj < oi; oj++) {
            int z = order[oj];
            Laurent_t pyz = KLCoefficient(n, count, lengths, entries, entryCount, z, y);
            if (!HasNonZero(&pyz) || !HasNonZero(&dualCoeff[z])) {
                continue;
            }

            Laurent_t term = MultiplyLaurentSafe(pyz, dualCoeff[z]);
            sum = SumLaurent(sum, term);
        }

        Laurent_t minusSum = ZeroInitializeLaurent();
        for (int k = 0; k < 57; k++) {
            minusSum.coeff[k] = -sum.coeff[k];
        }
        dualCoeff[y] = SumLaurent(rhs, minusSum);
    }

    PrintDualElement(n, x, count, dualCoeff);

    free(entries);
    free(lengths);
    free(order);
    free(dualCoeff);
    return 0;
}
