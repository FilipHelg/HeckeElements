#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lehmer.h"
#include "laurent.h"

static Laurent_t NegateLaurent(Laurent_t p) {
    Laurent_t result = ZeroInitializeLaurent();
    for (int i = 0; i < 57; i++) {
        result.coeff[i] = -p.coeff[i];
    }
    return result;
}

static int ParsePermutationString(int n, const char *token, int *outIndex) {
    if ((int)strlen(token) != n) {
        return 0;
    }

    char perm[8] = {0};
    int seen[8] = {0};
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

static int ParseElementArg(int n, const char *arg, int *outIndex) {
    int maybePerm = 1;
    if ((int)strlen(arg) != n) {
        maybePerm = 0;
    } else {
        for (int i = 0; i < n; i++) {
            if (!isdigit((unsigned char)arg[i])) {
                maybePerm = 0;
                break;
            }
        }
    }

    if (maybePerm && ParsePermutationString(n, arg, outIndex)) {
        return 1;
    }

    char *endptr = NULL;
    long value = strtol(arg, &endptr, 10);
    if (endptr == arg || *endptr != '\0') {
        return 0;
    }
    *outIndex = (int)value;
    return 1;
}

static void PrintPermutation(int n, int index) {
    char perm[8] = {0};
    IndexToPerm(n, index, perm);
    for (int i = 0; i < n; i++) {
        printf("%d", perm[i]);
    }
}

static int GetStructureConstantPos(int exponent) {
    int pos = 28 + exponent;
    if (pos < 0 || pos >= 57) {
        return -1;
    }
    return pos;
}

static void AddLaurentAtExponent(Laurent_t *poly, int exponent, int coeff) {
    int pos = GetStructureConstantPos(exponent);
    if (pos >= 0) {
        poly->coeff[pos] += coeff;
    }
}

static Laurent_t *MultiplySimpleHeckeLeft(int n, int s, const Laurent_t H[], const int lengths[]) {
    int count = fac(n);
    Laurent_t *product = (Laurent_t *)calloc((size_t)count, sizeof(Laurent_t));
    Laurent_t lowerTermFactor = ZeroInitializeLaurent();
    lowerTermFactor.coeff[27] = 1;
    lowerTermFactor.coeff[29] = -1;

    for (int i = 0; i < count; i++) {
        if (!HasNonZero(&H[i])) {
            continue;
        }

        int composition = MultiplyIndex(n, s, i);
        product[composition] = SumLaurent(product[composition], H[i]);

        if (lengths[composition] != lengths[i] + 1) {
            Laurent_t lowerTerm = MultiplyLaurent(lowerTermFactor, H[i]);
            product[i] = SumLaurent(product[i], lowerTerm);
        }
    }

    return product;
}

static Laurent_t *CopyHeckeElement(int count, const Laurent_t source[]) {
    Laurent_t *result = (Laurent_t *)calloc((size_t)count, sizeof(Laurent_t));
    for (int i = 0; i < count; i++) {
        result[i] = source[i];
    }
    return result;
}

static Laurent_t *LeftMultiplyByIndex(int n, int index, const Laurent_t H[], const int lengths[]) {
    int count = fac(n);
    Laurent_t *term = CopyHeckeElement(count, H);
    char expression[28] = {0};
    int len = ReducedExpression(n, index, expression);

    for (int k = 0; k < len; k++) {
        int generator = (int)expression[k];
        int s = fac(n - generator);
        Laurent_t *next = MultiplySimpleHeckeLeft(n, s, term, lengths);
        free(term);
        term = next;
    }

    return term;
}

static Laurent_t *MultiplyHeckeElements(int n, const Laurent_t H1[], const Laurent_t H2[], const int lengths[]) {
    int count = fac(n);
    Laurent_t *product = (Laurent_t *)calloc((size_t)count, sizeof(Laurent_t));

    for (int i = 0; i < count; i++) {
        if (!HasNonZero(&H1[i])) {
            continue;
        }

        Laurent_t *leftMultiplied = LeftMultiplyByIndex(n, i, H2, lengths);
        Laurent_t scalar = H1[i];
        for (int j = 0; j < count; j++) {
            Laurent_t scaledTerm = MultiplyLaurent(scalar, leftMultiplied[j]);
            product[j] = SumLaurent(product[j], scaledTerm);
        }
        free(leftMultiplied);
    }

    return product;
}

static int LoadKLBasisFromFile(int n, const char *filename, Laurent_t *klBasis, const int lengths[]) {
    int count = fac(n);

    for (int w = 0; w < count; w++) {
        klBasis[w * count + w].coeff[28] = 1;
        for (int y = 0; y < count; y++) {
            if (y == w) {
                continue;
            }
            if (!BruhatSmaller2(n, y, w)) {
                continue;
            }

            int diff = lengths[w] - lengths[y];
            AddLaurentAtExponent(&klBasis[w * count + y], diff, 1);
        }
    }

    FILE *fptr = fopen(filename, "r");
    if (!fptr) {
        fprintf(stderr, "Could not open %s\n", filename);
        return 1;
    }

    char line[256];
    while (fgets(line, sizeof(line), fptr) != NULL) {
        char yPerm[16] = {0};
        char wPerm[16] = {0};
        char coeffToken[128] = {0};

        if (sscanf(line, "%15s %15s %127s", yPerm, wPerm, coeffToken) != 3) {
            continue;
        }

        int y = -1;
        int w = -1;
        if (!ParsePermutationString(n, yPerm, &y) || !ParsePermutationString(n, wPerm, &w)) {
            fprintf(stderr, "Skipping malformed line: %s", line);
            continue;
        }

        if (!BruhatSmaller2(n, y, w)) {
            fprintf(stderr, "Skipping non-Bruhat pair in data: ");
            PrintPermutation(n, y);
            printf(" ");
            PrintPermutation(n, w);
            printf("\n");
            continue;
        }

        int diff = lengths[w] - lengths[y];
        AddLaurentAtExponent(&klBasis[w * count + y], diff, -1);

        char coeffBuffer[128] = {0};
        snprintf(coeffBuffer, sizeof(coeffBuffer), "%s", coeffToken);
        char *part = strtok(coeffBuffer, ",");
        int k = 0;
        while (part != NULL) {
            int coeff = atoi(part);
            int exponent = diff - 2 * k;
            AddLaurentAtExponent(&klBasis[w * count + y], exponent, coeff);
            part = strtok(NULL, ",");
            k++;
        }
    }

    fclose(fptr);
    return 0;
}

static Laurent_t *DecomposeHeckeToKLBasis(int n, const Laurent_t heckeElement[], const Laurent_t *klBasis, const int lengths[]) {
    int count = fac(n);
    int maxLength = n * (n - 1) / 2;
    Laurent_t *residual = CopyHeckeElement(count, heckeElement);
    Laurent_t *coeffs = (Laurent_t *)calloc((size_t)count, sizeof(Laurent_t));

    for (int currentLength = maxLength; currentLength >= 0; currentLength--) {
        for (int z = 0; z < count; z++) {
            if (lengths[z] != currentLength) {
                continue;
            }

            Laurent_t top = residual[z];
            if (!HasNonZero(&top)) {
                continue;
            }

            coeffs[z] = top;
            Laurent_t minusTop = NegateLaurent(top);
            const Laurent_t *basisElement = klBasis + z * count;

            for (int y = 0; y < count; y++) {
                Laurent_t correction = MultiplyLaurent(minusTop, basisElement[y]);
                residual[y] = SumLaurent(residual[y], correction);
            }
        }
    }

    for (int i = 0; i < count; i++) {
        if (HasNonZero(&residual[i])) {
            fprintf(stderr, "Warning: non-zero residual remained after KL decomposition.\n");
            break;
        }
    }

    free(residual);
    return coeffs;
}

static void PrintKLEquation(int n, int x, int y, const Laurent_t coeffs[]) {
    int count = fac(n);
    printf("C_{");
    PrintPermutation(n, x);
    printf("} * C_{");
    PrintPermutation(n, y);
    printf("} = ");

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
        printf(")C_{");
        PrintPermutation(n, z);
        printf("}");
        first = 0;
    }

    if (first) {
        printf("0");
    }
    printf("\n");
}

static void PrintUsage(const char *programName) {
    printf("Usage:\n");
    printf("  %s n\n", programName);
    printf("  %s n x y\n", programName);
    printf("\n");
    printf("n must be 4 or 5.\n");
    printf("x and y can be Lehmer indices or permutation strings (e.g. 0213 for S4).\n");
}

int main(int argc, char *argv[]) {
    if (argc != 2 && argc != 4) {
        PrintUsage(argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    if (n != 4 && n != 5) {
        fprintf(stderr, "This proof-of-concept currently supports only S4 and S5.\n");
        return 1;
    }

    const char *dataFile = (n == 4) ? "S4.txt" : "S5.txt";
    int count = fac(n);
    int *lengths = (int *)calloc((size_t)count, sizeof(int));
    Laurent_t *klBasis = (Laurent_t *)calloc((size_t)count * (size_t)count, sizeof(Laurent_t));
    if (!lengths || !klBasis) {
        fprintf(stderr, "Out of memory while allocating KL basis data.\n");
        free(lengths);
        free(klBasis);
        return 1;
    }

    for (int i = 0; i < count; i++) {
        lengths[i] = IndexToLength(n, i);
    }

    if (LoadKLBasisFromFile(n, dataFile, klBasis, lengths) != 0) {
        free(lengths);
        free(klBasis);
        return 1;
    }

    int x = fac(n - 1);
    int y = fac(n - 1);
    if (argc == 4) {
        if (!ParseElementArg(n, argv[2], &x) || !ParseElementArg(n, argv[3], &y)) {
            fprintf(stderr, "Could not parse x/y. Use Lehmer index or permutation string.\n");
            free(lengths);
            free(klBasis);
            return 1;
        }
    }

    if (x < 0 || x >= count || y < 0 || y >= count) {
        fprintf(stderr, "x and y must be in [0, %d].\n", count - 1);
        free(lengths);
        free(klBasis);
        return 1;
    }

    const Laurent_t *Cx = klBasis + x * count;
    const Laurent_t *Cy = klBasis + y * count;
    Laurent_t *productInH = MultiplyHeckeElements(n, Cx, Cy, lengths);
    Laurent_t *structureConstants = DecomposeHeckeToKLBasis(n, productInH, klBasis, lengths);

    printf("Loaded KL data from %s for S_%d (%d elements).\n", dataFile, n, count);
    PrintKLEquation(n, x, y, structureConstants);

    free(productInH);
    free(structureConstants);
    free(lengths);
    free(klBasis);
    return 0;
}