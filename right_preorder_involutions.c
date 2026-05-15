#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

#include <string.h>

typedef struct {
    int n;
    int d;
    const char *klDataPath;
    int useAllMultipliers;
} RightPreorderArgs;

static int IndexToLengthSafeLocal(int n, int index) {
    return IndexToLength(n, index);
}

static int IsInvolution(int n, int w) {
    return MultiplyIndex(n, w, w) == 0;
}

static void PrintRightPreorderUsage(const char *prog) {
    printf("Usage:\n");
    printf("  %s n d [--kl-data PATH]\n", prog);
    printf("\n");
    printf("Print all involutions d' in S_n satisfying d' <=_R d.\n");
    printf("n must be in {3,4,5,6,7,8,9}. d is a Lehmer index in [0, n!-1].\n");
    printf("Default KL data path: S<n>.txt\n");
}

static int ParseArgs(int argc, char *argv[], RightPreorderArgs *args) {
    if (argc < 3) {
        return 0;
    }

    args->n = atoi(argv[1]);
    args->d = atoi(argv[2]);
    args->klDataPath = NULL;

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--kl-data") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--kl-data requires a path.\n");
                return 0;
            }
            args->klDataPath = argv[++i];
            continue;
        }

        if (strcmp(argv[i], "--all-multipliers") == 0) {
            args->useAllMultipliers = 1;
            continue;
        }

        fprintf(stderr, "Unknown option: %s\n", argv[i]);
        return 0;
    }

    return 1;
}

int main(int argc, char *argv[]) {
    RightPreorderArgs args;
    if (!ParseArgs(argc, argv, &args)) {
        PrintRightPreorderUsage(argv[0]);
        return 1;
    }

    if (args.n < 3 || args.n > 9) {
        fprintf(stderr, "n must be between 3 and 9.\n");
        return 1;
    }

    int count = fac(args.n);
    if (args.d < 0 || args.d >= count) {
        fprintf(stderr, "d must be a valid Lehmer index in [0, %d].\n", count - 1);
        return 1;
    }

    if (!IsInvolution(args.n, args.d)) {
        char perm[8] = {0};
        IndexToPerm(args.n, args.d, perm);
        fprintf(stderr, "Error: d=%d is not an involution in S_%d.\n", args.d, args.n);
        fprintf(stderr, "Permutation: ");
        for (int i = 0; i < args.n; i++) {
            fprintf(stderr, "%d ", (int)perm[i]);
        }
        fprintf(stderr, "\n");
        return 1;
    }

    char defaultKlData[32];
    snprintf(defaultKlData, sizeof(defaultKlData), "S%d.txt", args.n);
    if (!args.klDataPath) {
        args.klDataPath = defaultKlData;
    }

    FastContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = args.n;
    ctx.count = count;
    ctx.maxLength = args.n * (args.n - 1) / 2;

    ctx.lengths = (int *)calloc((size_t)ctx.count, sizeof(int));
    if (!ctx.lengths) {
        fprintf(stderr, "Out of memory for length table.\n");
        return 1;
    }

    for (int i = 0; i < ctx.count; i++) {
        ctx.lengths[i] = IndexToLengthSafeLocal(args.n, i);
    }

    if (!BuildLengthOrder(&ctx) || !BuildGeneratorTables(&ctx) || !BuildExpressionTables(&ctx)) {
        fprintf(stderr, "Failed to build context tables.\n");
        FreeContext(&ctx);
        return 1;
    }

    if (args.n >= 4) {
        if (!LoadNonTrivialData(&ctx, args.klDataPath)) {
            FreeContext(&ctx);
            return 1;
        }
    }

    if (!MaybeBuildKLBasisCache(&ctx)) {
        fprintf(stderr, "KL basis cache build failed; continuing without cache.\n");
    }

    SparseHecke Ccur;
    SparseHecke product;
    SparseHecke tmp1;
    SparseHecke tmp2;
    SparseHecke residual;
    SparseHecke basisTmp;
    DenseAccum accum;

    if (!SparseInit(&Ccur, ctx.count) ||
        !SparseInit(&product, ctx.count) ||
        !SparseInit(&tmp1, ctx.count) ||
        !SparseInit(&tmp2, ctx.count) ||
        !SparseInit(&residual, ctx.count) ||
        !SparseInit(&basisTmp, ctx.count) ||
        !DenseAccumInit(&accum, ctx.count)) {
        fprintf(stderr, "Failed to allocate sparse scratch storage.\n");
        SparseFree(&Ccur);
        SparseFree(&product);
        SparseFree(&tmp1);
        SparseFree(&tmp2);
        SparseFree(&residual);
        SparseFree(&basisTmp);
        DenseAccumFree(&accum);
        FreeContext(&ctx);
        return 1;
    }

    int multiplierCount = args.useAllMultipliers ? ctx.count : ctx.genCount;
    SparseHecke *Cs = (SparseHecke *)calloc((size_t)multiplierCount, sizeof(SparseHecke));
    if (!Cs) {
        fprintf(stderr, "Out of memory for simple generators.\n");
        SparseFree(&Ccur);
        SparseFree(&product);
        SparseFree(&tmp1);
        SparseFree(&tmp2);
        SparseFree(&residual);
        SparseFree(&basisTmp);
        DenseAccumFree(&accum);
        FreeContext(&ctx);
        return 1;
    }

    for (int m = 0; m < multiplierCount; m++) {
        if (!SparseInit(&Cs[m], ctx.count)) {
            fprintf(stderr, "Failed to allocate multiplier storage.\n");
            for (int j = 0; j < m; j++) {
                SparseFree(&Cs[j]);
            }
            free(Cs);
            SparseFree(&Ccur);
            SparseFree(&product);
            SparseFree(&tmp1);
            SparseFree(&tmp2);
            SparseFree(&residual);
            SparseFree(&basisTmp);
            DenseAccumFree(&accum);
            FreeContext(&ctx);
            return 1;
        }
        int idx = args.useAllMultipliers ? m : ctx.genLehmer[m];
        BuildKLElementSparse(&ctx, idx, &Cs[m]);
    }

    unsigned char *visited = (unsigned char *)calloc((size_t)ctx.count, sizeof(unsigned char));
    int *queue = (int *)calloc((size_t)ctx.count, sizeof(int));
    int *pred = (int *)calloc((size_t)ctx.count, sizeof(int));
    int *predGen = (int *)calloc((size_t)ctx.count, sizeof(int));
    int *constTerm = (int *)calloc((size_t)ctx.count, sizeof(int));
    Poly57 *coeffs = (Poly57 *)calloc((size_t)ctx.count, sizeof(Poly57));
    if (!visited || !queue || !coeffs) {
        fprintf(stderr, "Out of memory for traversal.\n");
        free(coeffs);
        free(queue);
        free(visited);
        for (int g = 0; g < ctx.genCount; g++) {
            SparseFree(&Cs[g]);
        }
        free(Cs);
        SparseFree(&Ccur);
        SparseFree(&product);
        SparseFree(&tmp1);
        SparseFree(&tmp2);
        SparseFree(&residual);
        SparseFree(&basisTmp);
        DenseAccumFree(&accum);
        FreeContext(&ctx);
        return 1;
    }

    int head = 0;
    int tail = 0;
    for (int i = 0; i < ctx.count; i++) { pred[i] = -1; predGen[i] = -1; constTerm[i] = 0; }
    visited[args.d] = 1;
    queue[tail++] = args.d;

    while (head < tail) {
        int cur = queue[head++];

        BuildKLElementSparse(&ctx, cur, &Ccur);

        for (int m = 0; m < multiplierCount; m++) {
            SparseClear(&product);
            /* multiply C(s) * C(cur) where s is either a simple generator
               or any group element depending on options */
            int mulIdx = args.useAllMultipliers ? m : ctx.genLehmer[m];
            MultiplyHeckeSparse(&ctx, mulIdx, &Ccur, &Cs[m], &product, &tmp1, &tmp2, &accum, NULL, 0);

            DecomposeToKLBasis(&ctx, &product, coeffs, &residual, &basisTmp);

            for (int x = 0; x < ctx.count; x++) {
                /* Require the KL-basis coefficient polynomial to have a nonzero
                   constant term (v^0). This avoids counting terms that only
                   appear with positive or negative v-powers which do not
                   witness the right-preorder relation at q=1. */
                if (PolyIsZero(&coeffs[x]) || coeffs[x].coeff[28] == 0) {
                    continue;
                }
                if (!visited[x]) {
                    visited[x] = 1;
                    pred[x] = cur;
                    predGen[x] = m;
                    constTerm[x] = coeffs[x].coeff[28];
                    queue[tail++] = x;
                }
            }
        }
    }

    printf("# Involutions d' with d' <=_R d in S_%d\n", args.n);
    printf("# d = %d\n", args.d);
    printf("index,perm\n");

    int involutionCount = 0;
    for (int x = 0; x < ctx.count; x++) {
        if (!visited[x] || !IsInvolution(args.n, x)) {
            continue;
        }

        char perm[8] = {0};
        IndexToPerm(args.n, x, perm);
        printf("%d,", x);
        for (int i = 0; i < args.n; i++) {
            printf("%d", (int)perm[i]);
            if (i + 1 < args.n) {
                printf(" ");
            }
        }
        printf("\n");
        involutionCount++;
    }

    fprintf(stderr, "Reached %d elements in <=_R-closure from d=%d.\n", tail, args.d);
    fprintf(stderr, "Involutions in that closure: %d.\n", involutionCount);

    /* Debug: print reconstruction paths for each involution found. */
    fprintf(stderr, "\nPaths (d <- ... <- target):\n");
    for (int x = 0; x < ctx.count; x++) {
        if (!visited[x] || !IsInvolution(args.n, x)) continue;
        int stackSize = 0;
        int *stack = malloc(sizeof(int) * ctx.count);
        int cur = x;
        while (cur != -1) {
            stack[stackSize++] = cur;
            if (cur == args.d) break;
            cur = pred[cur];
        }
        for (int j = stackSize - 1; j >= 0; j--) {
            fprintf(stderr, "%d", stack[j]);
            if (j) {
                int next = stack[j-1];
                fprintf(stderr, " -[g=%d,c=%d]-> ", predGen[next], constTerm[next]);
            }
        }
        fprintf(stderr, "\n");
        free(stack);
    }

    free(coeffs);
    free(queue);
    free(visited);

    for (int g = 0; g < ctx.genCount; g++) {
        SparseFree(&Cs[g]);
    }
    free(Cs);

    SparseFree(&Ccur);
    SparseFree(&product);
    SparseFree(&tmp1);
    SparseFree(&tmp2);
    SparseFree(&residual);
    SparseFree(&basisTmp);
    DenseAccumFree(&accum);
    FreeContext(&ctx);

    return 0;
}
