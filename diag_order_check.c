#define KL_INVOLUTION_CHECK_SKIP_INCLUDE
#define main kl_involution_check_embedded_main
#include "kl_involution_check.c"
#undef main

static int CompactFromSparseSorted(const SparseHecke *src, CompactHecke *dst) {
    return CompactFromSparse(src, dst);
}

static void PrintCompact(const char *label, const FastContext *ctx, const CompactHecke *h, const char *permTable) {
    printf("%s: ", label);
    if (h->supportSize <= 0) {
        printf("0\n");
        return;
    }
    for (int i = 0; i < h->supportSize; i++) {
        int idx = h->support[i];
        if (i > 0) printf(" + ");
        FprintPoly(stdout, &h->coeff[i]);
        printf("*C_{%s}", permTable + (size_t)idx * (size_t)(ctx->n + 1));
    }
    printf("\n");
}

int main(void) {
    int n = 4;
    int w = 7;
    int xs[2] = {1, 18};

    FastContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = n;
    ctx.count = fac(n);
    ctx.maxLength = n * (n - 1) / 2;

    ctx.lengths = (int *)calloc((size_t)ctx.count, sizeof(int));
    for (int i = 0; i < ctx.count; i++) ctx.lengths[i] = IndexToLengthSafeLocal(n, i);
    if (!BuildLengthOrder(&ctx) || !BuildGeneratorTables(&ctx) || !BuildExpressionTables(&ctx)) return 1;
    if (!LoadNonTrivialData(&ctx, "S4.txt")) return 1;
    (void)MaybeBuildKLBasisCache(&ctx);

    DualCache dualCache;
    memset(&dualCache, 0, sizeof(dualCache));
    int dualN = 0;
    if (!LoadDualCache("dual_kl_bulk_n4.bin", ctx.count, &dualCache, &dualN)) return 1;
    if (!dualCache.present[w]) return 1;

    SparseHecke kx, productLR, productRL, tmp1, tmp2;
    DenseAccum accum;
    if (!SparseInit(&kx, ctx.count) || !SparseInit(&productLR, ctx.count) || !SparseInit(&productRL, ctx.count) ||
        !SparseInit(&tmp1, ctx.count) || !SparseInit(&tmp2, ctx.count) || !DenseAccumInit(&accum, ctx.count)) return 1;

    char *permTable = BuildPermStringTable(n, ctx.count);
    for (int t = 0; t < 2; t++) {
        int x = xs[t];
        SparseClear(&kx);
        BuildKLElementSparse(&ctx, x, &kx);

        SparseClear(&productLR);
        MultiplyHeckeSparse(&ctx, x, &dualCache.records[w], &kx, &productLR, &tmp1, &tmp2, &accum, NULL, 0);

        SparseClear(&productRL);
        MultiplyHeckeSparse(&ctx, w, &kx, &dualCache.records[w], &productRL, &tmp1, &tmp2, &accum, NULL, 0);

        CompactHecke cLR, cRL;
        memset(&cLR, 0, sizeof(cLR));
        memset(&cRL, 0, sizeof(cRL));
        CompactFromSparseSorted(&productLR, &cLR);
        CompactFromSparseSorted(&productRL, &cRL);

        printf("x=%d\n", x);
        PrintCompact("dw*kx", &ctx, &cLR, permTable);
        PrintCompact("kx*dw", &ctx, &cRL, permTable);
        CompactFree(&cLR);
        CompactFree(&cRL);
        printf("\n");
    }

    SparseFree(&kx);
    SparseFree(&productLR);
    SparseFree(&productRL);
    SparseFree(&tmp1);
    SparseFree(&tmp2);
    DenseAccumFree(&accum);
    free(permTable);
    DualCacheFree(&dualCache);
    FreeContext(&ctx);
    return 0;
}
