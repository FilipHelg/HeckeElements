#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

#define KL_INVOLUTION_CHECK_SKIP_INCLUDE
#define main kl_involution_check_embedded_main
#include "kl_involution_check.c"
#undef main

static int CompactFromSparseSorted(const SparseHecke *src, CompactHecke *dst) {
    return CompactFromSparse(src, dst);
}

static int SameQTableau(int n, int a, int b) {
    char P1[8][8] = {{0}}, Q1[8][8] = {{0}}, shape1[8] = {0};
    char P2[8][8] = {{0}}, Q2[8][8] = {{0}}, shape2[8] = {0};
    (void)RSTableaux(n, a, P1, Q1, shape1);
    (void)RSTableaux(n, b, P2, Q2, shape2);
    return memcmp(Q1, Q2, sizeof(Q1)) == 0;
}

static int SamePTableau(int n, int a, int b) {
    char P1[8][8] = {{0}}, Q1[8][8] = {{0}}, shape1[8] = {0};
    char P2[8][8] = {{0}}, Q2[8][8] = {{0}}, shape2[8] = {0};
    (void)RSTableaux(n, a, P1, Q1, shape1);
    (void)RSTableaux(n, b, P2, Q2, shape2);
    return memcmp(P1, P2, sizeof(P1)) == 0;
}

int main(void) {
    const int n = 5;
    const int w = 7;
    const int pairs[4][2] = {{1, 9}, {4, 8}, {6, 18}, {48, 96}};

    FastContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = n;
    ctx.count = fac(n);
    ctx.maxLength = n * (n - 1) / 2;

    ctx.lengths = (int *)calloc((size_t)ctx.count, sizeof(int));
    if (!ctx.lengths) return 1;
    for (int i = 0; i < ctx.count; i++) ctx.lengths[i] = IndexToLengthSafeLocal(n, i);

    if (!BuildLengthOrder(&ctx) || !BuildGeneratorTables(&ctx) || !BuildExpressionTables(&ctx)) return 1;
    if (!LoadNonTrivialData(&ctx, "S5.txt")) return 1;
    (void)MaybeBuildKLBasisCache(&ctx);

    DualCache dualCache;
    memset(&dualCache, 0, sizeof(dualCache));
    int dualN = 0;
    if (!LoadDualCache("dual_kl_bulk_n5.bin", ctx.count, &dualCache, &dualN)) return 1;
    if (dualN != n || !dualCache.present[w]) return 1;

    SparseHecke kx, ky, kw, prodLx, prodLy, prodRx, prodRy, prodKwLx, prodKwLy, prodKwRx, prodKwRy, tmp1, tmp2;
    DenseAccum accum;
    if (!SparseInit(&kx, ctx.count) || !SparseInit(&ky, ctx.count) || !SparseInit(&kw, ctx.count) ||
        !SparseInit(&prodLx, ctx.count) || !SparseInit(&prodLy, ctx.count) ||
        !SparseInit(&prodRx, ctx.count) || !SparseInit(&prodRy, ctx.count) ||
        !SparseInit(&prodKwLx, ctx.count) || !SparseInit(&prodKwLy, ctx.count) ||
        !SparseInit(&prodKwRx, ctx.count) || !SparseInit(&prodKwRy, ctx.count) ||
        !SparseInit(&tmp1, ctx.count) || !SparseInit(&tmp2, ctx.count) ||
        !DenseAccumInit(&accum, ctx.count)) return 1;

    SparseClear(&kw);
    BuildKLElementSparse(&ctx, w, &kw);

    printf("n=%d, w=%d\n", n, w);
    for (int p = 0; p < 4; p++) {
        int x = pairs[p][0];
        int y = pairs[p][1];

        SparseClear(&kx);
        SparseClear(&ky);
        BuildKLElementSparse(&ctx, x, &kx);
        BuildKLElementSparse(&ctx, y, &ky);

        SparseClear(&prodLx);
        SparseClear(&prodLy);
        SparseClear(&prodRx);
        SparseClear(&prodRy);
        SparseClear(&prodKwLx);
        SparseClear(&prodKwLy);
        SparseClear(&prodKwRx);
        SparseClear(&prodKwRy);

        /* "Left" order used by non-right variant: dw * kx */
        MultiplyHeckeSparse(&ctx, x, &dualCache.records[w], &kx, &prodLx, &tmp1, &tmp2, &accum, NULL, 0);
        MultiplyHeckeSparse(&ctx, y, &dualCache.records[w], &ky, &prodLy, &tmp1, &tmp2, &accum, NULL, 0);

        /* "Right" order used by right variant: kx * dw */
        MultiplyHeckeSparse(&ctx, w, &kx, &dualCache.records[w], &prodRx, &tmp1, &tmp2, &accum, NULL, 0);
        MultiplyHeckeSparse(&ctx, w, &ky, &dualCache.records[w], &prodRy, &tmp1, &tmp2, &accum, NULL, 0);

        /* Same checks but with kH_w in place of dKH_w */
        MultiplyHeckeSparse(&ctx, x, &kw, &kx, &prodKwLx, &tmp1, &tmp2, &accum, NULL, 0);
        MultiplyHeckeSparse(&ctx, y, &kw, &ky, &prodKwLy, &tmp1, &tmp2, &accum, NULL, 0);
        MultiplyHeckeSparse(&ctx, w, &kx, &kw, &prodKwRx, &tmp1, &tmp2, &accum, NULL, 0);
        MultiplyHeckeSparse(&ctx, w, &ky, &kw, &prodKwRy, &tmp1, &tmp2, &accum, NULL, 0);

        CompactHecke cLx, cLy, cRx, cRy, cKwLx, cKwLy, cKwRx, cKwRy;
        memset(&cLx, 0, sizeof(cLx));
        memset(&cLy, 0, sizeof(cLy));
        memset(&cRx, 0, sizeof(cRx));
        memset(&cRy, 0, sizeof(cRy));
        memset(&cKwLx, 0, sizeof(cKwLx));
        memset(&cKwLy, 0, sizeof(cKwLy));
        memset(&cKwRx, 0, sizeof(cKwRx));
        memset(&cKwRy, 0, sizeof(cKwRy));

        if (!CompactFromSparseSorted(&prodLx, &cLx) || !CompactFromSparseSorted(&prodLy, &cLy) ||
            !CompactFromSparseSorted(&prodRx, &cRx) || !CompactFromSparseSorted(&prodRy, &cRy) ||
            !CompactFromSparseSorted(&prodKwLx, &cKwLx) || !CompactFromSparseSorted(&prodKwLy, &cKwLy) ||
            !CompactFromSparseSorted(&prodKwRx, &cKwRx) || !CompactFromSparseSorted(&prodKwRy, &cKwRy)) {
            printf("pair (%d,%d): compact failed\n", x, y);
            continue;
        }

        printf("pair (%d,%d): Q=%s P=%s\n",
               x, y,
               SameQTableau(n, x, y) ? "same" : "diff",
               SamePTableau(n, x, y) ? "same" : "diff");

        printf("  dw*kx vs dw*ky: equal=%s, sizes=(%d,%d)\n",
               CompactEquals(&cLx, &cLy) ? "yes" : "no",
               cLx.supportSize, cLy.supportSize);

        printf("  kx*dw vs ky*dw: equal=%s, sizes=(%d,%d)\n",
               CompactEquals(&cRx, &cRy) ? "yes" : "no",
               cRx.supportSize, cRy.supportSize);

         printf("  kw*kx vs kw*ky: equal=%s, sizes=(%d,%d)\n",
             CompactEquals(&cKwLx, &cKwLy) ? "yes" : "no",
             cKwLx.supportSize, cKwLy.supportSize);

         printf("  kx*kw vs ky*kw: equal=%s, sizes=(%d,%d)\n",
             CompactEquals(&cKwRx, &cKwRy) ? "yes" : "no",
             cKwRx.supportSize, cKwRy.supportSize);

        CompactFree(&cLx);
        CompactFree(&cLy);
        CompactFree(&cRx);
        CompactFree(&cRy);
        CompactFree(&cKwLx);
        CompactFree(&cKwLy);
        CompactFree(&cKwRx);
        CompactFree(&cKwRy);
    }

    SparseFree(&kx);
    SparseFree(&ky);
    SparseFree(&kw);
    SparseFree(&prodLx);
    SparseFree(&prodLy);
    SparseFree(&prodRx);
    SparseFree(&prodRy);
    SparseFree(&prodKwLx);
    SparseFree(&prodKwLy);
    SparseFree(&prodKwRx);
    SparseFree(&prodKwRy);
    SparseFree(&tmp1);
    SparseFree(&tmp2);
    DenseAccumFree(&accum);
    DualCacheFree(&dualCache);
    FreeContext(&ctx);
    return 0;
}