/* Minimal test to compute one S5 product manually */
#ifndef KL_TRIPLE_CHECK_NONZERO_SKIP_INCLUDE
#define KL_TRIPLE_CHECK_NONZERO_SKIP_INCLUDE 1
#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int index;
    SparseHecke poly;
} InvolutionRecord;

static int LoadInvolutions(const char *path, int count, InvolutionRecord **outRecords, int *outCount) {
    FILE *in = fopen(path, "rb");
    if (!in) {
        fprintf(stderr, "Could not open %s\n", path);
        return 0;
    }
    
    char magic[8];
    uint32_t version, n, cnt, inv_cnt;
    fread(magic, 1, 8, in);
    fread(&version, 4, 1, in);
    fread(&n, 4, 1, in);
    fread(&cnt, 4, 1, in);
    fread(&inv_cnt, 4, 1, in);
    
    printf("Loaded header: n=%u, count=%u, involutions=%u\n", n, cnt, inv_cnt);
    
    unsigned char *present = (unsigned char *)malloc(cnt);
    fread(present, 1, cnt, in);
    
    InvolutionRecord *records = (InvolutionRecord *)calloc(inv_cnt, sizeof(InvolutionRecord));
    int idx = 0;
    
    for (uint32_t i = 0; i < inv_cnt; i++) {
        uint32_t x, support_count;
        fread(&x, 4, 1, in);
        fread(&support_count, 4, 1, in);
        
        records[idx].index = x;
        SparseInit(&records[idx].poly, count);
        
        for (uint32_t s = 0; s < support_count; s++) {
            uint32_t z;
            uint8_t term_count;
            fread(&z, 4, 1, in);
            fread(&term_count, 1, 1, in);
            
            Poly57 poly;
            PolyZero(&poly);
            
            for (uint8_t t = 0; t < term_count; t++) {
                int8_t exp;
                int32_t coeff;
                fread(&exp, 1, 1, in);
                fread(&coeff, 4, 1, in);
                
                int pos = 28 + (int)exp;
                if (pos >= 0 && pos < 57) {
                    poly.coeff[pos] += coeff;
                    if (poly.minPos > pos) poly.minPos = pos;
                    if (poly.maxPos < pos) poly.maxPos = pos;
                }
            }
            PolyNormalize(&poly);
            SparseAddPoly(&records[idx].poly, (int)z, &poly);
        }
        
        if (idx < 5) {
            printf("  Record %d: index=%u, support_count=%u -> poly support=%d\n",
                   idx, x, support_count, records[idx].poly.supportSize);
        }
        
        idx++;
    }
    
    free(present);
    fclose(in);
    
    *outRecords = records;
    *outCount = idx;
    return 1;
}

int main(int argc, char *argv[]) {
    FastContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = 5;
    ctx.count = fac(5);
    ctx.maxLength = 5 * 4 / 2;
    
    ctx.lengths = (int *)calloc((size_t)ctx.count, sizeof(int));
    for (int i = 0; i < ctx.count; i++) {
        ctx.lengths[i] = IndexToLength(5, i);
    }
    
    if (!BuildLengthOrder(&ctx) || !BuildGeneratorTables(&ctx) || !BuildExpressionTables(&ctx)) {
        fprintf(stderr, "Failed to build context\n");
        return 1;
    }
    
    if (!LoadNonTrivialData(&ctx, "S5.txt")) {
        fprintf(stderr, "Failed to load S5.txt\n");
        return 1;
    }
    
    // Load involutions
    InvolutionRecord *involutions = NULL;
    int inv_count = 0;
    if (!LoadInvolutions("dual_kl_bulk_n5.bin", ctx.count, &involutions, &inv_count)) {
        fprintf(stderr, "Failed to load involutions\n");
        return 1;
    }
    
    printf("\nTesting products for first involution (w=%d):\n", involutions[0].index);
    int w = involutions[0].index;
    
    // Take first two elements in same left cell
    int x = 0, y = 1;  
    
    printf("Computing dKH_%d * KH_%d and dKH_%d * KH_%d\n", w, x, w, y);
    
    // Build KL elements
    SparseHecke kx, ky, product1, product2;
    SparseHecke tmp1, tmp2;
    DenseAccum accum;
    memset(&accum, 0, sizeof(accum));
    DenseAccumInit(&accum, ctx.count);
    
    SparseInit(&kx, ctx.count);
    SparseInit(&ky, ctx.count);
    SparseInit(&product1, ctx.count);
    SparseInit(&product2, ctx.count);
    SparseInit(&tmp1, ctx.count);
    SparseInit(&tmp2, ctx.count);
    
    BuildKLElementSparse(&ctx, x, &kx);
    BuildKLElementSparse(&ctx, y, &ky);
    
    printf("KH_%d support: %d\n", x, kx.supportSize);
    printf("KH_%d support: %d\n", y, ky.supportSize);
    printf("dKH_%d support: %d\n", w, involutions[0].poly.supportSize);
    
    // Compute products
    MultiplyHeckeSparse(&ctx, x, &involutions[0].poly, &kx, &product1, &tmp1, &tmp2, &accum, NULL, 0);
    MultiplyHeckeSparse(&ctx, y, &involutions[0].poly, &ky, &product2, &tmp1, &tmp2, &accum, NULL, 0);
    
    printf("\ndKH_%d * KH_%d support: %d\n", w, x, product1.supportSize);
    printf("dKH_%d * KH_%d support: %d\n", w, y, product2.supportSize);
    
    if (product1.supportSize > 0) {
        printf("First product nonzero terms (first 5): ");
        for (int i = 0; i < 5 && i < product1.supportSize; i++) {
            int idx = product1.support[i];
            printf("[%d]", idx);
        }
        printf("\n");
    }
    
    DenseAccumFree(&accum);
    SparseFree(&kx);
    SparseFree(&ky);
    SparseFree(&product1);
    SparseFree(&product2);
    SparseFree(&tmp1);
    SparseFree(&tmp2);
    FreeContext(&ctx);
    
    return 0;
}
