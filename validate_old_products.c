/* Validator: Load old checker CSV and verify which products are actually nonzero */

#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    unsigned char bytes[72];
} TableauKey;

typedef struct {
    int index;
    TableauKey key;
} TableauEntry;

typedef struct {
    SparseHecke *records;
    unsigned char *present;
    int count;
} DualCache;



static int LoadDualCache(const char *path, int count, DualCache *cache, int *loadedNOut) {
    FILE *in = fopen(path, "rb");
    if (!in) {
        fprintf(stderr, "Could not open dual KL cache %s\n", path);
        return 0;
    }

    typedef struct {
        char magic[8];
        uint32_t version;
        uint32_t n;
        uint32_t count;
        uint32_t involutionCount;
        uint32_t recordFormat;
        uint32_t reserved[3];
    } DualBulkHeader;

    DualBulkHeader header;
    if (fread(&header, sizeof(header), 1, in) != 1) {
        fprintf(stderr, "Failed to read dual KL header\n");
        fclose(in);
        return 0;
    }

    if ((int)header.count != count) {
        fprintf(stderr, "Dual KL cache count mismatch\n");
        fclose(in);
        return 0;
    }

    cache->records = (SparseHecke *)calloc((size_t)count, sizeof(SparseHecke));
    cache->present = (unsigned char *)calloc((size_t)count, sizeof(unsigned char));
    cache->count = count;
    if (!cache->records || !cache->present) {
        fclose(in);
        return 0;
    }

    for (uint32_t rec = 0; rec < header.involutionCount; rec++) {
        uint32_t x = 0;
        uint32_t supportCount = 0;
        if (!fread(&x, sizeof(x), 1, in) || !fread(&supportCount, sizeof(supportCount), 1, in)) {
            fprintf(stderr, "Failed reading record %u\n", rec);
            fclose(in);
            return 0;
        }

        if (x >= (uint32_t)count) {
            fprintf(stderr, "Record index %u out of range\n", x);
            fclose(in);
            return 0;
        }

        if (!SparseInit(&cache->records[x], count)) {
            fprintf(stderr, "Out of memory\n");
            fclose(in);
            return 0;
        }
        cache->present[x] = 1;

        for (uint32_t s = 0; s < supportCount; s++) {
            uint32_t z = 0;
            uint8_t termCount = 0;
            if (!fread(&z, sizeof(z), 1, in) || !fread(&termCount, sizeof(termCount), 1, in)) {
                fprintf(stderr, "Failed reading support\n");
                fclose(in);
                return 0;
            }

            Poly57 poly;
            PolyZero(&poly);
            for (uint8_t t = 0; t < termCount; t++) {
                int8_t exponent = 0;
                int32_t coeff = 0;
                if (!fread(&exponent, sizeof(exponent), 1, in) || 
                    !fread(&coeff, sizeof(coeff), 1, in)) {
                    fprintf(stderr, "Failed reading poly term\n");
                    fclose(in);
                    return 0;
                }
                int pos = 28 + (int)exponent;
                if (pos >= 0 && pos < 57) {
                    poly.coeff[pos] += (int)coeff;
                    if (poly.minPos > pos) poly.minPos = pos;
                    if (poly.maxPos < pos) poly.maxPos = pos;
                }
            }
            PolyNormalize(&poly);
            SparseAddPoly(&cache->records[x], (int)z, &poly);
        }
    }

    if (loadedNOut) *loadedNOut = (int)header.n;
    fclose(in);
    return 1;
}

static void DualCacheFree(DualCache *cache) {
    if (!cache) return;
    if (cache->records && cache->present) {
        for (int i = 0; i < cache->count; i++) {
            if (cache->present[i]) {
                SparseFree(&cache->records[i]);
            }
        }
    }
    free(cache->records);
    free(cache->present);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s n csv_file\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    const char *csvPath = argv[2];

    if (n < 3 || n > 8) {
        fprintf(stderr, "n must be 3-8\n");
        return 1;
    }

    FastContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = n;
    ctx.count = fac(n);
    ctx.maxLength = n * (n - 1) / 2;

    ctx.lengths = (int *)calloc((size_t)ctx.count, sizeof(int));
    if (!ctx.lengths) return 1;

    for (int i = 0; i < ctx.count; i++) {
        ctx.lengths[i] = IndexToLength(n, i);
    }

    if (!BuildLengthOrder(&ctx) || !BuildGeneratorTables(&ctx) || !BuildExpressionTables(&ctx)) {
        fprintf(stderr, "Failed to build context\n");
        FreeContext(&ctx);
        return 1;
    }

    if (n >= 4) {
        char klPath[32];
        snprintf(klPath, sizeof(klPath), "S%d.txt", n);
        if (!LoadNonTrivialData(&ctx, klPath)) {
            fprintf(stderr, "Failed to load KL data\n");
            FreeContext(&ctx);
            return 1;
        }
    }

    if (!MaybeBuildKLBasisCache(&ctx)) {
        fprintf(stderr, "KL cache build failed\n");
    }

    DualCache dualCache;
    memset(&dualCache, 0, sizeof(dualCache));
    char dualPath[128];
    snprintf(dualPath, sizeof(dualPath), "dual_kl_bulk_n%d.bin", n);
    int dualN = 0;
    if (!LoadDualCache(dualPath, ctx.count, &dualCache, &dualN)) {
        fprintf(stderr, "Failed to load dual cache\n");
        FreeContext(&ctx);
        return 1;
    }

    /* Allocate working space */
    SparseHecke kxTemp, tmp1, tmp2, product;
    DenseAccum accum;
    LeftMulCache leftCache;
    if (!SparseInit(&kxTemp, ctx.count) || !SparseInit(&tmp1, ctx.count) ||
        !SparseInit(&tmp2, ctx.count) || !SparseInit(&product, ctx.count) || 
        !DenseAccumInit(&accum, ctx.count) || !LeftMulCacheInit(&leftCache, ctx.count)) {
        fprintf(stderr, "Failed to allocate space\n");
        DualCacheFree(&dualCache);
        FreeContext(&ctx);
        return 1;
    }

    /* Load and validate CSV */
    FILE *csvIn = fopen(csvPath, "r");
    if (!csvIn) {
        fprintf(stderr, "Could not open %s\n", csvPath);
        DualCacheFree(&dualCache);
        FreeContext(&ctx);
        return 1;
    }

    char line[512];
    if (!fgets(line, sizeof(line), csvIn)) {
        fprintf(stderr, "CSV is empty\n");
        fclose(csvIn);
        DualCacheFree(&dualCache);
        FreeContext(&ctx);
        return 1;
    }

    int64_t totalRows = 0;
    int64_t zeroRows = 0;
    int64_t nonzeroRows = 0;

    printf("Validating %s:\n", csvPath);
    printf("Checking products for zero/nonzero status...\n\n");

    while (fgets(line, sizeof(line), csvIn)) {
        int w, x, y;
        if (sscanf(line, "%d,%d,%d", &w, &x, &y) != 3) continue;

        totalRows++;

        /* Compute product dKH_w * kH_x */
        SparseClear(&kxTemp);
        BuildKLElementSparse(&ctx, x, &kxTemp);
        SparseClear(&product);
        MultiplyHeckeSparse(&ctx, x, &dualCache.records[w], &kxTemp, &product, 
                           &tmp1, &tmp2, &accum, NULL, 0);

        if (product.supportSize <= 0) {
            zeroRows++;
        } else {
            nonzeroRows++;
            if (nonzeroRows <= 20) {  /* Print first 20 nonzero */
                printf("Nonzero: (%d, %d, %d)\n", w, x, y);
            }
        }
    }

    fclose(csvIn);

    printf("\n=== SUMMARY ===\n");
    printf("Total rows: %lld\n", (long long)totalRows);
    printf("Zero products: %lld (%.1f%%)\n", (long long)zeroRows, 100.0*zeroRows/totalRows);
    printf("Nonzero products: %lld (%.1f%%)\n", (long long)nonzeroRows, 100.0*nonzeroRows/totalRows);
    printf("\nIf nonzero > 0: Old code finds real matches\n");
    printf("If nonzero = 0: Old code only finds zero=zero matches (bug in old code!)\n");

    SparseFree(&kxTemp);
    SparseFree(&tmp1);
    SparseFree(&tmp2);
    SparseFree(&product);
    DenseAccumFree(&accum);
    LeftMulCacheFree(&leftCache);
    DualCacheFree(&dualCache);
    FreeContext(&ctx);
    return 0;
}
