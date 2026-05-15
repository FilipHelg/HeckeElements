#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

#include <stdio.h>
#include <string.h>

typedef struct {
    SparseHecke *records;
    unsigned char *present;
    int count;
} DualCache;

static int LoadDualCache(int n, const char *path, DualCache *cache) {
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    
    char magic[8];
    fread(magic, 1, 8, f);
    if (strncmp(magic, "DKBULK1\0", 8) != 0) {
        fclose(f);
        return 0;
    }
    
    int version, fileN, count;
    fread(&version, sizeof(int), 1, f);
    fread(&fileN, sizeof(int), 1, f);
    fread(&count, sizeof(int), 1, f);
    
    if (fileN != n || count <= 0) {
        fclose(f);
        return 0;
    }
    
    cache->count = count;
    cache->present = calloc(count, sizeof(unsigned char));
    cache->records = calloc(count, sizeof(SparseHecke));
    
    if (!cache->present || !cache->records) {
        fclose(f);
        free(cache->present);
        free(cache->records);
        return 0;
    }
    
    for (int i = 0; i < count; i++) {
        SparseInit(&cache->records[i], count);
    }
    
    for (int i = 0; i < count; i++) {
        unsigned char present;
        fread(&present, 1, 1, f);
        cache->present[i] = present;
        
        if (!present) continue;
        
        int supportSize;
        fread(&supportSize, sizeof(int), 1, f);
        
        if (supportSize > 0) {
            cache->records[i].support = malloc(supportSize * sizeof(int));
            cache->records[i].coeff = malloc(supportSize * sizeof(Poly57));
            fread(cache->records[i].support, sizeof(int), supportSize, f);
            fread(cache->records[i].coeff, sizeof(Poly57), supportSize, f);
            cache->records[i].supportSize = supportSize;
        }
    }
    
    fclose(f);
    return 1;
}

int main(int argc, char *argv[]) {
    int n = 5;
    int count = 120;
    
    FastContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.n = n;
    ctx.count = count;
    ctx.maxLength = n * (n - 1) / 2;
    
    ctx.lengths = calloc(count, sizeof(int));
    for (int i = 0; i < count; i++) {
        ctx.lengths[i] = IndexToLength(n, i);
    }
    
    if (!BuildLengthOrder(&ctx) || !BuildGeneratorTables(&ctx) || !BuildExpressionTables(&ctx)) {
        fprintf(stderr, "Failed to build context\n");
        return 1;
    }
    
    DualCache dualCache;
    memset(&dualCache, 0, sizeof(dualCache));
    if (!LoadDualCache(n, "dual_kl_bulk_n5.bin", &dualCache)) {
        fprintf(stderr, "Failed to load dual cache\n");
        return 1;
    }
    
    printf("Loaded dual cache: %d involutions present\n\n", dualCache.count);
    
    /* Find a sample left cell with multiple elements in S5 */
    printf("Finding left cells in S5...\n");
    
    typedef struct {
        int index;
        unsigned char P[64];
    } Entry;
    
    Entry *entries = malloc(count * sizeof(Entry));
    for (int i = 0; i < count; i++) {
        char P[8][8] = {{0}};
        char Q[8][8] = {{0}};
        char shape[8] = {0};
        RSTableaux(n, i, P, Q, shape);
        entries[i].index = i;
        memcpy(entries[i].P, P, 64);
    }
    
    typedef int (*compare_fn)(const void*, const void*);
    qsort(entries, count, sizeof(Entry), (compare_fn)memcmp);
    
    /* Find first left cell with 2+ elements */
    for (int i = 0; i < count; ) {
        int j = i + 1;
        while (j < count && memcmp(entries[i].P, entries[j].P, 64) == 0) j++;
        
        if (j - i >= 2) {
            printf("Found left cell at position %d with %d elements:\n", i, j - i);
            int *cell_members = (int*)(entries + i);
            for (int k = i; k < j; k++) {
                printf("  x = %d\n", entries[k].index);
            }
            
            /* Try w = 0 (identity, always an involution) */
            int w = 0;
            if (!dualCache.present[w]) {
                printf("w=%d is not an involution\n", w);
                i = j;
                continue;
            }
            
            printf("\nTesting w=%d with this left cell:\n", w);
            
            /* Compute dkH_w * kH_x for x in this cell */
            DenseAccum accum;
            SparseHecke tmp1, tmp2, product, kx;
            
            DenseAccumInit(&accum, count);
            SparseInit(&tmp1, count);
            SparseInit(&tmp2, count);
            SparseInit(&product, count);
            SparseInit(&kx, count);
            
            for (int k = i; k < j && k < i + 3; k++) {
                int x = entries[k].index;
                printf("\n  x=%d:\n", x);
                
                /* Build kH_x */
                SparseClear(&kx);
                if (!BuildKLElementSparse(&ctx, x, &kx)) {
                    printf("    Failed to build kH_%d\n", x);
                    continue;
                }
                printf("    kH_%d support size: %d\n", x, kx.supportSize);
                
                /* Compute dkH_w * kH_x */
                SparseClear(&product);
                MultiplyHeckeSparse(&ctx, x, &dualCache.records[w], &kx, &product, &tmp1, &tmp2, &accum, NULL, 0);
                printf("    (dkH_%d * kH_%d) support size: %d\n", w, x, product.supportSize);
                
                if (product.supportSize > 0 && product.supportSize <= 5) {
                    printf("      Support: ");
                    for (int s = 0; s < product.supportSize; s++) {
                        printf("%d ", product.support[s]);
                    }
                    printf("\n");
                }
            }
            
            SparseFree(&kx);
            SparseFree(&product);
            SparseFree(&tmp1);
            SparseFree(&tmp2);
            DenseAccumFree(&accum);
            
            break;  /* Just analyze first multi-element cell */
        }
        
        i = j;
    }
    
    free(entries);
    FreeContext(&ctx);
    
    return 0;
}
