#include <stdio.h>
#include "lehmer.h"

int main() {
    int n = 3;
    
    printf("Verifying generator setup:\n");
    printf("n=%d, genCount=%d\n", n, n-1);
    
    // Compute what generators should be
    for (int g = 0; g < n-1; g++) {
        long gen_lehmer = 1;
        for (int i = 0; i < (n - 1 - g); i++) gen_lehmer *= (i + 1);
        printf("Generator %d: fac(%d) = %ld\n", g, n-1-g, gen_lehmer);
    }
    
    printf("\nComputing leftAction values manually:\n");
    printf("leftAction[g][w] = MultiplyIndex(genLehmer[g], w)\n\n");
    
    for (int g = 0; g < n-1; g++) {
        // Compute generator index
        long gen_idx = 1;
        for (int i = 0; i < (n - 1 - g); i++) gen_idx *= (i + 1);
        
        printf("Generator %d (genLehmer=%ld):\n", g, gen_idx);
        printf("  w=0: %d\n", MultiplyIndex(n, (int)gen_idx, 0));
        printf("  w=1: %d\n", MultiplyIndex(n, (int)gen_idx, 1));
        printf("  w=2: %d\n", MultiplyIndex(n, (int)gen_idx, 2));
        printf("  w=3: %d\n", MultiplyIndex(n, (int)gen_idx, 3));
        printf("  w=4: %d\n", MultiplyIndex(n, (int)gen_idx, 4));
        printf("  w=5: %d\n", MultiplyIndex(n, (int)gen_idx, 5));
    }
    
    return 0;
}
