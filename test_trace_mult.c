#include <stdio.h>
#include "lehmer.h"

int main() {
    int n = 3;
    
    printf("=== Manual trace of H_3 * H_4 ===\n\n");
    
    // Index 3 and 4
    char perm3[8], perm4[8];
    IndexToPerm(n, 3, perm3);
    IndexToPerm(n, 4, perm4);
    printf("Index 3 (perm ");
    for (int i = 0; i < n; i++) printf("%d", perm3[i]);
    printf("): length = %d\n", IndexToLength(n, 3));
    
    printf("Index 4 (perm ");
    for (int i = 0; i < n; i++) printf("%d", perm4[i]);
    printf("): length = %d\n", IndexToLength(n, 4));
    
    // Reduced expressions
    char expr3[28], expr4[28];
    int len3 = ReducedExpression(n, 3, expr3);
    int len4 = ReducedExpression(n, 4, expr4);
    printf("\nReduced expression for 3: ");
    for (int i = 0; i < len3; i++) printf("%d ", expr3[i]);
    printf("(1-based generators)\n");
    
    printf("Reduced expression for 4: ");
    for (int i = 0; i < len4; i++) printf("%d ", expr4[i]);
    printf("(1-based generators)\n");
    
    // What generator indices are these?
    printf("\nGenerator conversions (gen_idx = fac(n - generator_label)):\n");
    for (int i = 0; i < len3; i++) {
        int gen_label = expr3[i];
        long gen_idx = 1;
        for (int j = 0; j < n - gen_label; j++) gen_idx *= (j+1);
        printf("  expr3[%d] = %d -> fac(%d) = %ld\n", i, gen_label, n - gen_label, gen_idx);
    }
    
    printf("\nComputing H_3 * H_4 as: apply generators of 3 (reversed) to H_4\n");
    printf("We apply the reduced expression for 3 in reverse order: ");
    for (int i = len3-1; i >= 0; i--) printf("%d ", expr3[i]);
    printf("\n\n");
    
    // Simulate applying generators
    printf("Start with index 4 (perm 201)\n");
    int current = 4;
    
    for (int i = len3-1; i >= 0; i--) {
        int gen_label = expr3[i];
        long gen_idx = 1;
        for (int j = 0; j < n - gen_label; j++) gen_idx *= (j+1);
        
        printf("Apply generator %d (fac(%d) = index %ld):\n", gen_label, n - gen_label, gen_idx);
        printf("  leftAction[%d][%d] = MultiplyIndex(%ld, %d) = ", gen_label-1, current, gen_idx, current);
        
        int next = MultiplyIndex(n, (int)gen_idx, current);
        printf("%d\n", next);
        
        char perm_next[8];
        IndexToPerm(n, next, perm_next);
        printf("  Result: perm ");
        for (int k = 0; k < n; k++) printf("%d", perm_next[k]);
        printf(", length %d\n\n", IndexToLength(n, next));
        
        current = next;
    }
    
    printf("Final permutation result: index %d\n", current);
    char perm_final[8];
    IndexToPerm(n, current, perm_final);
    printf("Permutation: ");
    for (int i = 0; i < n; i++) printf("%d", perm_final[i]);
    printf("\n");
    
    return 0;
}
