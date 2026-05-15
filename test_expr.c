#include <stdio.h>
#include "lehmer.c"

int main() {
    printf("Testing ReducedExpression for S_3:\n");
    printf("Index 3 (permutation 120):\n");
    char expr3[8];
    int len3 = ReducedExpression(3, 3, expr3);
    printf("  Length: %d\n", len3);
    printf("  Generators (forward): ");
    for (int i = 0; i < len3; i++) printf("%d ", expr3[i]);
    printf("\n");
    printf("  Generators (reverse): ");
    for (int i = len3 - 1; i >= 0; i--) printf("%d ", expr3[i]);
    printf("\n");

    printf("\nIndex 4 (permutation 201):\n");
    char expr4[8];
    int len4 = ReducedExpression(3, 4, expr4);
    printf("  Length: %d\n", len4);
    printf("  Generators (forward): ");
    for (int i = 0; i < len4; i++) printf("%d ", expr4[i]);
    printf("\n");
    printf("  Generators (reverse): ");
    for (int i = len4 - 1; i >= 0; i--) printf("%d ", expr4[i]);
    printf("\n");

    printf("\nIndex 0 (identity):\n");
    char expr0[8];
    int len0 = ReducedExpression(3, 0, expr0);
    printf("  Length: %d\n", len0);
    printf("  Generators (forward): ");
    for (int i = 0; i < len0; i++) printf("%d ", expr0[i]);
    printf("\n");

    printf("\nIndex 1 (permutation 010):\n");
    char expr1[8];
    int len1 = ReducedExpression(3, 1, expr1);
    printf("  Length: %d\n", len1);
    printf("  Generators (forward): ");
    for (int i = 0; i < len1; i++) printf("%d ", expr1[i]);
    printf("\n");
    
    return 0;
}
