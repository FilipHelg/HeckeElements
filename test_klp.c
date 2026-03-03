#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lehmer.h"
#include "laurent.h"

Laurent_t FindR(int n, int w1, int w2);
Laurent_t FindKLP(int n, int x, int y);

int main(){
    printf("=== Testing small cases ===\n");
    
    // Test P(0,1) in S_3 - should be 1
    printf("\nFindKLP(3, 0, 1):\n");
    Laurent_t p1 = FindKLP(3, 0, 1);
    printf("Result: "); DisplayLaurentPoly(p1); printf("\n");
    
    // Test lengths to understand the structure
    printf("\nLengths in S_3:\n");
    for(int i = 0; i < 6; i++){
        printf("l(%d) = %d\n", i, IndexToLength(3, i));
    }
    
    // Test what ElementsBetween2 returns
    printf("\nElementsBetween2(3, 0, 1):\n");
    int* elts = ElementsBetween2(3, 0, 1);
    printf("Count: %d\n", elts[0]);
    for(int i = 1; i <= elts[0]; i++) printf("Element %d: %d\n", i, elts[i]);
    
    // Test R-polynomials
    printf("\nR-polynomials:\n");
    printf("R(0,0): "); DisplayLaurentPoly(FindR(3, 0, 0)); printf("\n");
    printf("R(0,1): "); DisplayLaurentPoly(FindR(3, 0, 1)); printf("\n");
    printf("R(0,2): "); DisplayLaurentPoly(FindR(3, 0, 2)); printf("\n");
    
    return 0;
}
