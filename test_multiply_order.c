#include <stdio.h>
#include "lehmer.h"

void PrintPerm(char *perm, int n) {
    for (int i = 0; i < n; i++) printf("%d", perm[i]);
}

int main() {
    int n = 3;
    
    printf("Testing permutation composition order:\n\n");
    
    // Test case 1: x=1, y=2 (two simple transpositions in S_3)
    printf("Index 1 (perm 021 - transposition (1,2)):\n");
    char perm1[8];
    IndexToPerm(n, 1, perm1);
    printf("  Permutation: ");
    PrintPerm(perm1, n);
    printf("\n");
    
    printf("Index 2 (perm 102 - transposition (0,1)):\n");
    char perm2[8];
    IndexToPerm(n, 2, perm2);
    printf("  Permutation: ");
    PrintPerm(perm2, n);
    printf("\n");
    
    // What does current MultiplyIndex(n, 1, 2) produce?
    int comp12 = MultiplyIndex(n, 1, 2);
    char comp12_perm[8];
    IndexToPerm(n, comp12, comp12_perm);
    printf("MultiplyIndex(n=3, x=1, y=2) = index %d (perm ", comp12);
    PrintPerm(comp12_perm, n);
    printf(")\n");
    
    // Manual composition right-to-left (apply y=2 first, then x=1):
    // Apply perm2: (0,1,2) -> (1,0,2)
    // Apply perm1: (1,0,2) -> (1,2,0)
    // Result: 120 = index 3
    printf("Expected if right-to-left (y first, then x): perm 120 (index 3)\n");
    
    // Manual composition left-to-right (apply x=1 first, then y=2):
    // Apply perm1: (0,1,2) -> (0,2,1)
    // Apply perm2: (0,2,1) -> (2,0,1)
    // Result: 201 = index 4
    printf("Expected if left-to-right (x first, then y): perm 201 (index 4)\n");
    
    printf("\n");
    printf("Test 2: x=3 (120), y=4 (201) - from your original failing test\n");
    char perm3[8], perm4[8];
    IndexToPerm(n, 3, perm3);
    IndexToPerm(n, 4, perm4);
    printf("  Index 3: ");
    PrintPerm(perm3, n);
    printf("\n  Index 4: ");
    PrintPerm(perm4, n);
    printf("\n");
    
    int comp34 = MultiplyIndex(n, 3, 4);
    char comp34_perm[8];
    IndexToPerm(n, comp34, comp34_perm);
    printf("MultiplyIndex(n=3, x=3, y=4) = index %d (perm ", comp34);
    PrintPerm(comp34_perm, n);
    printf(")\n");
    
    // Right-to-left: apply 4 first (201), then 3 (120)
    // (0,1,2) -> (apply 201) -> (1,2,0) -> (apply 120) -> (2,0,1)
    // Wait, let me recalculate:
    // 201 means: 0->2, 1->0, 2->1
    // 120 means: 0->1, 1->2, 2->0
    // Apply 201: (0,1,2) positions get values (201)
    // Position 0 receives perm4[0]=2, Position 1 receives perm4[1]=0, Position 2 receives perm4[2]=1
    // Result: (2,0,1)
    // Now apply 120 to (2,0,1):
    // Position 0 receives perm3[2]=1, Position 1 receives perm3[0]=1... wait
    // perm3[i] tells us what value is at position i after applying perm3
    // If input is (2,0,1), applying perm3=(1,2,0) means:
    // Position 0 gets input[perm3[0]?]... no that's wrong
    // perm3 as a permutation means: perm3[i] is where input position i goes to
    // So if we have input (2,0,1):
    //   input[0]=2 goes to position perm3[0]=1
    //   input[1]=0 goes to position perm3[1]=2
    //   input[2]=1 goes to position perm3[2]=0
    // Result output: position 0 gets 1, position 1 gets 2, position 2 gets 0 = (1,2,0)
    // That's index 3!
    printf("Expected if right-to-left (y first, then x): should give index 3 (perm 120)\n");
    
    // Let me also compute left-to-right manually:
    // Apply 3 (120) first:
    // (0,1,2) -> (apply 120) -> (1,2,0)
    // Then apply 4 (201):
    // (1,2,0) -> (apply 201): position 0 gets 0, position 1 gets 1, position 2 gets 2... 
    // No wait, 201 means perm[i] = what ends up at position i after the permutation
    // If we apply 201=(2,0,1) to (1,2,0):
    //   input[0]=1 goes to perm[0]=2... that's confusing
    //
    // Let me use function notation instead. If perm=(2,0,1), that's the function f where f(0)=2, f(1)=0, f(2)=1
    // So f composed with identity: f(identity) = (f(0), f(1), f(2)) = (2,0,1) ✓
    // Now if we have current state (1,2,0) and apply function g=(1,2,0):
    // Result = (g(1), g(2), g(0)) = (2, 0, 1) = 201 = index 4
    printf("Expected if left-to-right (x first, then y): should give index 4 (perm 201)\n");
    
    return 0;
}
