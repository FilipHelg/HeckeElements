#include <stdio.h>

int main() {
    printf("S_3 permutations:\n");
    printf("  012: identity\n");
    printf("  021: (12)\n");
    printf("  102: (01)\n");
    printf("  120: (012)\n");
    printf("  201: (021)\n");
    printf("  210: (01)(12)\n");
    
    printf("\nFor test 1: x=3 (perm 120), y=4 (perm 201)\n");
    printf("Computing H_120 * H_201 in Hecke algebra:\n");
    printf("This is left-multiply by H_120 on H_201\n");
    printf("In standard basis, this gives H_(120*201) + (v^-1-v)*H_(...)\n");
    printf("Where 120*201 means: apply 201 then 120 (right-to-left)\n");
    printf("  0 -> 2 (by 120) -> 1 (by 201) = 0 -> 1\n");
    printf("  1 -> 0 (by 120) -> 2 (by 201) = 1 -> 2\n");
    printf("  2 -> 1 (by 120) -> 0 (by 201) = 2 -> 0\n");
    printf("  Result perm: 120 (same!)\n");
    printf("\nWait, that means H_120 * H_201 starts with H_120 as first term!\n");
    printf("But that doesn't seem right for Hecke algebra multiplication...\n");
    
    return 0;
}
