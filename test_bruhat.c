#include <stdio.h>
#include <stdlib.h>
#include "lehmer.h"
#include "laurent.h"

// Correct Bruhat order using rank matrix criterion:
// x <= y iff for all 1 <= i, j <= n: r_x(i,j) <= r_y(i,j)
// where r_w(i,j) = #{k <= i : w(k) >= j}
int BruhatSmallerCorrect(int n, int x, int y) {
    char permx[8], permy[8];
    IndexToPerm(n, x, permx);
    IndexToPerm(n, y, permy);
    
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            // Count #{k <= i : perm[k] >= j}
            int rx = 0, ry = 0;
            for(int k = 0; k <= i; k++) {
                if(permx[k] >= j) rx++;
                if(permy[k] >= j) ry++;
            }
            if(rx > ry) return 0;
        }
    }
    return 1;
}

// FindR using CORRECT Bruhat order
Laurent_t FindR_Correct(int n, int w1, int w2){
    Laurent_t result = ZeroInitializeLaurent();
    if(w1 == w2){
        result.coeff[28] = 1;
        return result;
    }else if(BruhatSmallerCorrect(n, w1, w2)){  // Use correct Bruhat!
        char* exp = (char*)calloc(28, sizeof(char));
        int len = ReducedExpression(n, w2, exp);
        int s = fac(n - exp[0]);
        int sw1 = MultiplyIndex(n, s, w1);
        int sw2 = MultiplyIndex(n, s, w2);
        
        if(IndexToLength(n, sw1) < IndexToLength(n, w1)){
            free(exp);
            return FindR_Correct(n, sw1, sw2);
        }else{
            Laurent_t vminus1 = ZeroInitializeLaurent();
            Laurent_t v = ZeroInitializeLaurent();
            v.coeff[29] = 1;
            vminus1.coeff[28] = -1;
            vminus1.coeff[29] = 1;
            result = SumLaurent(MultiplyLaurent(v, FindR_Correct(n, sw1, sw2)), MultiplyLaurent(vminus1, FindR_Correct(n, w1, sw2)));
            free(exp);
            return result;
        }
    }else{
        return result;
    }
}

int main() {
    printf("=== Compare BruhatSmaller vs correct rank criterion ===\n");
    printf("Looking for discrepancies in S_4...\n\n");
    
    int discrepancies = 0;
    for(int x = 0; x < 24; x++) {
        for(int y = 0; y < 24; y++) {
            int subword = BruhatSmaller(4, x, y);
            int correct = BruhatSmallerCorrect(4, x, y);
            if(subword != correct) {
                char px[8], py[8];
                IndexToPerm(4, x, px);
                IndexToPerm(4, y, py);
                printf("DISCREPANCY: x=%d [%d%d%d%d], y=%d [%d%d%d%d]\n",
                       x, px[0], px[1], px[2], px[3],
                       y, py[0], py[1], py[2], py[3]);
                printf("  Subword check: %d, Correct: %d\n", subword, correct);
                discrepancies++;
            }
        }
    }
    
    if(discrepancies == 0) {
        printf("No discrepancies found - BruhatSmaller is correct for S_4!\n");
    } else {
        printf("\nTotal discrepancies: %d\n", discrepancies);
    }
    
    // Check interval [11, 23] with correct criterion
    printf("\n=== Interval [11, 23] using correct Bruhat criterion ===\n");
    printf("Elements z with 11 <= z <= 23:\n");
    for(int z = 0; z < 24; z++) {
        int le_23 = BruhatSmallerCorrect(4, z, 23);
        int ge_11 = BruhatSmallerCorrect(4, 11, z);
        if(le_23 && ge_11) {
            char pz[8];
            IndexToPerm(4, z, pz);
            char exp[28] = {0};
            int len = ReducedExpression(4, z, exp);
            printf("  %d: [%d,%d,%d,%d] len=%d expr=", z, pz[0], pz[1], pz[2], pz[3], len);
            for(int i=0; i<len; i++) printf("s%d", exp[i]);
            printf("\n");
        }
    }
    
    // Compute P(11, 23) with CORRECTED interval
    printf("\n=== Compute P(11, 23) with correct interval ===\n");
    printf("Interval: {11, 17, 21, 23}\n");
    printf("Sum = R(11,17)*P(17,23) + R(11,21)*P(21,23) + R(11,23)*P(23,23)\n");
    
    // Note: P(17,23), P(21,23), P(23,23) = 1 for adjacent elements
    printf("\nR(11, 17) = "); DisplayLaurentPoly(FindR_Correct(4, 11, 17)); printf("\n");
    printf("R(11, 21) = "); DisplayLaurentPoly(FindR_Correct(4, 11, 21)); printf(" (WAS ZERO, NOW CORRECT)\n");
    printf("R(11, 23) = "); DisplayLaurentPoly(FindR_Correct(4, 11, 23)); printf("\n");
    
    // Sum them up
    Laurent_t sum = ZeroInitializeLaurent();
    sum = SumLaurent(sum, FindR_Correct(4, 11, 17));
    sum = SumLaurent(sum, FindR_Correct(4, 11, 21));
    sum = SumLaurent(sum, FindR_Correct(4, 11, 23));
    printf("\nTotal sum = "); DisplayLaurentPoly(sum); printf("\n");
    printf("Expected for P=1: v^2 - 1 = -1 + v^2\n");
    
    return 0;
}
