#include <stdio.h>
#include "lehmer.h"

int main(void) {
    int n = 4;
    for (int x = 0; x < fac(n); x++) {
        char expr[40] = {0};
        int len = ReducedExpression(n, x, expr);

        int left_rev = 0;
        for (int i = len - 1; i >= 0; i--) {
            int s = fac(n - (int)expr[i]);
            left_rev = MultiplyIndex(n, s, left_rev);
        }

        int left_fwd = 0;
        for (int i = 0; i < len; i++) {
            int s = fac(n - (int)expr[i]);
            left_fwd = MultiplyIndex(n, s, left_fwd);
        }

        int right_fwd = 0;
        for (int i = 0; i < len; i++) {
            int s = fac(n - (int)expr[i]);
            right_fwd = MultiplyIndex(n, right_fwd, s);
        }

        if (left_rev != x) {
            printf("Mismatch at x=%d: left_rev=%d left_fwd=%d right_fwd=%d len=%d\n", x, left_rev, left_fwd, right_fwd, len);
            return 1;
        }
    }
    printf("All elements reconstructed by left multiplication with reverse expression order.\n");
    return 0;
}
