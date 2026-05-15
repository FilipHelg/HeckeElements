#include <stdio.h>
#include "lehmer.h"

void PrintTableau(char tab[8][8], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%2d ", (int)tab[i][j]);
        }
        printf("\n");
    }
}

int main() {
    int n = 4;
    int indices[] = {1, 18, 3, 12, 6, 9};
    
    for (int idx = 0; idx < 6; idx++) {
        int i = indices[idx];
        char P[8][8] = {{0}}, Q[8][8] = {{0}}, shape[8] = {0};
        RSTableaux(n, i, P, Q, shape);
        
        printf("Index %d:\n", i);
        printf("P:\n");
        PrintTableau(P, n);
        printf("Q:\n");
        PrintTableau(Q, n);
        printf("\n");
    }
    return 0;
}
