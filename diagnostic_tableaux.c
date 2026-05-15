/* Diagnostic: check tableau key generation for S5 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

static void PrintTableau(char P[8][8], int n) {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", P[i][j]);
        }
        printf("| ");
    }
}

int main(int argc, char *argv[]) {
    int n = 5;
    int count = 120;  // 5!
    
    printf("S5 tableaux and their grouping:\n\n");
    
    /* Group by tableau key */
    typedef struct {
        int index;
        char P[8][8];
        char shape[8];
    } Entry;
    
    Entry *entries = malloc(count * sizeof(Entry));
    
    for (int i = 0; i < count; i++) {
        char P[8][8] = {{0}};
        char Q[8][8] = {{0}};
        char shape[8] = {0};
        
        RSTableaux(n, i, P, Q, shape);
        
        entries[i].index = i;
        memcpy(entries[i].P, P, sizeof(P));
        memcpy(entries[i].shape, shape, sizeof(shape));
    }
    
    /* Sort by tableau key (P + shape) */
    typedef struct {
        unsigned char bytes[72];
    } TableauKey;
    
    TableauKey *keys = malloc(count * sizeof(TableauKey));
    for (int i = 0; i < count; i++) {
        memset(keys[i].bytes, 0, sizeof(keys[i].bytes));
        memcpy(keys[i].bytes, entries[i].P, 64);
        memcpy(keys[i].bytes + 64, entries[i].shape, 8);
    }
    
    /* Count unique keys */
    int uniqueKeys = 1;
    for (int i = 1; i < count; i++) {
        if (memcmp(keys[i].bytes, keys[i-1].bytes, sizeof(TableauKey)) != 0) {
            uniqueKeys++;
        }
    }
    
    printf("Total permutations: %d\n", count);
    printf("Unique left cells (tableau keys): %d\n", uniqueKeys);
    
    /* Print grouping */
    printf("\nLeft cell grouping:\n");
    int cellNum = 0;
    for (int i = 0; i < count; ) {
        int j = i + 1;
        while (j < count && memcmp(keys[i].bytes, keys[j].bytes, sizeof(TableauKey)) == 0) {
            j++;
        }
        int cellSize = j - i;
        
        if (cellSize > 1) {
            printf("Cell %d (size %d): ", cellNum++, cellSize);
            for (int k = i; k < j; k++) {
                printf("%d ", entries[k].index);
            }
            printf("\n");
        }
        
        i = j;
    }
    
    free(keys);
    free(entries);
    return 0;
}
