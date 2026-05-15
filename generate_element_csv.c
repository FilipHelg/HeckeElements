#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lehmer.h"

// Format a permutation as one-line notation
void FormatOneLine(int n, char perm[8], char *buffer, int bufsize) {
    int pos = 0;
    for (int i = 0; i < n; i++) {
        pos += snprintf(buffer + pos, bufsize - pos, "%d", (int)perm[i]);
    }
}

// Format reduced expression as a string (generator indices)
void FormatReducedExpression(int n, int idx, char *buffer, int bufsize) {
    char expression[28] = {0};
    int len = ReducedExpression(n, idx, expression);
    
    int pos = 0;
    for (int i = 0; i < len; i++) {
        if (i > 0) pos += snprintf(buffer + pos, bufsize - pos, ",");
        pos += snprintf(buffer + pos, bufsize - pos, "%d", (int)expression[i]);
    }
    if (pos == 0) {
        snprintf(buffer, bufsize, "empty");
    }
}

// Format a tableau as a string (row-by-row, comma-separated)
void FormatTableau(char tab[8][8], int n, char rowLen[8], char *buffer, int bufsize) {
    int pos = 0;
    int rows = 0;
    
    // Count rows
    for (int i = 0; i < n; i++) {
        if (rowLen[i] > 0) rows++;
        else break;
    }
    
    if (rows == 0) {
        snprintf(buffer, bufsize, "empty");
        return;
    }
    
    for (int i = 0; i < rows; i++) {
        if (i > 0) pos += snprintf(buffer + pos, bufsize - pos, "|");
        for (int j = 0; j < rowLen[i]; j++) {
            if (j > 0) pos += snprintf(buffer + pos, bufsize - pos, ",");
            pos += snprintf(buffer + pos, bufsize - pos, "%d", (int)tab[i][j]);
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s n\n", argv[0]);
        fprintf(stderr, "Generate CSV of all elements in S_n with:\n");
        fprintf(stderr, "  Lehmer index, one-line notation, reduced expression, P-tableau, Q-tableau\n");
        fprintf(stderr, "n must be in {3,4,5,6,7,8,9}\n");
        return 1;
    }

    int n = atoi(argv[1]);
    if (n < 3 || n > 9) {
        fprintf(stderr, "Error: n must be between 3 and 9.\n");
        return 1;
    }

    int count = fac(n);

    // Print CSV header
    printf("lehmer_index,one_line_notation,reduced_expression,p_tableau,q_tableau\n");

    // Generate data for each element
    for (int idx = 0; idx < count; idx++) {
        char perm[8] = {0};
        char P[8][8] = {{0}};
        char Q[8][8] = {{0}};
        char shape[8] = {0};
        char rowLen[8] = {0};

        // Get permutation
        IndexToPerm(n, idx, perm);

        // Format one-line notation
        char oneline[32];
        FormatOneLine(n, perm, oneline, sizeof(oneline));

        // Format reduced expression
        char expr_str[256];
        FormatReducedExpression(n, idx, expr_str, sizeof(expr_str));

        // Get tableaux
        int rows = RSTableaux(n, idx, P, Q, shape);
        for (int i = 0; i < rows; i++) {
            rowLen[i] = shape[i];
        }

        // Format tableaux
        char p_str[256], q_str[256];
        FormatTableau(P, n, rowLen, p_str, sizeof(p_str));
        FormatTableau(Q, n, rowLen, q_str, sizeof(q_str));

        // Print CSV row (escape quotes if any appear)
        printf("%d,%s,\"%s\",\"%s\",\"%s\"\n", idx, oneline, expr_str, p_str, q_str);
    }

    return 0;
}
