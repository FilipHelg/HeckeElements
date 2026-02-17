#include "lehmer.h"

typedef struct{
    int n; // Dimension
    int x; // Lehmer index
    char code[8]; // Lehmer code [Hard-coded max length 8]
    char perm[8]; // Permutation result
    char reducedExpression[28]; // Expression in generators [Hard-coded max length 28]
} perm_t;

int fac(int n){
    if(n == 1) return 1;
    return n * fac(n-1);
}

// Convert integer x (0 <= x < n!) to its Lehmer code.
// Optimization: the last Lehmer digit is always 0
char* LehmerCode(int n, int x, char code[]){
    // produce n-1 numeric digits (0..)
    for(int i = 0; i < n - 1; i++){
        int fact = 1;
        for(int j = 1; j <= n - 1 - i; j++) fact *= j; // (n-1-i)!
        int digit = x / fact;
        x = x % fact;
        code[i] = (char)digit; // store numeric digit (not ASCII)
    }
    // last slot is always zero
    code[n - 1] = 0;
    return code;
}

int CodeSum(int n, char code[]){
    int sum = 0;
    for (int i = 0; i < n; i++){
        sum += code[i];
    }
    return sum;
}

int CodeToIndex(int n, char code[]){
    int sum = 0;
    for(int i = 0; i < n; i++){
        sum += code[i] * (int)tgamma(n - i);
    }
    return sum;
}

int IndexToLength(int n, int x){
    char code[8];
    return(CodeSum(n, LehmerCode(n, x, code)));
}

int MultiplyIndex(int n, int x, int y){
    char permx[8], permy[8], permp[8];
    IndexToPerm(n, x, permx);
    IndexToPerm(n, y, permy);
    for(int i = 0; i < n; i++){
        permp[i] = permy[(int)permx[i]];
    }
    return PermToIndex(n, permp);
}

int PermToIndex(int n, char* perm){
    char code[8];
    for(int i = 0; i < n; i++){
        int sum = 0;
        for(int j = i + 1; j < n; j++) if(perm[j] < perm[i]) sum++;
        code[i] = sum;
    }
    return CodeToIndex(n, code);
}

void IndexToPerm(int n, int x, char* perm){
    char code[8];
    LehmerCode(n, x, code);
    char available[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    for(int i = 0; i < n; i++){
        int j = code[i];
        while(1){
            perm[i] = available[j];
            for(int k = j; k < n-1-i; k++) available[k] = available[k+1];
            break;
            j++;
        }
    }
    return;
}

int FirstGeneratorIndex(int n, int x){
    int s_i = x/(int)tgamma(n);
    return s_i;
}

// Produce a reduced expression by moving units right; return length.
int ReducedExpression(int n, int x, char expression[]){
    char code[8]; // [Hard-coded max length 8]
    int counter = 0;
    LehmerCode(n, x, code);
    while(CodeSum(n, code) != 0){
        for(int i = n - 2; i >= 0; i--){
            if (code[i] > 0){
                expression[counter++] = (char)(i + 1);
                // Move units to the right and reduce total sum by 1
                code[i+1] = code[i] - 1;
                code[i] = 0;
                break;
            }
        }
    }
    return counter;
}
