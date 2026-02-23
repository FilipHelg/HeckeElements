#include "lehmer.h"

// Unused struct. Elements are now instead indexed as integers by their Lehmer index.
typedef struct{
    int n; // Dimension
    int x; // Lehmer index
    char code[8]; // Lehmer code [Hard-coded max length 8]
    char perm[8]; // Permutation result
    char reducedExpression[28]; // Expression in generators [Hard-coded max length 28]
} perm_t;

// Returns factorial of n. Should be used instead of tgamma(n+1).
int fac(int n){
    if(n == 1) return 1;
    return n * fac(n-1);
}

int powint(int x, int n){
    if(n == 1) return x;
    return x * pow(x, n - 1);
}

// Converts index x to its Lehmer code.
// Optimization: the last Lehmer digit is always 0. Maybe reasonable to replace by '\0'?
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

// Gets sum from Lehmer code, or equivalently length from Lehmer code
int CodeSum(int n, char code[]){
    int sum = 0;
    for (int i = 0; i < n; i++){
        sum += code[i];
    }
    return sum;
}

// Gets index from Lehmer code
int CodeToIndex(int n, char code[]){
    int sum = 0;
    for(int i = 0; i < n; i++){
        sum += code[i] * (int)tgamma(n - i);
    }
    return sum;
}

// Gets length of element from index
int IndexToLength(int n, int x){
    char code[8];
    return(CodeSum(n, LehmerCode(n, x, code)));
}

// Returns index of the composition x*y. Not y*x.
int MultiplyIndex(int n, int x, int y){
    char permx[8], permy[8], permp[8];
    IndexToPerm(n, x, permx);
    IndexToPerm(n, y, permy);
    for(int i = 0; i < n; i++){
        permp[i] = permy[(int)permx[i]];
    }
    return PermToIndex(n, permp);
}

// Returns index of permutation
int PermToIndex(int n, char* perm){
    char code[8];
    for(int i = 0; i < n; i++){
        int sum = 0;
        for(int j = i + 1; j < n; j++) if(perm[j] < perm[i]) sum++;
        code[i] = sum;
    }
    return CodeToIndex(n, code);
}

// Puts the permutation of index x into perm
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

// Gets one generator index [UNUSED?]
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

int BruhatSmaller(int n, int w1, int w2){
    char* exp1 = (char*)calloc(28, sizeof(char));
    char* exp2 = (char*)calloc(28, sizeof(char));
    ReducedExpression(n, w1, exp1);
    ReducedExpression(n, w2, exp2);

    //printf("Reduced expression for %d: ", w1); for(int i = 0; i < 28; i++) printf("%d ", exp1[i]); printf("\n");
    //printf("Reduced expression for %d: ", w2); for(int i = 0; i < 28; i++) printf("%d ", exp2[i]); printf("\n");

    int i = 0;
    for(int j = 0; j < 28; j++){
        //printf("exp1: %d\n", exp1[i]);
        if(exp1[j] == 0){
            free(exp1);
            free(exp2);
            //printf("Returned since exp1[%d] = 0\n", j);
            return 1;
        } 
        
        int found = 0;
        for(int k = i; k < 28; k++){
            if(exp2[k] == exp1[j]){
                //printf("Found %d at k = %d and j = %d!\n", exp1[j], k, j);
                found = 1;
                i = k + 1;
                break;
            }
        }

        if(found == 0){
            //printf("Returned due to not finding %d at j = %d\n", exp1[j], j);
            free(exp1);
            free(exp2);
            return 0;
        } 

    }
    free(exp1);
    free(exp2);
    //printf("Base case return!\n");
    return 1;
}

int TrueBruhatSmaller(int n, int w1, int w2){
    if(w1 == w2){
        return 0;
    }else{
        return BruhatSmaller(n, w1, w2);
    }
}