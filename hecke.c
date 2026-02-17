#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lehmer.h"
#include "laurent.h"

// Zero initializes all coefficients of a Hecke element [UNUSED?]
void SetToZero(int n, Laurent_t H[]){
    for (int i = 0; i < tgamma(n + 1); i++){
        H[i] = ZeroInitializeLaurent();
    }
}

// Prints out a Hecke element in the terminal
void DisplayHecke(int count, Laurent_t* element){
    int first = 1;
    for(int i = 0; i < count; i++){
        if(HasNonZero(&element[i])){
            if(!first) printf(" + ");
            printf("(");
            DisplayLaurentPoly(element[i]);
            printf(")H[%d]", i);
            first = 0;
        }
    }
    printf("\n");
}

// Adds two Hecke elements and returns the sum
Laurent_t* AddHecke(int n, Laurent_t H1[], Laurent_t H2[]){
    Laurent_t* sum = (Laurent_t*)calloc((size_t)tgamma(n + 1), sizeof(Laurent_t));
    for (int i = 0; i < tgamma(n + 1); i++){
        sum[i] = SumLaurent(H1[i], H2[i]);
    }
    return sum;
}

// Multiply with simple reflection on the left
Laurent_t* MultiplySimpleHecke(int n, int s, Laurent_t factor, Laurent_t H[]){
    Laurent_t* product = (Laurent_t*)calloc((size_t)tgamma(n + 1), sizeof(Laurent_t));
    if (s == 0) {for(int i = 0; i < tgamma(n+1); i++) product[i] = MultiplyLaurent(factor, H[i]); return product;}
    
    // Create a factor (v^-1 - v^1) ready to use (loop invariant)
    Laurent_t lowerTermFactor = ZeroInitializeLaurent();
    lowerTermFactor.coeff[27] = 1;
    lowerTermFactor.coeff[29] = -1;

    for(int i = 0; i < tgamma(n + 1); i++){
        int composition = MultiplyIndex(n, s, i);
        if(IndexToLength(n, composition) == IndexToLength(n, i) + 1){
            product[composition] = SumLaurent(MultiplyLaurent(H[i], factor), product[composition]);  // H_sw
        }else{
            product[composition] = SumLaurent(MultiplyLaurent(H[i], factor), product[composition]);  // H_sw
            product[i] = SumLaurent(MultiplyLaurent(lowerTermFactor, MultiplyLaurent(H[i], factor)), product[i]);    // (v^-1 - v^1)H_w
        }
    }
    return product;
}

// Multiplies two Hecke elements. Specifically calculates H1*H2, not H2*H1.
Laurent_t* MultiplyHecke(int n, Laurent_t H1[], Laurent_t H2[]){
    Laurent_t* product = (Laurent_t*)calloc((size_t)tgamma(n + 1), sizeof(Laurent_t));
    Laurent_t* term = (Laurent_t*)calloc((size_t)tgamma(n + 1), sizeof(Laurent_t));
    Laurent_t id = ZeroInitializeLaurent();
    id.coeff[28] = 1;
    for(int i = 0; i < tgamma(n+1); i++){
        if(HasNonZero(&H1[i]) == 0){/*printf("H[%d] has coefficient zero\n", i);*/ continue;} //Skip when coefficient is zero

        char* generators = (char*)calloc(28, sizeof(char));
        int expression[28];
        ReducedExpression(n, i, generators);
        for(int l = 0; l < 28; l++) if(generators[l] != 0) expression[l] = fac(n - generators[l]); else expression[l] = 0; // Convert expression from generator indexation to Lehmer indexation
        
        // Debug
        /*printf("Expression for %d is ", i);
        for(int j = 0; j < 28; j++) if (expression[j] != 0) printf("%d", expression[j]); printf("\n");*/
        
        term = H2; // start out the next term of the product as all of H2
        for(int k = 27; k >= 0; k--){ // Multiply successively with generators
            int generator = expression[k];

            //Debug
            //printf("Current generator is %d\n", generator);

            if(generator == 0) continue; //Skip when generator is identity

            //Previous solution, works!
            /*Laurent_t* new_term = (Laurent_t*)calloc((size_t)tgamma(n + 1), sizeof(Laurent_t)); //[UNNECESSARY MEMORY USAGE? NEVER FREED?]
            new_term = MultiplySimpleHecke(n, generator, id, term);
            term = new_term;*/

            //New solution, also seems to work! Less memory allocation.
            term = MultiplySimpleHecke(n, generator, id, term); //[REASONABLE CHANGE? OR PREVIOUS CORRECT FOR SOME REASON?]

            //Debug
            /*printf("Midterm at k = %d: ", k);
            //DisplayHecke(fac(n),term);*/
        }
        //Debug
        /*printf("Term from multiplication with %d is ", i);
        DisplayHecke(fac(n),term);*/

        //Add term to total product
        product = AddHecke(n, MultiplySimpleHecke(n, 0, H1[i], term), product);
    }
    //free(term);
    return product;
}

/*Laurent_t* MultiplyHecke(int n, Laurent_t H1[], Laurent_t H2[]){
    int count = fac(n);
    Laurent_t* product = (Laurent_t*)calloc(count, sizeof(Laurent_t));

    for(int i = 0; i < count; i++){
        char generators[28];
        int expression[28];
        ReducedExpression(n, i, generators);
        for(int l = 0; l < 28; l++) expression[l] = fac(n - generators[l]);
    }
}*/

void SetKL(int n, int x, Laurent_t KLelement[]){
    char l = IndexToLength(n, x);
    if(l == 0){
        KLelement[0].coeff[28] = 1; // H_e
    }else if(l==1){
        KLelement[0].coeff[29] = 1; // vH_e
        KLelement[x].coeff[28] = 1; // H_n
    }else{
        //int s = FirstGeneratorIndex(n, x);
    }
}

Laurent_t* MultiplySimpleHecke2(int n, int s, Laurent_t H[]){
    int count = (int)tgamma(n + 1); // n! elements
    Laurent_t* product = calloc(count, sizeof(Laurent_t));

    // Factor for lower term in Hecke algebra: v^-1 - v^1
    Laurent_t lowerTermFactor = ZeroInitializeLaurent();
    lowerTermFactor.coeff[27] = 1;
    lowerTermFactor.coeff[29] = -1;

    for(int i = 0; i < count; i++){
        if(!HasNonZero(&H[i])) continue;

        int composition = MultiplyIndex(n, s, i);

        if(IndexToLength(n, composition) == IndexToLength(n, i) + 1){
            // Length increased, simple multiplication
            product[composition] = SumLaurent(H[i], product[composition]);
        } else {
            // Length did not increase, add lower term
            product[composition] = SumLaurent(H[i], product[composition]);
            product[i] = SumLaurent(MultiplyLaurent(lowerTermFactor, H[i]), product[i]);
        }
    }

    return product;
}

Laurent_t* MultiplyHecke2(int n, Laurent_t H1[], Laurent_t H2[]){
    int count = (int)tgamma(n + 1);
    Laurent_t* product = calloc(count, sizeof(Laurent_t));

    for(int i = 0; i < count; i++){
        if(!HasNonZero(&H2[i])) continue;

        // Start term correctly
        Laurent_t* term = calloc(count, sizeof(Laurent_t));
        term[i] = H2[i];

        // Find reduced expression in generators for i
        char expression[28];
        int len = ReducedExpression(n, i, expression);
        for(int l = 0; l < 28; l++) expression[l] = (int)tgamma(n - expression[l] + 1); // Convert expression from generator indexation to Lehmer indexation

        // Multiply succesively with each generator
        for(int k = 0; k < len; k++){
            int s = expression[k];
            term = MultiplySimpleHecke2(n, s, term);
        }

        // Add term to product
        Laurent_t* new_product = AddHecke(n, product, term);
        free(product);
        free(term);
        product = new_product;
    }

    return product;
}

#define TEST_INDEXTOPERM 0
#define TEST_REDUCEDEXPRESSION 0
#define TEST_MULTIPLYINDEX 0
#define TEST_LAURENT 0
#define TEST_DISPLAYHECKE 0
#define TEST_MULTSIMPLEHECKE 0
#define TEST_MULTHECKE 1
#define TEST_ALTERNATE 0

int main(){
    /// TESTS ///

    if(TEST_INDEXTOPERM){
        int n = 3;
        char perm[8];
        for(int i = 0; i < (int)tgamma(n + 1); i++) {IndexToPerm(n, i, perm); printf("Index %d has perm ", i); for (int j = 0; j < n; j++) printf("%d", perm[j]); printf("\n");}
    }

    if(TEST_REDUCEDEXPRESSION){
        char* generators = (char*)calloc(28, sizeof(int));
        int* expression = (int*)calloc(28, sizeof(int));
        int n = 3;
        int x = 3;
        ReducedExpression(n, x, generators);
        printf("The reduced expression for Lehmer index %d is ", x);
        for(int i = 0; i < 28; i++) if(generators[i] != 0) printf("%d", generators[i]);
        printf("\nNote that the first element in the array is the last to be multiplied with when multiplying from the left\n");

        for(int l = 0; l < 28; l++) if(generators[l] != 0) expression[l] = fac(n - generators[l]); else expression[l] = 0;
        for(int i = 0; i < 28; i++) if(expression != 0) printf("%d ", expression[i]);
    }

    if(TEST_MULTIPLYINDEX){
        int n, x, y, z;
        n = 3;
        x = 0;
        y = 2;
        z = MultiplyIndex(n, x, y);

        printf("%d * %d = %d in S_%d\n", x, y, z, n);
    }


    /// Test of Laurent_t ///

    if(TEST_LAURENT){
        Laurent_t l1 = ZeroInitializeLaurent(), l2 = ZeroInitializeLaurent();

        l1.coeff[29] = 1;
        l1.coeff[5] = 2;

        l2.coeff[30] = 1;
        l2.coeff[36] = 5;

        DisplayLaurentPoly(l1); printf("\n");
        DisplayLaurentPoly(l2); printf("\n");
        DisplayLaurentPoly(MultiplyLaurent(l1, l2)); printf("\n");
        DisplayLaurentPoly(SumLaurent(l1, l2)); printf("\n");
    }


    /// Test of DisplayHecke ///

    if(TEST_DISPLAYHECKE){
        int n = 8;
        int count = (int)tgamma(n+1);
        Laurent_t *KLelement = (Laurent_t*)calloc((size_t)count, sizeof(Laurent_t));
        Laurent_t factor = ZeroInitializeLaurent();
        factor.coeff[28] = 1;
        SetKL(n, 1, KLelement);
        DisplayHecke(count, KLelement);

        char perm[8];
        IndexToPerm(8, 24, perm);
        for(int i = 0; i < 8; i++) printf("%d ", perm[i]);
        printf("\n");

        int y = MultiplyIndex(n, 25489, 2442);
        printf("%d", y);

        DisplayHecke(n, KLelement);
        DisplayHecke(n, MultiplySimpleHecke(8, 2, factor, KLelement));
    }


    /// Test of MultiplySimpleHecke ///

    if(TEST_MULTSIMPLEHECKE){
        int n = 8;
        int count = (int)tgamma(n+1);
        Laurent_t *h1 = (Laurent_t*)calloc((size_t)count, sizeof(Laurent_t)), *h2 = (Laurent_t*)calloc((size_t)count, sizeof(Laurent_t));
        Laurent_t factor = ZeroInitializeLaurent();
        factor.coeff[28] = 1;
        factor.coeff[30] = 2;

        h1[0].coeff[29] = 1;
        h1[1].coeff[28] = 1;

        h2[2].coeff[29] = 1;
        h2[5].coeff[27] = -2;

        DisplayHecke(count, h1);
        DisplayHecke(count, h2);
        DisplayHecke(count, MultiplySimpleHecke(n, 1, factor, h1));
        DisplayHecke(count, MultiplySimpleHecke(n, 5, factor, h2));
    }


    /// Test of MultiplyHecke ///

    if(TEST_MULTHECKE){
        int n = 8;
        int count = (int)tgamma(n+1);
        Laurent_t *h1 = (Laurent_t*)calloc((size_t)count, sizeof(Laurent_t)), *h2 = (Laurent_t*)calloc((size_t)count, sizeof(Laurent_t));

        //Initialize 1
        h1[0].coeff[29] = 1;
        h1[1].coeff[28] = 1;

        h2[2].coeff[29] = 1;
        h2[5].coeff[27] = -2;

        //Calculations and prints 1
        printf("\nTEST MULTIPLICATION 1\n");
        printf("Left factor: ");
        DisplayHecke(count, h1);
        printf("Right factor: ");
        DisplayHecke(count, h2);
        printf("Product: ");
        DisplayHecke(count, MultiplyHecke(n, h1, h2));

        //Reset 1
        h1[0] = ZeroInitializeLaurent();
        h1[1] = ZeroInitializeLaurent();
        h2[2] = ZeroInitializeLaurent();
        h2[5] = ZeroInitializeLaurent();

        //Initialize 2
        h1[0].coeff[29] = 1;
        h1[1].coeff[28] = 1;

        h2[0].coeff[29] = 1;
        h2[2].coeff[28] = 1;

        //Calculations and prints 2
        printf("\nTEST MULTIPLICATION 2\n");
        printf("Left factor: ");
        DisplayHecke(count, h1);
        printf("Right factor: ");
        DisplayHecke(count, h2);
        printf("Product: ");
        DisplayHecke(count, MultiplyHecke(n, h1, h2));

        //Reset 2 
        h1[0] = ZeroInitializeLaurent();
        h1[1] = ZeroInitializeLaurent();
        h2[0] = ZeroInitializeLaurent();
        h2[2] = ZeroInitializeLaurent();

        //Initialize 3 
        h1[0].coeff[29] = 1;
        h1[2].coeff[28] = 1;

        h2[0].coeff[29] = 1;
        h2[1].coeff[28] = 1;

        //Calculations and prints 3
        printf("\nTEST MULTIPLICATION 3\n");
        printf("Left factor: ");
        DisplayHecke(count, h1);
        printf("Right factor: ");
        DisplayHecke(count, h2);
        printf("Product: ");
        DisplayHecke(count, MultiplyHecke(n, h1, h2));
    }

    if(TEST_ALTERNATE){
        int n = 3;
        int count = (int)tgamma(n + 1);
        Laurent_t* h1 = calloc(count, sizeof(Laurent_t));
        Laurent_t* h2 = calloc(count, sizeof(Laurent_t));

        h1[0].coeff[28] = 1; // H_e
        h2[2].coeff[29] = 1; // H_w

        Laurent_t* prod = MultiplyHecke2(n, h1, h2);
        DisplayHecke(count, h1);
        DisplayHecke(count, h2);
        DisplayHecke(count, prod);
        free(h1);
        free(h2);
        free(prod);
    }

    printf("\n Remember: Measure-Command { ./hecke } to test run time!\n");
}