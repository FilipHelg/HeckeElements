#include <stdio.h>
#include "laurent.h"

// Returns a new Laurent polynomial with all coefficients set to zero
Laurent_t ZeroInitializeLaurent(){
    Laurent_t poly;
    for(int i = 0; i < 57; i++){
        poly.coeff[i] = 0;
    }
    return poly;
}

// Checks if a Laurent polynomial is non-zero. Returns 1/0 for true/false.
int HasNonZero(const Laurent_t *p){
    for(int j = 0; j < 57; j++){
        if(p->coeff[j] != 0) return 1;
    }
    return 0;
}

// Prints a Laurent polynomial in the terminal. Mostly used inline for DisplayHecke.
void DisplayLaurentPoly(Laurent_t poly){
    int first = 1;
    if(HasNonZero(&poly) == 0) printf("Zero polynomial!");
    for(int i = 0; i < 57; i++){
        if(poly.coeff[i] != 0){
            if(!first) printf(" + ");
            printf("%d*v^%d", poly.coeff[i], i - 28);
            first = 0;
        }
    }
}

// Returns the sum of two Laurent polynomials
Laurent_t SumLaurent(Laurent_t l1,  Laurent_t l2){
    Laurent_t sum = ZeroInitializeLaurent();
    for(int i = 0; i < 57; i++){
        sum.coeff[i] = l1.coeff[i] + l2.coeff[i];
    }
    return sum;
}

// Returns the product of two Laurent polynomials
Laurent_t MultiplyLaurent(Laurent_t l1, Laurent_t l2){
    Laurent_t product = ZeroInitializeLaurent();
    for(int i = 0; i < 57; i++){
        for(int j = 0; j < 57; j++){
            int index = i + j - 28;
            if(index < 57) product.coeff[index] += l1.coeff[i]*l2.coeff[j]; // Undefined behaviour when i + j - 28 > 56, but this should not happen as max length of word in S_8 is 28    
        }
    }
    return product;
}

Laurent_t Cutoff(Laurent_t l, int lastDeg){
    Laurent_t result = ZeroInitializeLaurent();
    for(int i = 28; i <= 28 + lastDeg; i++){
        result.coeff[i] = l.coeff[i];
    }
    return result;
}