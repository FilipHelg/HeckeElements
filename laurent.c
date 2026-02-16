#include <stdio.h>
#include "laurent.h"

Laurent_t ZeroInitializeLaurent(){
    Laurent_t poly;
    for(int i = 0; i < 57; i++){
        poly.coeff[i] = 0;
    }
    return poly;
}

int HasNonZero(const Laurent_t *p){
    for(int j = 0; j < 57; j++){
        if(p->coeff[j] != 0) return 1;
    }
    return 0;
}

void DisplayLaurentPoly(Laurent_t poly){
    int first = 1;
    for(int i = 0; i < 57; i++){
        if(poly.coeff[i] != 0){
            if(!first) printf(" + ");
            printf("%d*v^%d", poly.coeff[i], i - 28);
            first = 0;
        }
    }
}

Laurent_t SumLaurent(Laurent_t l1,  Laurent_t l2){
    Laurent_t sum = ZeroInitializeLaurent();
    for(int i = 0; i < 57; i++){
        sum.coeff[i] = l1.coeff[i] + l2.coeff[i];
    }
    return sum;
}

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