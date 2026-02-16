#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int CalculateNumInvolutions(int n){
    int sum = 1; // Identiteten är alltid en involution
    for(int i = 1; i <= n/2; i++){
        // addera antal involutioner med i transpositioner
        sum += (tgamma(n+1)/(tgamma(2*i+1)*tgamma(n-2*i+1)))*(tgamma(2*i+1)/(tgamma(i+1)*pow(2, i)));
    }
    return sum;
}

int main(){
    int n = 0;
    printf("Enter value of n for symmetric group S_n: ");
    scanf("%d", &n);

    int numInvs = CalculateNumInvolutions(n);
    printf("%d", numInvs);
}