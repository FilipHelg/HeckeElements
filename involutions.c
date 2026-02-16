#include "stdio.h"
#include "stdlib.h"

void DisplayInvolution(char* involution, int n){
    for(int i = 0; i < n; i++){
        printf("%c", involution[i]);
    }
    printf("\n");
}

int main(){
    int n = 0;
    printf("Enter value of n for symetirc group S_n: ");
    scanf(" %d", &n);
    char *inv = (char*)malloc(n*sizeof(char));
    printf("Enter involution in full one-line format: ");

    char currentChar = 0;
    int counter = 0;
    scanf(" %c", &currentChar);
    while(currentChar != '\n'){
        inv[counter] = currentChar;
        counter++;
        scanf("%c", &currentChar);
    }
    DisplayInvolution(inv, n);
    free(inv);
}

