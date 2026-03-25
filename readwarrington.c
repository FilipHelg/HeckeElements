#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv){
    FILE *fptr;
    
    fptr = fopen("S8.txt", "r");

    // Logic for getting # of characters in .txt file
    size_t pos = ftell(fptr);
    fseek(fptr, 0, SEEK_END);
    size_t length = ftell(fptr);
    fseek(fptr, pos, SEEK_SET);
    printf("Length of .txt file = %d\n", length);

    // Read file into string
    char* string = (char*)calloc(length, sizeof(char));
    fread(string, sizeof(char), length, fptr);

    // Count elements in file
    int elements = 1;
    for(int i = 0; i < length; i++) if(string[i] == '\n') elements++;
    printf("# of non-trivial pairs = %d\n", elements);
    
    char* firstperm = (char*)calloc(length*10, sizeof(char));
    char* secondperm = (char*)calloc(length*10, sizeof(char));
    char* coeffdata = (char*)calloc(length*20, sizeof(char));

    int i = 0, count = 0;
    while(i < length){
        int m = 0;
        while(string[i] != '\n' && i < length){
            char* data;
            if(m%3==0)data = firstperm + count*10;
            else if(m%3==1)data = secondperm + count*10;
            else {data = coeffdata + count*20; count++;}
            m++;
            int k = 0;
            while(string[i] != ' ' && i < length){
                data[k] = string[i];
                i++;k++;
            }
            data[k] = '\0';
            //printf(data);printf("\n");
            if(i < length) i++;
        }
        if(i < length) i++;
    }
    
    int whati = -1, foundcomma = 0;
    int max = 0;
    for(int i = 0; i < elements; i++){
        char* data = coeffdata + i*20;
        int coeffcount = 1;
        foundcomma = 0;
        for(int j = 0; j < 20; j++){
            //printf("%c\n", data[j]);
            if(data[j] == '\0') break;
            else if(data[j] == ',') {foundcomma = 1; coeffcount++;}
            if(max < coeffcount) {max = coeffcount; whati = i;}
        }
        if(foundcomma == 0) printf("No comma at i = %d\n", i);
    }
    //max++; max = max - max/2;
    printf("Maximal # of coefficients in a KL-polynomial = %d at i = %d\n", max, whati);

    fclose(fptr);
    free(string); free(firstperm); free(secondperm); free(coeffdata);
}