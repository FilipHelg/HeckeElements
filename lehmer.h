#include <stdio.h>
#include <stdlib.h>
#include "math.h"

int fac(int n);
int powint(int x, int n);
char* LehmerCode(int n, int x, char code[]);
int CodeSum(int n, char code[]);
int CodeToIndex(int n, char code[]);
int IndexToLength(int n, int x);
void IndexToPerm(int n, int x, char* perm);
int PermToIndex(int n, char* perm);
int MultiplyIndex(int n, int x, int y);
int FirstGeneratorIndex(int n, int x);
int ReducedExpression(int n, int x, char expression[]);
int BruhatSmaller(int n, int w1, int w2);
int TrueBruhatSmaller(int n, int w1, int w2);
int* ElementsBetween(int n, int x, int y);
int* ElementsBetween2(int n, int x, int y);