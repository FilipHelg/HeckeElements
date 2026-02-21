typedef struct
{
    int coeff[57]; //28 terms v^n, 1 constant term, 28 terms v^-n
}Laurent_t;

Laurent_t ZeroInitializeLaurent();
int HasNonZero(const Laurent_t *p);
void DisplayLaurentPoly(Laurent_t poly);
Laurent_t SumLaurent(Laurent_t l1,  Laurent_t l2);
Laurent_t MultiplyLaurent(Laurent_t l1, Laurent_t l2);
Laurent_t Cutoff(Laurent_t l, int lastDeg);