#include "kl_structure_bulk_opt_parallel.c"

int main(int argc, int* argv){
    FastContext *ctx;
    SparseHecke *left, *right, *prod, *temp1, *temp2;
    DenseAccum *accum;
    BenchStats *bench;
    int xid;



    MultiplyHeckeSparse(ctx, xid, left, right, prod, temp1, temp2, accum, bench, 0);

}