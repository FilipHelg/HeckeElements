/* Build the required KL context implementation in this TU, then include the
    nonzero triple checker implementation under a renamed main so we can call
    it through a compatibility wrapper. */

/* Bring in core context implementation (rename its main to avoid collision). */
#define main kl_structure_bulk_opt_parallel_main
#include "kl_structure_bulk_opt_parallel.c"
#undef main

/* When including the nonzero checker, skip its internal re-include of the
    context (we already compiled it above), and rename its main so we can call
    it directly. */
#define KL_TRIPLE_CHECK_NONZERO_SKIP_INCLUDE
#define main kl_triple_check_nonzero_included_main
#include "kl_triple_check_nonzero.c"
#undef main
#undef KL_TRIPLE_CHECK_NONZERO_SKIP_INCLUDE

#include <stdio.h>
#include <stdlib.h>

static void Usage(const char *prog) {
    printf("Usage: %s n [--dual-bin PATH] [--kl-data PATH] [--out PATH]\n", prog);
    printf("This wrapper forwards directly to kl_triple_check_nonzero.\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) { Usage(argv[0]); return 1; }

    int n = atoi(argv[1]);
    if (n < 3 || n > 8) { fprintf(stderr, "n must be between 3 and 8.\n"); return 1; }

    return kl_triple_check_nonzero_included_main(argc, argv);
}
