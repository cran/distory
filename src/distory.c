#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP phycpp_bin_trees(SEXP treelist);
SEXP phycpp_compute_tree_distance_set(SEXP trees, SEXP verbose);
SEXP gromov_distmatrix(SEXP distmatrix, SEXP bDeltas, SEXP scale_method);

static R_CallMethodDef Call_entries[] = {
    {"phycpp_bin_trees", (DL_FUNC) &phycpp_bin_trees, 1},
    {"phycpp_compute_tree_distance_set", (DL_FUNC) &phycpp_compute_tree_distance_set, 2},
    {"gromov_distmatrix", (DL_FUNC) &gromov_distmatrix, 3},
    {NULL, NULL, 0}
};

void R_init_distory(DllInfo *info)
{
    R_registerRoutines(info, NULL, Call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
