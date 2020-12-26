#ifndef RACUSUM_MC_H
#define RACUSUM_MC_H

#include "vlad.h"
#include "misc.h"
#include <math.h>
#include <algorithm>

using namespace Rcpp;

arma::vec subset_vec(const arma::vec& x) {
  return x.elem(find(x > 0));
}

arma::vec so_subset_params(arma::vec x, arma::uvec pos, arma::vec vals) {
  // Subset and set equal to
  x.elem(pos) = vals;
  return x;
}

double racusum_arl_mc(NumericMatrix pmix, double RA, double RQ, double h, double scaling = 600, int rounding = 1, int method = 1);
double racusum_crit_mc(NumericMatrix pmix, double L0, double RA, double RQ, double scaling = 600, int rounding = 1, int method = 1, int jmax = 4, bool verbose = false);

#endif
