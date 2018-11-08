#ifndef RACUSUM_MC_H
#define RACUSUM_MC_H

#include "vlad.h"
#include "misc.h"
#include <math.h>
#include <algorithm>

using namespace Rcpp;

double gettherisk(int parsonnetscore, NumericVector coeff);
double racusum_arl_mc(NumericMatrix pmix, double RA, double RQ, double h, double scaling = 600, bool rounding = true, int method = 1);
double racusum_crit_mc(NumericMatrix pmix, double L0, double RA, double RQ, double scaling = 600, bool rounding = true, int method = 1, bool verbose = false);

#endif
