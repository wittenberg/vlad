#ifndef RACUSUM_DISCRETEBETA_SIM_H
#define RACUSUM_DISCRETEBETA_SIM_H

#include "vlad.h"
#include <math.h>
#include <algorithm>

using namespace Rcpp;

int racusum_discretebeta_arl_sim(int r, double shape1, double shape2, NumericVector coeff, double h, double RA, int rs, double RQ);

#endif
