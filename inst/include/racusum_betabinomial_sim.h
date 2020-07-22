#ifndef RACUSUM_BETA_SIM_H
#define RACUSUM_BETA_SIM_H

#include "vlad.h"
#include <math.h>
#include <algorithm>

using namespace Rcpp;

int racusum_betabinoial_arl_sim(int r, double shape1, double shape2, NumericVector coeff, double h, double RA = 2, int rs = 71, double RQ = 1);

#endif
