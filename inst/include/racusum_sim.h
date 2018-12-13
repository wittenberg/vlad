#ifndef RACUSUM_SIM_H
#define RACUSUM_SIM_H

#include "vlad.h"
#include <math.h>
#include <algorithm>

using namespace Rcpp;

double gettherisk(int parsonnetscore, NumericVector coeff) {
double logitp, risk;
logitp = coeff[0] + parsonnetscore * coeff[1];
risk = exp(logitp)/(1 + exp(logitp));
return risk;
}
double loglikelihood(DataFrame df, NumericVector coeff, double R0 = 1, double RA = 2, bool yemp = true);
int bcusum_arl_sim(int r, double h, DataFrame df, double R0 = 1, double RA = 2);

int racusum_arl_sim(int r, DataFrame pmix, double h, double RA, double RQ = 1, bool yemp = false);
int racusum_ad_sim(int r, DataFrame pmix, double h, double RA, double RQ = 1, int m = 50, int type = 1);
#endif
