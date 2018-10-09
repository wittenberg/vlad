#ifndef RACUSUM_SIM_H
#define RACUSUM_SIM_H

#include "vlad.h"
#include <math.h>
#include <algorithm>

using namespace Rcpp;

double gettherisk(int parsonnetscore, NumericVector coeff);

double loglikelihood(DataFrame df, NumericVector coeff, double R0 = 1, double RA = 2, bool yemp = true);
double loglikelihood2(DataFrame df, NumericVector coeff, NumericVector coeff2, double R0 = 1, double RA = 2, double RQ = 1);
double loglikelihood3(DataFrame df, double R0 = 1, double RA = 2);
int cusum_arl_sim(int r, double h, DataFrame df, double R0 = 1, double RA = 2);
int racusum_arl_sim(int r, NumericVector coeff, double h, DataFrame df, double R0 = 1, double RA = 2, bool yemp = true);
int racusum_arloc_sim(int r, NumericVector coeff, NumericVector coeff2, double h, DataFrame df, double R0 = 1, double RA = 2, double RQ = 1);
int racusum_adoc_sim(int r, NumericVector coeff, NumericVector coeff2, double h, DataFrame df, double R0 = 1, double RA = 2,  double RQ = 1, int m = 5, int type = 1);

#endif
