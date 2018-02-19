#ifndef EOCUSUM_SIM_H
#define EOCUSUM_SIM_H

#include "vlad.h"
#include "misc.h"
#include <math.h>
#include <algorithm>

using namespace Rcpp;

int eocusum_adoc_sim(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS = 1, int side = 1, int type= 1, int m = 5);
int eocusum_adoc_sim11(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m);
int eocusum_adoc_sim12(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m);
int eocusum_adoc_sim21(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m);
int eocusum_adoc_sim22(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m);

#endif
