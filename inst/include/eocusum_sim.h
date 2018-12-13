#ifndef EOCUSUM_SIM_H
#define EOCUSUM_SIM_H

#include "vlad.h"
#include <math.h>
#include <algorithm>

using namespace Rcpp;

double optimal_k(DataFrame pmix, double RA, bool yemp = FALSE);

int eocusum_arl_sim(int r, DataFrame pmix, double k, double h, double RQ = 1, bool yemp = FALSE, int side = 1);

int eocusum_ad_sim(int r, DataFrame pmix, double k, double h, double RQ, int side, int type, int m);
int eocusum_ad_sim11(int r, DataFrame pmix, double k, double h, double RQ, int m);
int eocusum_ad_sim12(int r, DataFrame pmix, double k, double h, double RQ, int m);
int eocusum_ad_sim21(int r, DataFrame pmix, double k, double h, double RQ, int m);
int eocusum_ad_sim22(int r, DataFrame pmix, double k, double h, double RQ, int m);

#endif
