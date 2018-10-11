#ifndef EOCUSUM_SIM_H
#define EOCUSUM_SIM_H

#include "vlad.h"
#include <math.h>
#include <algorithm>

using namespace Rcpp;

double gettherisk(int parsonnetscore, NumericVector coeff);

//double optimal_k(double QA, IntegerVector parsonnetscores, NumericVector coeff);
double optimal_k(double QA, DataFrame df, NumericVector coeff, bool yemp = true);
// Rcpp::Nullable<Rcpp::NumericVector> coeff = R_NilValue

double calceo(DataFrame df, NumericVector coeff, bool yemp = true);

double calceo2(DataFrame df, NumericVector coeff, NumericVector coeff2, double QS = 1);

int eocusum_arl_sim(int r, double k, double h, DataFrame df, NumericVector coeff, bool yemp = true, int side = 1);

int eocusum_arloc_sim(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS = 1, int side = 1);

int eocusum_ad_sim(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS = 1, int side = 1, int type= 1, int m = 5);
int eocusum_ad_sim11(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m);
int eocusum_ad_sim12(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m);
int eocusum_ad_sim21(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m);
int eocusum_ad_sim22(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m);

#endif
