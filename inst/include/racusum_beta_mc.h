#ifndef RACUSUM_BETA_MC_H
#define RACUSUM_BETA_MC_H

#include "vlad.h"
#include <math.h>
#include <algorithm>
#include <stdlib.h> // for NULL
#include <stdio.h>
#include <R.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using namespace Rcpp;
using namespace std;
using namespace boost::math::quadrature;

double hS(double s, double g0, double g1, double QS);
double gX(double x, double g0, double g1, double QS);
double s0(double w, double QA, double g0, double g1, double QS);
double s1(double w, double QA, double g0, double g1, double QS);
double FX(double s, double g0, double g1, double shape1, double shape2, double QS);

double FWT2(double w, double QA, double g0, double g1, double shape1, double shape2, double QS);

double Qi(int i, double w, double QA, double g0, double g1, double shape1, double shape2, double QS);
double qij(int i, int j, double w, double QA, double g0, double g1, double shape1, double shape2, double QS);
double racusum_beta_arl_mc(double h, double QA, double g0, double g1, double shape1, double shape2, int r = 600, int method = 1, double QS = 1);
double racusum_beta_crit_mc(double L0, double QA, double g0, double g1, double shape1, double shape2, int method = 1, int r = 600, int jmax = 4, bool verbose = true, double QS = 1);

#endif
