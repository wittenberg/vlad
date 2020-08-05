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

double hS(double s, double g0, double g1, double RQ);
double gX(double x, double g0, double g1, double RQ);
double s0(double w, double RA, double g0, double g1, double RQ);
double s1(double w, double RA, double g0, double g1, double RQ);
double FX(double s, double g0, double g1, double shape1, double shape2, double RQ);

double FWT2(double w, double RA, double g0, double g1, double shape1, double shape2, double RQ);

double Qi(int i, double w, double RA, double g0, double g1, double shape1, double shape2, double RQ);
double qij(int i, int j, double w, double RA, double g0, double g1, double shape1, double shape2, double RQ);
double racusum_beta_arl_mc(double h, double RA, double g0, double g1, double shape1, double shape2, int r = 600, int method = 1, double RQ = 1);
double racusum_beta_crit_mc(double L0, double RA, double g0, double g1, double shape1, double shape2, int method = 1, int r = 600, int jmax = 4, bool verbose = true, double RQ = 1);

#endif
