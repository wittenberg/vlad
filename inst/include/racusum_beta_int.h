#ifndef RACUSUM_BETA_INT_H
#define RACUSUM_BETA_INT_H

#include "vlad.h"
#include "misc.h"
#include <math.h>
#include <algorithm>
#include <stdlib.h> // for NULL
#include <stdio.h>
#include <R.h>

#include <boost/math/quadrature/gauss_kronrod.hpp>

using namespace Rcpp;
using namespace std;
using namespace boost::math::quadrature;

double Tn(double z, int j);
double FWT2(double w, double RA, double g0, double g1, double shape1, double shape2, double RQ);
double hS(double s, double g0, double g1, double RQ);
double s0(double w, double RA, double g0, double g1, double RQ);
double s1(double w, double RA, double g0, double g1, double RQ);

double gX(double x, double g0, double g1, double RQ);
double gXp(double x, double g1) {return(  1/(g1*x*(1-x)) );}
double s0p(double w, double RA, double g0, double g1, double RQ) {return( gXp(( exp(-w) - 1 ) / ( RA - 1), g1) * -(exp(-w)/(RA-1)) );}
double s1p(double w, double RA, double g0, double g1, double RQ) {return( gXp(( RA * exp(-w) - 1 ) / ( RA - 1), g1) * -(RA*exp(-w)/(RA-1)) );}

double f2(double w, double RA, double RQ, double g0, double g1, double shape1, double shape2);
double integ_t62(double xl, double xu, int j, double loc, double RA, double RQ, double g0, double g1, double shape1, double shape2);
double racusum_beta_arl_int(double h, int N, double RA, double RQ, double g0, double g1, double shape1, double shape2, bool pw);

#endif
