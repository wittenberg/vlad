#include "racusum_beta_mc.h"

// helper functions
// double hS(double s, double g0, double g1) {return( 1 / ( 1 + exp(-g0-g1*s) ) );}
double hS(double s, double g0, double g1, double QS) {
  double pt = 1 / ( 1 + exp(-g0-g1*s) );
  return( (QS * pt) / (1 - pt + QS * pt) );
  }

double gX(double x, double g0, double g1, double QS) {return( -g0/g1 - log( QS * (1/x - 1) ) / g1 );}
double s0(double w, double QA, double g0, double g1, double QS) {return( gX( ( exp(-w) - 1 ) / ( QA - 1), g0, g1, QS) );}
double s1(double w, double QA, double g0, double g1, double QS) {return( gX( ( QA * exp(-w) - 1 ) / ( QA - 1), g0, g1, QS) );}
double FX(double s, double g0, double g1, double shape1, double shape2, double QS) {return(R::pbeta(gX(s, g0, g1, QS), shape1, shape2, true, false));}

double luFW2(double w, double QA, double g0, double g1, double shape1, double shape2, double QS, int lu) {
  double ires = 0, x0a,  x1a;
  auto f1 = [g0, g1, shape1, shape2, QS](double x) { return FX(exp(x), g0, g1, shape1, shape2, QS) * exp(x);};
  /* apply GL quadrature + integration by parts + change of variables */
  switch (lu) {
    case 1:{  /* left side CDF(w) (deterioration) */
      x0a = hS( s0(w, QA, g0, g1, 1), g0, g1, QS );      /* lower limit no trans. */
      x1a = hS( 1, g0, g1, QS );                         /* upper limit no trans. */
      ires = gauss_kronrod<double, 31>::integrate(f1, log(x0a), log(x1a), 5, 1e-9);
      ires += 1 - x1a - (1 - x0a) * FX(x0a, g0, g1, shape1, shape2, QS);
      break;
    }
    case 2:{  /* right side CDF(w) (deterioration) */
      x0a = hS( s1(w, QA, g0, g1, 1), g0, g1, QS ) ;      /* lower limit no trans. */
      x1a = hS( 1, g0, g1, QS ) ;                         /* upper limit no trans. */
      ires = gauss_kronrod<double, 31>::integrate(f1, log(x0a), log(x1a), 5, 1e-9);
      ires = x1a - x0a*FX(x0a, g0, g1, shape1, shape2, QS) - ires;
      break;
    }
    case 3:{  /* left side CDF(w) (improvement) */
      x0a = hS( 0, g0, g1, QS);                            /* lower limit no trans. */
      x1a = hS( s1(w, QA, g0, g1, 1), g0, g1, QS );        /* upper limit no trans. */
      ires = gauss_kronrod<double, 31>::integrate(f1, log(x0a), log(x1a), 5, 1e-9);
      ires = x1a*FX(x1a, g0, g1, shape1, shape2, QS) - ires;
      // ires = x1a*FX(x1a, g0, g1, shape1, shape2, QS) - x0a*FX(x0a, g0, g1, shape1, shape2, QS) - ires;
      break;
    }
    case 4:{  /* right side CDF(w) (improvement) */
      x0a = hS( 0, g0, g1, QS );                          /* lower limit no trans. */
      x1a = hS( s0(w, QA, g0, g1, 1), g0, g1, QS );       /* upper limit no trans. */
      ires = gauss_kronrod<double, 31>::integrate(f1, log(x0a), log(x1a), 5, 1e-9);
      ires += (1 - x1a) * FX(x1a, g0, g1, shape1, shape2, QS);
      // ires = (1 - x1a) * FX(x1a, g0, g1, shape1, shape2, QS) - (1 - x0a) * FX(x0a, g0, g1, shape1, shape2, QS) + ires;
    break;
    }
  }
  return(ires);
}

/* Comlete CDF(w) */
// [[Rcpp::export(.FWT2)]]
double FWT2(double w, double QA, double g0, double g1, double shape1, double shape2, double QS) {

  double w0, w1, w2, w3, res, logQA = log(QA);

  if (QA > 1) {   /* Decting deterioration QA > 1 */
    w0 = -log(1 + (QA-1) * hS(1, g0, g1, 1));           /* lower limit left side of FW */
    w1 = -log(1 + (QA-1) * hS(0, g0, g1, 1));           /* upper limit left side of FW */
    w2 = -log(1 + (QA-1) * hS(1, g0, g1, 1)) + logQA; /* lower limit right side of FW */
    w3 = -log(1 + (QA-1) * hS(0, g0, g1, 1)) + logQA; /* upper limit right side of FW */
    if (w0 > w) {res = 0;}
    else if ((w0 <= w) & (w <= w1)) {res = luFW2(w,  QA, g0, g1, shape1, shape2, QS, 1);}
    else if ((w1 <= w) & (w <= w2)) {res = luFW2(w1, QA, g0, g1, shape1, shape2, QS, 1);}
    else if ((w2 <= w) & (w <= w3)) {res = luFW2(w1, QA, g0, g1, shape1, shape2, QS, 1) + luFW2(w, QA, g0, g1, shape1, shape2, QS, 2);}
    else {res = 1;}
  } else if ((QA < 1) & (QA >0)){   /* Decting improvement QA > 0 and QA < 1 */
    w0 = -log(1 + (QA-1) * hS(0, g0, g1, 1)) + logQA; /* lower limit left side of FW */
    w1 = -log(1 + (QA-1) * hS(1, g0, g1, 1)) + logQA; /* upper limit left side of FW */
    w2 = -log(1 + (QA-1) * hS(0, g0, g1, 1));           /* lower limit right side of FW */
    w3 = -log(1 + (QA-1) * hS(1, g0, g1, 1));           /* upper limit right side of FW */
    if (w0 > w) {res = 0;}
    else if ((w0 <= w) & (w <= w1)) {res = luFW2(w,  QA, g0, g1, shape1, shape2, QS, 3);}
    else if ((w1 <= w) & (w <= w2)) {res = luFW2(w1, QA, g0, g1, shape1, shape2, QS, 3);}
    else if ((w2 <= w) & (w <= w3)) {res = luFW2(w1, QA, g0, g1, shape1, shape2, QS, 3) + luFW2(w, QA, g0, g1, shape1, shape2, QS, 4);}
    else {res = 1;}
  } else{
    res = 0;
  }
  return(res);
}

double Qi(int i, double w, double QA, double g0, double g1, double shape1, double shape2, double QS) {
  return( FWT2(-i*w + w/2, QA, g0, g1, shape1, shape2, QS) );
}

double qij(int i, int j, double w, double QA, double g0, double g1, double shape1, double shape2, double QS) {
  return( FWT2((j-i)*w + w/2, QA, g0, g1, shape1, shape2, QS) -
          FWT2((j-i)*w - w/2, QA, g0, g1, shape1, shape2, QS) );
}

// [[Rcpp::export(.racusum_beta_arl_mc)]]
double racusum_beta_arl_mc(double h, double QA, double g0, double g1, double shape1, double shape2, int r, int method, double QS) {
  double w, value = 0;

  switch (method) {
  case 1: {   /* ARL calculation: Toeplitz SPRT approach */
    double al, be, et, ga, de;
    int N = r, N1, i, j;
    NumericVector x(N), y(N), z(N), phi(N), psi(N), a(2*N - 1), b1(N, 1.), b2(N);

    N1 = N - 1;
    w = (2. * h)/(2. * r - 1.);

    for (i = 0; i < a.size(); i++) {
      a[i] = -( FWT2(-w*(i-N1) + w/2, QA, g0, g1, shape1, shape2, QS) -
                FWT2(-w*(i-N1) - w/2, QA, g0, g1, shape1, shape2, QS) );
    }

    a[N1] += 1;

    for (i = 0; i < N; i++ ) b2[i] = FWT2(-i*w - w/2, QA, g0, g1, shape1, shape2, QS);

    x[0] = 1 / a[N1];
    y[0] = 1 / a[N1];
    phi[0] = b1[0] / a[N1];
    psi[0] = b2[0] / a[N1];

    for (i = 1; i < N; i++) {
      al = 0;
      for (j = 0; j < i; j++) al += a[N1 + i - j] * x[j];
      ga = 0;
      for (j = 0; j < i; j++) ga += a[N1 - 1 - j] * y[j];
      et = -b1[i];
      for (j = 0; j < i; j++) et += a[N1 + i - j] * phi[j];
      de = -b2[i];
      for (j = 0; j < i; j++) de += a[N1 + i - j] * psi[j];

      be = 1 - al*ga;

      z[0] = ( -ga * x[0] ) / be;
      for (j = 1; j < i; j++) z[j] = ( y[j - 1] - ga * x[j] ) / be;
      z[i] = y[i - 1] / be;

      x[0] = x[0] / be;
      for (j = 1; j < i; j++) x[j] = ( x[j] - al * y[j - 1] ) / be;
      x[i] = ( -al*y[i - 1] ) / be;

      for (j = 0; j <= i; j++) y[j] = z[j];

      for (j = 0; j < i; j++) {
        phi[j] = ( phi[j] - et*z[j] );
        psi[j] = ( psi[j] - de*z[j] );
      }

      phi[i] = -et * z[i];
      psi[i] = -de * z[i];
    }
    be = phi[0] / ( 1 - psi[0] );
    value = be;
    break;
  }
  case 2: {   /* ARL calculation: Brook and Evans 1972 approach */
    int t = r, i, j;
    NumericVector arl;
    arma::mat rrr, id;
    arma::colvec b;

    rrr.zeros(t, t);  /* transition probability matrix */
    id.eye(t, t);     /* identity matrix */
    b.ones(t);        /* vector of ones  */

    w = (2. * h)/(2. * r - 1.);

    /* first column of transition probability matrix */
    for (i = 0; i < t; i++ ) rrr(i, 0) = Qi(i, w, QA, g0, g1, shape1, shape2, QS);
    /* first row (without first entry) of transition probability matrix */
    for (j = 1; j < t; j++ ) rrr(0, j) = qij(0, j, w, QA, g0, g1, shape1, shape2, QS);
    /* fill rest of transition probability matrix */
    for (i = 1; i < t; i++ ) {
      for (j = 1; j < t; j++ ) {
        rrr(i, j) = qij(i, j, w, QA, g0, g1, shape1, shape2, QS);
      }
    }

    /* solve system of linear equations */
    arl = arma::solve(id - rrr, b, arma::solve_opts::fast + arma::solve_opts::no_approx);
    value = arl[0];
    break;
  }
  }
  return(value);
}

// [[Rcpp::export(.racusum_beta_crit_mc)]]
double racusum_beta_crit_mc(double L0, double QA, double g0, double g1, double shape1, double shape2, int method, int r, int jmax, bool verbose, double QS) {
  double L1, h, h1;
  int i, j, dh;
  for (i = 1; i < 10; i++ ) {
    L1 = racusum_beta_arl_mc(double(i), QA, g0, g1, shape1, shape2, r, method, QS);
    if ( verbose ) Rcpp::Rcout << "h = " <<  i << "\t" << "ARL = " << L1 << std::endl;
    if ( L1 > L0 ) break;
  }
  h1 = i;

  for (j = 1; j <= jmax; j++ ) {
    for (dh = 1; dh <= 19; dh++ ) {
      h = h1 + pow(-1., j) * dh / pow(10., j);
      L1 = racusum_beta_arl_mc(h, QA, g0, g1, shape1, shape2, r, method, QS);
      if ( verbose ) Rcpp::Rcout << "h = " <<  h << "\t" << "ARL = " << L1 << std::endl;
      if ( ((j % 2 == 1) & (L1 < L0) ) | ((j % 2 == 0) & (L1 > L0)) ) break;
    }
    h1 = h;
  }
  if ( L1 < L0 ) h = h + 1 / pow(10., jmax);
  return h;
}
