#include "racusum_mc.h"

// [[Rcpp::export(.racusum_arl_mc)]]
double racusum_arl_mc(NumericMatrix pmix, double RA, double RQ, double h, double scaling, int rounding, int method){

  NumericVector mp, mr1, mr2, z1, z2, pf, p1, p2;
  double value = 0;

  mp  = pmix(_, 0);                                       // Parsonnet distribution
  mr1 = pmix(_, 1);                                       // model
  mr2 = pmix(_, 2);                                       // either model or model-free

  pf = RQ * mr2 / ( 1 - mr2 + RQ * mr2 );                 // actual probability to die

  if (RA > 1) {                                           // Deterioration
    z1 = rev( log( 1  / ( 1 + (RA-1) * mr1 ) ) );         // survival
    z2 = rev( log( RA / ( 1 + (RA-1) * mr1 ) ) );         // death
    p1 = rev( (1-pf) * mp );
    p2 = rev( pf * mp );
  }
  else if ( (RA > 0) & (RA < 1) ) {                       // Improvement
    z1 = log( RA / ( 1 + (RA-1) * mr1 ) );                // death
    z2 = log( 1  / ( 1 + (RA-1) * mr1 ) );                // survival
    p1 = pf * mp;
    p2 = (1-pf) * mp;
  }
  else {
    Rcpp::Rcout << "Select a positive value for RA (RA > 0)." << "\n" << std::endl;
  }
  // paired rounding implementation
  if (rounding == 1) {
    arma::vec sz, iza, izb, pp, pa, pb, zz, ppp, zzz;

    sz = join_cols(as<arma::vec>(z1), as<arma::vec>(z2)) * scaling;
    iza = floor(sz);
    izb = ceil(sz);
    pp = join_cols(as<arma::vec>(p1), as<arma::vec>(p2));
    pa = pp % ( izb - sz );
    pb = pp % ( sz - iza );
    zz = join_cols(iza, izb);
    pp = join_cols(pa, pb);

    ppp = tapply3arma(pp, arma::conv_to<arma::ivec>::from(zz));
    zzz = unique(zz);

    switch (method) {
    // ARL calculation: Toeplitz SPRT approach
    case 1: {
      double al, be, et, ga, de;
      int i, j, N, N1;
      N = round(h * scaling);
      NumericVector x(N), y(N), z(N), phi(N), psi(N), a(2*N - 1), b1(N, 1.), b2(N);
      arma::vec jj;
      arma::uvec j0;

      N1 = N - 1;

      for (i = 0; i < zzz.n_elem; i++) a[N1 - zzz[i]] = -ppp[i];
      a[N1] = 1 + a[N1];

      for (i = 0; i < N1; i++ ) {
        jj = i + zzz;
        j0 = find(jj < 0);
        b2[i] = sum(ppp(j0));
      }

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
        for (j = 1; j < i; j++) z[j] = ( y[j-1] - ga * x[j] ) / be;
        z[i] = y[i-1] / be;

        x[0] = x[0] / be;
        for (j = 1; j < i; j++) x[j] = ( x[j] - al * y[j-1] ) / be;
        x[i] = ( -al*y[i-1] ) / be;

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
      // ARL calculation: Toeplitz full matrix inversion
    case 2: {
      int d, i, j, m;
      double tau, gamma, phi, f1, f2, f3, arl;
      arma::colvec b1, avg;
      arma::vec r, b, jj, x, y, z;
      arma::mat N;
      arma::uvec j0;

      d = round(h * scaling);
      r.zeros(d);
      x.zeros(d);
      y.zeros(d);
      z.zeros(d);
      N.zeros(d, d);
      b1.ones(d);

      for (i = 0; i < d-1; i++) {
        jj = i + zzz;
        j0 = find(jj <= 0);
        r(i, 0) = -sum(ppp(j0));
      }

      b.zeros(2 * d - 1);
      b( arma::conv_to<arma::uvec>::from( d + arma::conv_to<arma::ivec>::from(-zzz) - 1)  ) = -ppp;
      b[d-1] = 1 + b[d-1];
      r[0] = 1 + r[0];

      // Initial values
      x[0]  = 1 / r[0];
      y[0]  = 0;
      z[0]  = -b[d-2] * x[0];
      tau   = z[0];
      gamma = b[d-1];
      phi   = 2*b[d-2];

      // Iterations
      for (m = 1; m <= d-1; m++) {
        double d1, d2, psi1, psi2;
        d1 = r[m]*x[0];
        d2 = r[m]*(y[0] - 1);
        if ( m > 1 ) {
          for (i = 1; i <= m-1; i++) {
            d1 += b[m-1-i+d]*x[i];
            d2 += b[m-1-i+d]*y[i];
          }
        }
        gamma += -d1*phi + d2*tau;
        phi = (m < d-1 ? b[d-m-2] : 0);
        for (i = 0; i <= m-1; i++) phi += b[d-i-2] * z[i];
        psi1   = d1/gamma;
        psi2   = d2/gamma;
        tau    = z[0];
        z[m] = 1;

        for (arma::sword j = m; j >= 0; j--) {
          x[j] -= psi1*z[j];
          y[j] -= psi2*z[j];
          z[j] = (j > 0 ? z[j-1] - phi*x[j] + z[0]*y[j] : - phi*x[j] + z[0]*y[j]);
        }
      }

      // Build inverse Matrix
      N.col(0) = x;
      for (j = 1; j <= d-1; j++) {
        f1 = 0;
        f2 = N(d-1, j-1);
        f3 = N(0, j-1);
        for (i = 0; i <= d-2; i++ ) f1 += b[d-i-2]*N(i, j-1);
        N(0, j) = - f1*x[0] - f2*z[0] + f3*y[0];
        for (i = 1; i <= d-1; i++) N(i, j) = N(i-1, j-1) - f1*x[i] - f2*z[i] + f3*y[i];
      }

      // // Solve ARL-equation system
      avg = ( N * b1 );
      arl = avg[0];
      value = arl;
      break;
    }
      // ARL calculation: Brook and Evans 1972 approach
    case 3: {
      // arma::uword t, k, i;
      int i, k, t;
      arma::mat rrr, id;
      arma::colvec b, arl;
      arma::vec jj;
      arma::uvec j0, j1, j2;
      arma::ivec j3;

      t = round(h * scaling);
      rrr.zeros(t, t);               // system matrix
      id.eye(t, t);                  // identity matrix
      b.ones(t);                     // vector of ones

      for (i = 0; i < t; i++ ) {     // fill transition probability matrix
        jj = i + zzz;
        j0 = find(jj <= 0);
        rrr(i, 0) = sum(ppp(j0));
        j1 = find(0 < jj);
        j2 = find(jj < t);
        j3 = vintersection(arma::conv_to<arma::ivec>::from(j1), arma::conv_to<arma::ivec>::from(j2));
        for (k = 0; k < j3.size(); k++) rrr(i, jj(j3(k))) = ppp(j3(k));
      }

      arl = arma::solve(id - rrr, b, arma::solve_opts::fast + arma::solve_opts::no_approx);
      value = arl[0];
      break;
    }
    }
  }
  // simple rounding, Steiner et al. (2000) implementation
  else {
    arma::vec zz, pp, ppp, zzz;

    zz = join_cols(as<arma::vec>(z1), as<arma::vec>(z2));
    pp = join_cols(as<arma::vec>(p1), as<arma::vec>(p2));
    zz = round(zz*scaling);

    ppp = tapply3arma(pp, arma::conv_to<arma::ivec>::from(zz));
    zzz = unique(zz);

    switch (method) {
    // ARL calculation: Toeplitz SPRT approach
    case 1: {
      double al, be, et, ga, de;
      int i, j, N, N1;
      N = round(h * scaling);
      NumericVector x(N), y(N), z(N), phi(N), psi(N), a(2*N - 1), b1(N, 1.), b2(N);
      arma::vec jj;
      arma::uvec j0;

      N1 = N - 1;

      for (i = 0; i < zzz.n_elem; i++) a[N1 - zzz[i]] = -ppp[i];
      a[N1] = 1 + a[N1];

      for (i = 0; i < N1; i++ ) {
        jj = i + zzz;
        j0 = find(jj < 0);
        b2[i] = sum(ppp(j0));
      }

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
        for (j = 1; j < i; j++) z[j] = ( y[j-1] - ga * x[j] ) / be;
        z[i] = y[i-1] / be;

        x[0] = x[0] / be;
        for (j = 1; j < i; j++) x[j] = ( x[j] - al * y[j-1] ) / be;
        x[i] = ( -al*y[i-1] ) / be;

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
      // ARL calculation: Toeplitz full matrix inversion
    case 2: {
    int d, i, j, m;
      double tau, gamma, phi, f1, f2, f3, arl;
      arma::colvec b1, avg;
      arma::vec r, b, jj, x, y, z;
      arma::mat N;
      arma::uvec j0;

      d = round(h * scaling);
      r.zeros(d);
      x.zeros(d);
      y.zeros(d);
      z.zeros(d);
      N.zeros(d, d);
      b1.ones(d);

      for (i = 0; i < d-1; i++) {
        jj = i + zzz;
        j0 = find(jj <= 0);
        r(i, 0) = -sum(ppp(j0));
      }

      b.zeros(2 * d - 1);
      b( arma::conv_to<arma::uvec>::from( d + arma::conv_to<arma::ivec>::from(-zzz) - 1)  ) = -ppp;
      b[d-1] = 1 + b[d-1];
      r[0] = 1 + r[0];

      // Initial values
      x[0]  = 1 / r[0];
      y[0]  = 0;
      z[0]  = -b[d-2] * x[0];
      tau   = z[0];
      gamma = b[d-1];
      phi   = 2*b[d-2];

      // Iterations
      for (m = 1; m <= d-1; m++) {
        double d1, d2, psi1, psi2;
        d1 = r[m]*x[0];
        d2 = r[m]*(y[0] - 1);
        if ( m > 1 ) {
          for (i = 1; i <= m-1; i++) {
            d1 += b[m-1-i+d]*x[i];
            d2 += b[m-1-i+d]*y[i];
          }
        }
        gamma += -d1*phi + d2*tau;
        phi = (m < d-1 ? b[d-m-2] : 0);
        for (i = 0; i <= m-1; i++) phi += b[d-i-2] * z[i];
        psi1   = d1/gamma;
        psi2   = d2/gamma;
        tau    = z[0];
        z[m] = 1;

        for (arma::sword j = m; j >= 0; j--) {
          x[j] -= psi1*z[j];
          y[j] -= psi2*z[j];
          z[j] = (j > 0 ? z[j-1] - phi*x[j] + z[0]*y[j] : - phi*x[j] + z[0]*y[j]);
        }
      }

      // Build inverse Matrix
      N.col(0) = x;
      for (j = 1; j <= d-1; j++) {
        f1 = 0;
        f2 = N(d-1, j-1);
        f3 = N(0, j-1);
        for (i = 0; i <= d-2; i++ ) f1 += b[d-i-2]*N(i, j-1);
        N(0, j) = - f1*x[0] - f2*z[0] + f3*y[0];
        for (i = 1; i <= d-1; i++) N(i, j) = N(i-1, j-1) - f1*x[i] - f2*z[i] + f3*y[i];
      }

      // // Solve ARL-equation system
      avg = ( N * b1 );
      arl = avg[0];
      value = arl;
      break;
    }
      // ARL calculation: Brook and Evans 1972 approach
    case 3: {
      int i, k, t;
      arma::mat rrr, id;
      arma::colvec b, arl;
      arma::vec jj;
      arma::uvec j0, j1, j2;
      arma::ivec j3;

      t = round(h * scaling);
      rrr.zeros(t, t);               // system matrix
      id.eye(t, t);                  // identity matrix
      b.ones(t);                     // vector of ones

      for (i = 0; i < t; i++ ) {                  // fill transition probability matrix
        jj = i + zzz;
        j0 = find(jj <= 0);
        rrr(i, 0) = sum(ppp(j0));
        j1 = find(0 < jj);
        j2 = find(jj < t);
        j3 = vintersection(arma::conv_to<arma::ivec>::from(j1), arma::conv_to<arma::ivec>::from(j2));
        for (k = 0; k < j3.size(); k++) rrr(i, jj(j3(k))) = ppp(j3(k));
      }
      arl = arma::solve(id - rrr, b, arma::solve_opts::fast + arma::solve_opts::no_approx);
      value = arl[0];
      break;
    }
    }
  }
  return (value);
}

// [[Rcpp::export(.racusum_crit_mc)]]
double racusum_crit_mc(NumericMatrix pmix, double L0, double RA, double R, double scaling, int rounding, int method, int jmax, bool verbose) {
  double L1, h, h1;
  int i;
  for (i = 1; i < 10; i++ ) {
    L1 = racusum_arl_mc(pmix, RA, 1, double(i), scaling, rounding, method);
    if ( verbose ) Rcpp::Rcout << "h = " <<  i << "\t" << "ARL = " << L1 << std::endl;
      if ( L1 > L0 ) break;
  }
  h1 = i;

  for (int j = 1; j <= jmax; j++ ) {
    for (int dh = 1; dh <= 19; dh++ ) {
      h = h1 + pow(static_cast<double>(-1), j) * dh / pow(static_cast<double>(10), j);
      L1 = racusum_arl_mc(pmix, RA, 1, h, scaling, rounding, method);
      if ( verbose ) Rcpp::Rcout << "h = " <<  h << "\t" << "ARL = " << L1 << std::endl;
      if ( ((j % 2 == 1) & (L1 < L0) ) | ((j % 2 == 0) & (L1 > L0)) ) break;
    }
    h1 = h;
  }
  if ( L1 < L0 ) h = h + 1 / pow(static_cast<double>(10), jmax);
  return h;
}
