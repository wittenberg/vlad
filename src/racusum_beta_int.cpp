#include "racusum_beta_int.h"

// [[Rcpp::export(.Tn)]]
double Tn(double z, int j) {return( cos( (j-1) * acos(z) ) );}

//* PDF of W */
// [[Rcpp::export(.f2)]]
double f2(double w, double RA, double RQ, double g0, double g1, double shape1, double shape2) {
  double w0, w1, w2, w3, res, logRA=log(RA);

  if (RA > 1) {   /* Decting deterioration QA > 1 */
    w0 = -log(1 + (RA-1) * hS(1, g0, g1, 1));           /* lower limit left side of fW */
    w1 = -log(1 + (RA-1) * hS(0, g0, g1, 1));           /* upper limit left side of fW */
    w2 = -log(1 + (RA-1) * hS(1, g0, g1, 1)) + logRA; /* lower limit right side of fW */
    w3 = -log(1 + (RA-1) * hS(0, g0, g1, 1)) + logRA; /* upper limit right side of fW */
    if ((w0 <= w) & (w <= w1)) {
      res = -(1 - hS(s0(w, RA, g0, g1, RQ), g0, g1, 1)) *
            R::dbeta(s0(w, RA, g0, g1, 1), shape1, shape2, false) *
            s0p(w, RA, g0, g1, 1) ;}
    else if ((w2 <= w) & (w <= w3)) {
      res = -hS(s1(w, RA, g0, g1, RQ), g0, g1, 1) *
            R::dbeta(s1(w, RA, g0, g1, 1), shape1, shape2, false) *
            s1p(w, RA, g0, g1, 1);}
    else {res = 0;}
  } else {   /* Decting improvement QA > 0 and QA < 1 */
    w0 = -log(1 + (RA-1) * hS(0, g0, g1, 1)) + logRA; /* lower limit left side of fW */
    w1 = -log(1 + (RA-1) * hS(1, g0, g1, 1)) + logRA; /* upper limit left side of fW */
    w2 = -log(1 + (RA-1) * hS(0, g0, g1, 1));           /* lower limit right side of fW */
    w3 = -log(1 + (RA-1) * hS(1, g0, g1, 1));           /* upper limit right side of fW */
    if ((w0 <= w) & (w <= w1)) {
      res = hS(s1(w, RA, g0, g1, RQ), g0, g1, 1) *
            R::dbeta(s1(w, RA, g0, g1, 1), shape1, shape2, false) *
            s1p(w, RA, g0, g1, 1);}
    else if ((w2 <= w) & (w <= w3)) {
      res = (1 - hS(s0(w, RA, g0, g1, RQ), g0, g1, 1)) *
            R::dbeta(s0(w, RA, g0, g1, 1), shape1, shape2, false) *
            s0p(w, RA, g0, g1, 1) ;}
    else {res = 0;}
  }
  return(res);
}

// [[Rcpp::export(.integ_t62)]]
double integ_t62(double xl, double xu, int j, double loc, double RA, double RQ, double g0, double g1, double shape1, double shape2) {
  double w0 = 0, w1 = 0, w2 = 0, w3 = 0, ires = 0, I1l, I1u, I2l, I2u, ires1=0, ires2=0, lower, upper, logRA=log(RA);

  if (RA > 1) {   /* Detecting deterioration QA > 1 */
    w0 = -log(1 + (RA-1) * hS(1, g0, g1, RQ));           /* lower limit left side of fW */
    w1 = -log(1 + (RA-1) * hS(0, g0, g1, RQ));           /* upper limit left side of fW */
    w2 = -log(1 + (RA-1) * hS(1, g0, g1, RQ)) + logRA;   /* lower limit right side of fW */
    w3 = -log(1 + (RA-1) * hS(0, g0, g1, RQ)) + logRA;   /* upper limit right side of fW */

    if ( loc + w3 < xl ) {return(ires = 0.);}
    if ( loc + w0 > xu ) {return(ires = 0.);}

    I1l = max((loc + w0), xl); /* lower integ. limit left side of f_W */
    I1u = min((loc + w1), xu); /* upper integ. limit left side of f_W */
    I2l = max((loc + w2), xl); /* lower integ. limit right side of f_W */
    I2u = min((loc + w3), xu); /* lower integ. limit right side of f_W */

    /* Gauss-Kronrod quadrature (adaptive) from boost  */
    /* apply quadratic transformation (sqrt(limits) and f(y*y) * 2*y) */
    auto f01 = [xl, xu, j, loc, RA, RQ, g0, g1, shape1, shape2, w1](double y) {
      return Tn((2.*(w1-pow(static_cast<double>(y), 2)+loc)-xu-xl)/(xu-xl), j)*f2((w1-pow(static_cast<double>(y), 2)), RA, RQ, g0, g1, shape1, shape2)*-2.*y;
    };
    auto f23 = [xl, xu, j, loc, RA, RQ, g0, g1, shape1, shape2, w3](double y) {
      return Tn((2.*(w3-pow(static_cast<double>(y), 2)+loc)-xu-xl)/(xu-xl), j)*f2((w3-pow(static_cast<double>(y), 2)), RA, RQ, g0, g1, shape1, shape2)*-2.*y;
    };
    I1l -= loc, I1u -= loc, I2l -= loc, I2u -= loc;
    if ((loc + w1 >= xl) & (I1u > I1l)) {   /* left side of f_W */
      if (w1 - I1u <= 0) {lower = 0;} else {lower = pow(static_cast<double>(w1 - I1u), .5);}
      upper = pow(static_cast<double>(w1 - I1l), .5);
      ires1 = gauss_kronrod<double, 61>::integrate(f01, lower, upper, 10, 1e-4);
    }
    if ((loc + w2 <= xu) & (I2u > I2l)) {   /* right side of f_W */
      if (w3 - I2u <= 0) {lower = 0;} else {lower = pow(static_cast<double>(w3 - I2u), .5);}
      upper = pow(static_cast<double>(w3 - I2l), .5);
      ires2 = gauss_kronrod<double, 61>::integrate(f23, lower, upper, 10, 1e-4);
    }
    ires -= ires1 + ires2;
  } else if (RA < 1) {   /* Detecting improvement QA > 0 and QA < 1 */
    w0 = -log(1 + (RA-1) * hS(0, g0, g1, RQ)) + logRA;   /* lower limit left side of fW */
    w1 = -log(1 + (RA-1) * hS(1, g0, g1, RQ)) + logRA;   /* upper limit left side of fW */
    w2 = -log(1 + (RA-1) * hS(0, g0, g1, RQ));           /* lower limit right side of fW */
    w3 = -log(1 + (RA-1) * hS(1, g0, g1, RQ));           /* upper limit right side of fW */

    if ( loc + w3 < xl ) {return(ires = 0.);}
    if ( loc + w0 > xu ) {return(ires = 0.);}

    I1l = max((loc + w0), xl); /* lower integ. limit left side of f_W */
    I1u = min((loc + w1), xu); /* upper integ. limit left side of f_W */
    I2l = max((loc + w2), xl); /* lower integ. limit right side of f_W */
    I2u = min((loc + w3), xu); /* lower integ. limit right side of f_W */

    /* Gauss-Kronrod quadrature (adaptive) from boost  */
    /* apply quadratic transformation (sqrt(limits) and f(y*y) * 2*y) */
    auto f01 = [xl, xu, j, loc, RA, RQ, g0, g1, shape1, shape2, w0](double y) {
      return Tn((2.*(w0+pow(static_cast<double>(y), 2)+loc)-xu-xl)/(xu-xl), j)*f2((w0+pow(static_cast<double>(y), 2)), RA, RQ, g0, g1, shape1, shape2)*2.*y;
    };
    auto f23 = [xl, xu, j, loc, RA, RQ, g0, g1, shape1, shape2, w2](double y) {
      return Tn((2.*(w2+pow(static_cast<double>(y), 2)+loc)-xu-xl)/(xu-xl), j)*f2((w2+pow(static_cast<double>(y), 2)), RA, RQ, g0, g1, shape1, shape2)*2.*y;
    };
    I1l -= loc, I1u -= loc, I2l -= loc, I2u -= loc;
    if ((loc + w1 >= xl) & (I1u > I1l)) {   /* left side of f_W */
      if (I1l - w0 <= 0) {lower = 0;} else {lower = pow(static_cast<double>(I1l - w0), .5);}
      upper = pow(static_cast<double>(I1u - w0), .5);
      ires1 = gauss_kronrod<double, 61>::integrate(f01, lower, upper, 10, 1e-4);
    }
    if ((loc + w2 <= xu) & (I2u > I2l)) {   /* right side of f_W */
      if (I2l - w2 <= 0) {lower = 0;} else {lower = pow(static_cast<double>(I2l - w2), .5);}
      upper = pow(static_cast<double>(I2u - w2), .5);
      ires2 = gauss_kronrod<double, 61>::integrate(f23, lower, upper, 10, 1e-4);
    }
    ires += ires1 + ires2;
  }
  return(ires);
}

// [[Rcpp::export(.racusum_beta_arl_int)]]
double racusum_beta_arl_int(double h, int N, double RA, double RQ, double g0, double g1, double shape1, double shape2, bool pw) {
  double arl = 0., w0=0, w3=0, M1=0, M2=0, M=0, za=0, zb=0, xl=0, xu=0, logRA=log(RA);
  int i=0, j=0, ii=0, jj=0, dN=0, NN=0, Ntilde=0, qi=0, qj=0;

  if (pw == TRUE) { /* Piece-wise collocation */

    /* Decting deterioration RA > 1 */
    if (RA > 1.) {
      w3 = abs(-log(1. + (RA-1) * hS(0, g0, g1, RQ)) + logRA);
      M1 = ceil(h/w3) + 1;
      NumericVector b11(M1);

      for (i = 0; i < M1; i++) b11[i] = h - i*w3;      /* interval borders death */
      b11[M1-1] = 0.;

      M = M1 - 1;                                  /* merge intervals */

      NumericVector b(M);
      b = sort(as<arma::vec>(b11));
      Ntilde = static_cast<int>(ceil( static_cast<double>(N) / M ));
      dN     = Ntilde;
      NN     = M*Ntilde;
      NumericVector c1(N);

      NumericMatrix zn(M, Ntilde);
      arma::mat m(NN, NN);
      m.zeros(N, N);
      arma::colvec g(NN);
      g.ones(NN);
      Rcout << "Piece-wise collocation: " << M << " intervals" << std::endl;

      /* Chebyshev nodes on [b_1,b_2],[b_2,b_3],...,[b_M,hu] */
      for (i = 1; i <= M; i++) {                 /* intervals */
        for (j = 1; j <= Ntilde; j++) {          /* order of polynom */
          zn(i-1, j-1) = b[i-1] + ( b[i]-b[i-1] ) / 2. * ( 1. + cos( (2.*(Ntilde-j+1.)-1.)*M_PI   /  2./ dN ) );
        }
      }
      /* Solve system of linear equations, see ... */
      for(i = 1; i <= M; i++) {                 /* 1. Intervals */
        for(j = 1; j <= Ntilde; j++) {          /* 2. Order of polynom */
          qi = (i-1)*Ntilde + j-1;
          zb = zn(i-1, j-1) + w3;
          for(ii = 1; ii <= M; ii++) {          /* 3. Intervals */
            xl = b[ii-1];
            xu = min(b[ii], zb);
            for(jj = 1; jj <= Ntilde; jj++) {   /* 4. Order of polynom */
              qj = (ii-1)*Ntilde + jj-1;
              if ( xu < xl ) {m(qi, qj) = 0.;} else {
              m(qi, qj) = -integ_t62(b[ii-1], b[ii], jj, zn(i-1, j-1), RA, RQ, g0, g1, shape1, shape2);
              }
            }
          }
        for (jj = 1; jj <= Ntilde; jj++) {
          qj = (i-1)*Ntilde + jj-1;
          m(qi, jj-1) -= Tn(-1., jj) * FWT2(-zn(i-1,j-1), RA, g0, g1, shape1, shape2, RQ);
          m(qi, qj)   += Tn( (2.*zn(i-1, j-1)-b[i]-b[i-1]) / (b[i]-b[i-1]), jj );
        }
      }
    }
    c1 = arma::solve(m, g);                                   /* Unknown constants */
    for (j = 0; j < Ntilde; j++) arl += c1(j) * Tn(-1., j+1); /* Approximate ARL by N indep. interpol. func. */
    return(arl);
    }
    /* Decting improvement 0 < RA < 1 */
    else if (RA < 1.) {
      w0 = abs(-log(1 + (RA-1) * hS(0, g0, g1, RQ)) + logRA);
      M1 = 0;
      M2 = ceil(h/w0)+1;
      NumericVector b11(M1), b12(M2);
      b12[M2-1] = h;

      for (i = 1; i < (M2-1); i++) b12[i] = i*w0;          /* interval borders death */

      M = M1 + M2 - 1;                                  /* merge intervals */
      NumericVector b(M);
      b = sort(as<arma::vec>(b12));
      Ntilde = static_cast<int>(ceil( static_cast<double>(N) / M ));
      dN     = Ntilde;
      NN     = M*Ntilde;
      NumericVector c1(N);

      NumericMatrix zn(M, Ntilde);
      arma::mat m(NN, NN);
      m.zeros(N, N);
      arma::colvec g(NN);
      g.ones(NN);
      Rcout << "Piece-wise collocation: " << M << " intervals" << std::endl;
      /* Chebyshev nodes on [b_1,b_2],[b_2,b_3],...,[b_M,hu] */
      for (i = 1; i <= M; i++) {                 /* intervals */
        for (j = 1; j <= Ntilde; j++) {          /* order of polynom */
          zn(i-1, j-1) = b[i-1] + ( b[i]-b[i-1] ) / 2. * ( 1. + cos( (2.*(Ntilde-j+1.)-1.)*M_PI   /  2./ dN ) );
        }
      }
      /* Solve system of linear equations, see ... */
      for(i = 1; i <= M; i++) {                 /* 1. Intervals */
        for(j = 1; j <= Ntilde; j++) {          /* 2. Order of polynom */
          qi = (i-1)*Ntilde + j-1;
          za = zn(i-1, j-1) - w0;
          for(ii = 1; ii <= M; ii++) {          /* 3. Intervals */
            xl = max(b[ii-1], za);
            xu = b[ii];
            for(jj = 1; jj <= Ntilde; jj++) {   /* 4. Order of polynom */
              qj = (ii-1)*Ntilde + jj-1;
              if ( xu < xl ) {m(qi, qj) = 0.;} else {
              m(qi, qj) = -integ_t62(b[ii-1], b[ii], jj, zn(i-1, j-1), RA, RQ, g0, g1, shape1, shape2);
              }
            }
          }
          for (jj = 1; jj <= Ntilde; jj++) {
            qj = (i-1)*Ntilde + jj-1;
            m(qi, jj-1) -= Tn(-1., jj) * FWT2(-zn(i-1,j-1), RA, g0, g1, shape1, shape2, RQ);
            m(qi, qj)   += Tn( (2.*zn(i-1, j-1)-b[i]-b[i-1]) / (b[i]-b[i-1]), jj );
          }
        }
      }
      c1 = arma::solve(m, g);                                   /* Unknown constants */
      for (j = 0; j < Ntilde; j++) arl += c1(j) * Tn(-1., j+1); /* Approximate ARL by N indep. interpol. func. */
      return(arl);
    }
  } else { /* Plain collocation */
    double arl = 0.;
    int i, j;
    NumericVector z(N), c1(N);
    arma::mat m(N, N);
    arma::colvec b(N);

    b.ones(N);
    m.zeros(N, N);
    Rcout << "Collocation: " << 1 << " interval" << std::endl;
    /* Chebyshev nodes in the domain [u0, h] */
    for (i = 0; i < N; i++) z[i] = h / 2. * ( 1.+cos((2.*(i+1)-1)*M_PI / (2.*N)) );
    /* Solve system of linear equations, see E.1 in Tang et al. (2015) p. 25 */
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        m(i, j) =                                                        /* System matrix */
          Tn((2.*z[i]-h-0)/(h-0), j+1) -                                 /* Chebyshev polynomials */
          Tn(-1., j+1) * FWT2(-z[i], RA, g0, g1, shape1, shape2, RQ) -   /* CDF of W */
          integ_t62(0, h, j+1, z[i], RA, RQ, g0, g1, shape1, shape2);    /* PDF of W */
      }
    }
    c1 = arma::solve(m, b);                              /* Unknown constants */
    for (j = 0; j < N; j++) arl += c1(j) * Tn(-1., j+1); /* Approximate ARL by N indep. interpol. func. */
    return(arl);
  }
  return(arl);
}
