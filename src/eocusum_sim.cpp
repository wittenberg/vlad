#include "eocusum_sim.h"

// [[Rcpp::export(.optimal_k)]]
double optimal_k(DataFrame pmix, double RA, bool yemp) {
  int n;
  double optk = 0, sum = 0, pbar;
  NumericVector y, pt;

  y = pmix[0];
  pt = pmix[1];
  n = pmix.nrows();

  if (yemp == true) {
    for (int i=0; i < n; ++i) {sum += y[i];}
    pbar = sum/n;
  } else if (yemp == false) {
    for (int i=0; i < n; ++i) {sum += pt[i];}
    pbar = sum/n;
  }
  if (RA >= 1) { optk = pbar * ( RA - 1 - log(RA) ) / log(RA); }   // Detrerioration
  else if ( (RA > 0) & (RA < 1) ) {
    optk = pbar * ( 1 - RA + log(RA) ) / log(RA); }                 // Improvement
  return optk;
}

// [[Rcpp::export(.eocusum_arl_sim)]]
int eocusum_arl_sim(int r, DataFrame pmix, double k, double h, double RQ, bool yemp, int side) {
  double z = 0, x, pistar, pt;
  int y, row, rl = 0, pmixnrows = pmix.nrows();
  NumericVector pi1, pi2, pi3;
  pi1 = pmix[0];                        // Empirical data
  pi2 = pmix[1];                        // 2st Predicted probability
  pi3 = pmix[2];                        //
  if (side == 1) {                                    // lower side (deterioration)
    double tn = 0;
    do{
      rl++;
      row = floor(runif(1, 0, pmixnrows)[0]);
      x = pi2[row];                    // probability of death
      pistar = (RQ * x) / (1 - x + RQ * x);
      if ((yemp == true) & (RQ == 1)) {
        y = pi1[row];
        pt = pistar;
      } else {                                  // Surgical outcome empirical
        y = (R::runif(0, 1) < pistar ? 1 : 0);  // Step 4: y=Surgical outcome
        pt = pi3[row];
      }
      z = pt - y;
      tn = fmin(0, tn + z + k);
    } while (-tn <= h);
    return rl;
  }
  else {                                // upper side (improvement)
    double qn = 0;
    do{
      rl++;
      row = floor(runif(1, 0, pmixnrows)[0]);
      x = pi2[row];                    // probability of death
      pistar = (RQ * x) / (1 - x + RQ * x);
      if (yemp == true) {
        y = pi1[row];
        pt = pistar;
      } else {                                  // Surgical outcome empirical
        y = (R::runif(0, 1) < pistar ? 1 : 0);  // Step 4: y=Surgical outcome
        pt = pi3[row];
      }
      z = pt - y;
      qn = fmax(0, qn + z - k);
    } while (qn <= h);
    return rl;
  }
}

// [[Rcpp::export(.eocusum_ad_sim)]]
int eocusum_ad_sim(int r, DataFrame pmix, double k, double h, double RQ, int side, int type, int m) {
  if (type == 1) {                                                  // conditional steady-state ARL (EO-CUSUM)
    if (side == 1) {
      return eocusum_ad_sim11(r, pmix, k, h, RQ, m); // lower side (deterioration)
    } else {
      return eocusum_ad_sim12(r, pmix, k, h, RQ, m); // upper side (improvement)
    }
  } else {                                                          // cyclical steady-state ARL (EO-CUSUM)
    if (side == 2) {
      return eocusum_ad_sim21(r, pmix, k, h, RQ, m); // lower side (deterioration)
    } else {
      return eocusum_ad_sim22(r, pmix, k, h, RQ, m); // upper side (improvement)
    }
  }
}

// conditional steady-state ARL (EO-CUSUM) -- m = #ic-observations with m >= 0
// lower side (deterioration)
int eocusum_ad_sim11(int r, DataFrame pmix, double k, double h, double RQ, int m) {
  NumericVector pi2, pi3;
  int success = 0, rl = 0, y, row, pmixnrows = pmix.nrows();
  double z = 0, x, pistar, pt;
  double tn = 0, R = 1;
  pi2 = pmix[1];                                 // 2st Predicted probability
  pi3 = pmix[2];                                 // 2st Predicted probability
  while ( !success ) {
    rl = 0;
    tn = 0;
    R = 1;
    do{
      rl++;
      if ( rl > m) R = RQ;
      row = floor(runif(1, 0, pmixnrows)[0]);
      x = pi2[row];                               // probability of death
      pistar = (R * x) / (1 - x + R * x);
      y = (R::runif(0, 1) < pistar ? 1 : 0);      // Step 4: y=Surgical outcome
      pt = pi3[row];
      z = pt - y;
      tn = fmin(0, tn + z + k);
    } while (-tn <= h);
    if ( rl > m ) {
      rl += -m;
      success = 1;
    }
  }
  return rl;
}

// conditional steady-state ARL (EO-CUSUM) -- m = #ic-observations with m >= 0
// upper side (improvement)
int eocusum_ad_sim12(int r, DataFrame pmix, double k, double h, double RQ, int m) {
  NumericVector pi2, pi3;
  int success = 0, rl = 0, y, row, pmixnrows = pmix.nrows();
  double z = 0, x, pistar, pt;
  double qn = 0, R = 1;
  pi2 = pmix[1];                                 // 2st Predicted probability
  pi3 = pmix[2];                                 // 2st Predicted probability
  while ( !success ) {
    rl = 0;
    qn = 0;
    R = 1;
    do{
      rl++;
      if ( rl > m) R = RQ;
      row = floor(runif(1, 0, pmixnrows)[0]);
      x = pi2[row];                               // probability of death
      pistar = (R * x) / (1 - x + R * x);
      y = (R::runif(0, 1) < pistar ? 1 : 0);             // Step 4: y=Surgical outcome
      pt = pi3[row];
      z = pt - y;
      qn = fmax(0, qn + z - k);
    } while (qn <= h);
    if ( rl > m ) {
      rl += -m;
      success = 1;
    }
  }
  return rl;
}

// cyclical steady-state ARL (EO-CUSUM) -- m = #ic-observations with m >= 0
// lower side (deterioration)
int eocusum_ad_sim21(int r, DataFrame pmix, double k, double h, double RQ, int m) {
  NumericVector pi2, pi3;
  int rl = 0, y, row, pmixnrows = pmix.nrows();
  double z = 0, x, pistar, pt;
  double tn = 0, R = 1;
  pi2 = pmix[1];                                 // 2st Predicted probability
  pi3 = pmix[2];                                 // 2st Predicted probability
  rl = 0;
  tn = 0;
  R = 1;
  do{
    rl++;
    if ( rl > m) R = RQ;
    row = floor(runif(1, 0, pmixnrows)[0]);
    x = pi2[row];                               // probability of death
    pistar = (R * x) / (1 - x + R * x);
    y = (R::runif(0, 1) < pistar ? 1 : 0);             // Step 4: y=Surgical outcome
    pt = pi3[row];
    z = pt - y;
    tn = fmin(0, tn + z + k);
    if ( rl <= m ) {
      if ( -tn > h) {
        tn = 0;
      }
    }
  } while (-tn <= h);
  rl += -m ;
  return rl;
}

// cyclical steady-state ARL (EO-CUSUM) -- m = #ic-observations with m >= 0
// upper side (improvement)
int eocusum_ad_sim22(int r, DataFrame pmix, double k, double h, double RQ, int m) {
  NumericVector pi2, pi3;
  int rl = 0, y, row, pmixnrows = pmix.nrows();
  double z = 0, x, pistar, pt;
  double qn = 0, R = 1;
  pi2 = pmix[1];                                 // 2st Predicted probability
  pi3 = pmix[2];                                 // 2st Predicted probability
  rl = 0;
  qn = 0;
  R = 1;
  do{
    rl++;
    if ( rl > m) R = RQ;
    row = floor(runif(1, 0, pmixnrows)[0]);
    x = pi2[row];                               // probability of death
    pistar = (R * x) / (1 - x + R * x);
    y = (R::runif(0, 1) < pistar ? 1 : 0);             // Step 4: y=Surgical outcome
    pt = pi3[row];
    z = pt - y;
    qn = fmax(0, qn + z - k);
    if ( rl <= m ) {
      if ( qn > h) {
        qn = 0;
      }
    }
  } while (qn <= h);
  rl += -m ;
  return rl;
}
