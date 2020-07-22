#include "racusum_sim.h"

// [[Rcpp::export(.llr_score)]]
double llr_score(DataFrame df, NumericVector coeff, double R0, double RA, bool yemp){
  int y, row, s;
  double wt, pt;
  NumericVector col1, col2, rnd;

  col1 = df[0];                                       // Parsonnet score stored in the dataframe
  col2 = df[1];                                       // Surgical outcome stored in the dataframe
  rnd = runif(1);
  row = rnd[0] * df.nrows();
  s = col1[row];                                      // Step 1: s=Parsonnet score
  pt = gettherisk(s, coeff);                          // Step 2
  if (yemp == true) {
    y = col2[row];
  }
  else {
    NumericVector rndm = runif(1);
    double rdm = as<double>(rndm);
    y = (rdm < pt ? 1 : 0);
  }
  wt = (  y == 1 ? log( ((1 - pt + R0*pt)*RA) / ((1 - pt + RA*pt)*R0) ) : log( (1 - pt + R0*pt) / (1 - pt + RA*pt) )  );
  return wt;
}


double llr_score_noadjust(DataFrame df, double R0, double RA){
  NumericVector col2, rnd;
  int row, y;
  double wt;

  col2 = df[1];                                       // Surgical outcome stored in the dataframe
  rnd = runif(1);
  row = rnd[0] * df.nrows();
  y = col2[row];
  wt = (  y == 1 ? log( RA / R0 ) : log( (1 - RA) / (1 - R0) )  );
  return wt;
}

// [[Rcpp::export(.bcusum_arl_sim)]]
int bcusum_arl_sim(int r, double h, DataFrame df, double R0, double RA) {
  double qn = 0, wt = 0;
  int rl = 0;
  do{
    rl++;
    wt = llr_score_noadjust(df, R0, RA);
    qn = fmax(0, qn + wt);
  } while (qn <= h);
  return rl;
}

// [[Rcpp::export(.racusum_ad_sim)]]
int racusum_ad_sim(int r, DataFrame pmix, double h, double RA, double RQ, int m, int type) {
  // conditional steady-state ARL (RA-CUSUM) -- m = #ic-observations with m >= 0
  int rl = 0, y, row, pmixnrows = pmix.nrows();
  NumericVector pi2, pi3;
  double x, pistar, pt, logRA = log(RA);
  pi2 = pmix[1];                                 // 2st Predicted probability
  pi3 = pmix[2];                                 // 2st Predicted probability
  if (type == 1) {
    int success = 0, rl = 0;
    double qn = 0, wt = 0, R = 1;
    while ( !success ) {
      rl = 0;
      qn = 0;
      R = 1;
      do {
        rl++;
        if ( rl > m) R = RQ;
        row = floor(runif(1, 0, pmixnrows)[0]);
        x = pi2[row];                           // True probability of death
        pistar = (R * x) / (1 - x + R * x);
        y = (R::runif(0, 1) < pistar ? 1 : 0);             // Step 4: y=Surgical outcome
        pt = pi3[row];
        wt = -log(1 - pt + RA * pt) + y * logRA;
        qn = fmax(0, qn + wt);
      } while ( qn <= h );
      if ( rl > m ) {
        rl += -m;
        success = 1;
      }
    }
    return rl;
  } else if (type == 2) {
    // cyclical steady-state ARL (RA-CUSUM) -- m = #ic-observations with m >= 0
    int rl = 0;
    double qn = 0, wt = 0, R = 1;
    do {
      rl++;
      if ( rl > m) R = RQ;
      row = floor(runif(1, 0, pmixnrows)[0]);
      x = pi2[row];                           // True probability of death
      pistar = (R * x) / (1 - x + R * x);
      y = (R::runif(0, 1) < pistar ? 1 : 0);             // Step 4: y=Surgical outcome
      pt = pi3[row];
      wt = -log(1 - pt + RA * pt) + y * logRA;
      qn = fmax(0, qn + wt);
      if ( rl <= m ) {
        if ( qn > h ) {
          qn = 0;
        }
      }
    } while (qn <= h);
    rl += -m;
    return rl;
  }
  return rl;
}

// [[Rcpp::export(.racusum_arl_sim)]]
int racusum_arl_sim(int r, DataFrame pmix, double h, double RA, double RQ, bool yemp) {
  double qn = 0, wt = 0, x, pistar, pt, logRA = log(RA);
  int y, rl = 0, pmixnrows = pmix.nrows(), row;
  NumericVector pi1, pi2, pi3;
  pi1 = pmix[0];                        // Empirical data
  pi2 = pmix[1];                        // 2st Predicted probability
  pi3 = pmix[2];                        // 2st Predicted probability
  do {                                   // In-/Out-of-Control ARL
    rl++;
    row = floor(runif(1, 0, pmixnrows)[0]);
    x = pi2[row];                       // True probability of death
    pistar = (RQ * x) / (1 - x + RQ * x);
    if ((yemp == true) & (RQ == 1)) {     // Only for In-control ARL (emp/model)
      y = pi1[row];
      pt = pistar;
      } else {
      y = (R::runif(0, 1) < pistar ? 1 : 0);       // Surgical outcome
      pt = pi3[row];
    }
    wt = -log(1 - pt + RA * pt) + y * logRA;
    qn = fmax(0, qn + wt);
  } while (qn <= h);
  return rl;
}
