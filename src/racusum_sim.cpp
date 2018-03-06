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

double llr_score_oc(DataFrame df, NumericVector coeff, NumericVector coeff2, double R0, double RA, double RQ){
  int y, row, s;
  double wt, x, Qstar, xstar, rdm, pt;
  NumericVector col1, col2, rnd, rndm;

  col1 = df[0];                                       // Parsonnet score stored in the dataframe
  col2 = df[1];                                       // Surgical outcome stored in the dataframe
  rnd = runif(1);
  row = rnd[0] * df.nrows();
  s = col1[row];                                      // Step 1: s=Parsonnet score
  x = gettherisk(s, coeff);                           // Step 2:
  Qstar = RQ;                                         // Step 3:
  xstar = (Qstar*x) / (1-x+Qstar*x);
  rndm = runif(1);
  rdm = as<double>(rndm);
  y = (rdm < xstar ? 1 : 0);                          // Step 4: y=Surgical outcome
  pt = gettherisk(s, coeff2);                         // Step 5: x=True probability of death
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

// [[Rcpp::export(.cusum_arl_sim)]]
int cusum_arl_sim(int r, double h, DataFrame df, double R0, double RA) {
  double qn = 0, wt = 0;
  int rl = 0;
  do{
    rl++;
    wt = llr_score_noadjust(df, R0, RA);
    qn = fmax(0, qn + wt);
  } while (qn <= h);
  return rl;
}

// [[Rcpp::export(.racusum_arl_sim)]]
int racusum_arl_sim(int r, NumericVector coeff, double h, DataFrame df, double R0, double RA, bool yemp) {
  double qn = 0, wt = 0;
  int rl = 0;
  do{
    rl++;
    wt = llr_score(df, coeff, R0, RA, yemp);
    qn = fmax(0, qn + wt);
  } while (qn <= h);
  return rl;
}

// [[Rcpp::export(.racusum_arloc_sim)]]
int racusum_arloc_sim(int r, NumericVector coeff, NumericVector coeff2, double h, DataFrame df, double R0, double RA, double RQ) {
  double qn = 0, wt = 0;
  int rl = 0;
  do{
    rl++;
    wt = llr_score_oc(df, coeff, coeff2, R0, RA, RQ);
    qn = fmax(0, qn + wt);
  } while (qn <= h);
  return rl;
}

// [[Rcpp::export(.racusum_adoc_sim)]]
int racusum_adoc_sim(int r, NumericVector coeff, NumericVector coeff2, double h, DataFrame df, double R0, double RA,  double RQ, int m, int type) {
  // conditional steady-state ARL (RA-CUSUM) -- m = #ic-observations with m >= 0
  int rl = 0;
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
        wt = llr_score_oc(df, coeff, coeff2, R0, RA, R);
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
      if ( rl > m ) R = RQ;
      wt = llr_score_oc(df, coeff, coeff2, R0, RA, R);
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
