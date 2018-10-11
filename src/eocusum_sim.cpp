#include "eocusum_sim.h"

// [[Rcpp::export(.gettherisk)]]
double gettherisk(int parsonnetscore, NumericVector coeff) {
  double logitp, risk;

  logitp = coeff[0] + parsonnetscore * coeff[1];
  risk = exp(logitp)/(1 + exp(logitp));
  return risk;
}

// [[Rcpp::export(.optimal_k)]]
double optimal_k(double QA, DataFrame df, NumericVector coeff, bool yemp) {
  int n;
  double optk = 0, sum = 0, pbar;
  IntegerVector parsonnetscores = df[0];
  NumericVector outcome = df[1];

  n = df.nrows();
  if (yemp == true) {
    for (int i=0; i < n; ++i) {sum += outcome[i];}
    pbar = sum/n;
  } else if (yemp == false) {
  //} else if (yemp == false & coeff.isNotNull()) { NULL not working yet
    for (int i=0; i < n; ++i) {sum += gettherisk(parsonnetscores[i], coeff);}
    pbar = sum/n;
  }
  if (QA > 1) {optk = pbar * ( QA - 1 - log(QA) ) / log(QA);}                     // Detrerioration
  else if ( (QA > 0) & (QA < 1) ){optk = pbar * ( 1 - QA + log(QA) ) / log(QA);}  // Improvement
  return optk;
}

// [[Rcpp::export(.calceo)]]
double calceo(DataFrame df, NumericVector coeff, bool yemp) {
  int y, row;
  double x;
  NumericVector col1, col2, rnd;

  col1 = df[0];
  col2 = df[1];
  rnd = runif(1);
  row = rnd[0] * df.nrows();
  x = gettherisk(col1[row], coeff);
  if (yemp == true) {
    y = col2[row];
  }
  else {
    NumericVector rndm = runif(1);
    double rdm = as<double>(rndm);
    y = (rdm < x ? 1 : 0);
  }
  return x - y;
}

double calceo2(DataFrame df, NumericVector coeff, NumericVector coeff2, double QS) {
  int y, row, s;
  double x, Qstar, xstar, rdm, x2;
  NumericVector col1, col2, rnd, rndm;

  col1 = df[0];                                       // Parsonnet score stored in the dataframe
  col2 = df[1];                                       // Surgical outcome stored in the dataframe
  rnd = runif(1);
  row = rnd[0] * df.nrows();
  s = col1[row];                                      // Step 1: s=Parsonnet score
  x = gettherisk(col1[row], coeff);                   // Step 2
  Qstar = QS;                                         // Step 3
  xstar = (Qstar*x) / (1 - x + Qstar * x);
  rndm = runif(1);
  rdm = as<double>(rndm);
  y = (rdm < xstar ? 1 : 0);                          // Step 4: y=Surgical outcome
  x2 = gettherisk(s, coeff2);                         // Step 5: x=True probability of death
  return x2 - y;
}

// [[Rcpp::export(.eocusum_arl_sim)]]
int eocusum_arl_sim(int r, double k, double h, DataFrame df, NumericVector coeff, bool yemp, int side) {
  int rl = 0;
  double z = 0;
  if (side == 1) {                                    // lower side (deterioration)
    double tn = 0;
    do{
      rl++;
      z = calceo(df, coeff, yemp);
      tn = fmin(0, tn + z + k);
    } while (-tn <= h);
    return rl;
  }
  else if (side == 2){                                // upper side (improvement)
    double qn = 0;
    do{
      rl++;
      z = calceo(df, coeff, yemp);
      qn = fmax(0, qn + z - k);
    } while (qn <= h);
    return rl;
  }
  return rl;
}


// [[Rcpp::export(.eocusum_arloc_sim)]]
int eocusum_arloc_sim(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int side) {
  int rl = 0;
  double z = 0;
  if (side == 1) {                                         // lower side (deterioration)
    double tn = 0;
    do{
      rl++;
      z = calceo2(df, coeff, coeff2, QS);
      tn = fmin(0, tn + z + k);
    } while (-tn <= h);
    return rl;
  }
  else if (side == 2){                                    // upper side (improvement)
    double qn = 0;
    do{
      rl++;
      z = calceo2(df, coeff, coeff2, QS);
      qn = fmax(0, qn + z - k);
    } while (qn <= h);
    return rl;
  }
  return rl;
}

// [[Rcpp::export(.eocusum_ad_sim)]]
int eocusum_ad_sim(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int side, int type, int m) {
  if (type == 1) {                                                  // conditional steady-state ARL (EO-CUSUM)
    if (side == 1) {
      return eocusum_ad_sim11(r, k, h, df, coeff, coeff2, QS, m); // lower side (deterioration)
    } else {
      return eocusum_ad_sim12(r, k, h, df, coeff, coeff2, QS, m); // upper side (improvement)
    }
  } else {                                                          // cyclical steady-state ARL (EO-CUSUM)
    if (side == 1) {
      return eocusum_ad_sim21(r, k, h, df, coeff, coeff2, QS, m); // lower side (deterioration)
    } else {
      return eocusum_ad_sim22(r, k, h, df, coeff, coeff2, QS, m); // upper side (improvement)
    }
  }
}

// conditional steady-state ARL (EO-CUSUM) -- m = #ic-observations with m >= 0
// lower side (deterioration)
int eocusum_ad_sim11(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m) {
  int success = 0, rl = 0;
  double z = 0;
  double tn = 0, Q = 1;
  while ( !success ) {
    rl = 0;
    tn = 0;
    Q = 1;
    do{
      rl++;
      if ( rl > m) Q = QS;
      z = calceo2(df, coeff, coeff2, Q);
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
int eocusum_ad_sim12(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS , int m) {
  int success = 0, rl = 0;
  double z = 0;
  double qn = 0, Q = 1;
  while ( !success ) {
    rl = 0;
    qn = 0;
    Q = 1;
    do{
      rl++;
      if ( rl > m) Q = QS;
      z = calceo2(df, coeff, coeff2, Q);
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
int eocusum_ad_sim21(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m) {
  int rl = 0;
  double z = 0;
  double tn = 0, Q = 1;
  do{
    rl++;
    if ( rl > m) Q = QS;
    z = calceo2(df, coeff, coeff2, Q);
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
int eocusum_ad_sim22(int r, double k, double h, DataFrame df, NumericVector coeff, NumericVector coeff2, double QS, int m) {
  int rl = 0;
  double z = 0;
  double qn = 0, Q = 1;
  do{
    rl++;
    if ( rl > m) Q = QS;
    z = calceo2(df, coeff, coeff2, Q);
    qn = fmax(0, qn + z - k);
    if ( rl <= m ) {
      if ( qn > h) {
        qn = 0;
      }
    }
  } while (qn <= h);
  rl += -m;
  return rl;
}
