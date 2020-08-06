#include "racusum_discretebeta_sim.h"

// [[Rcpp::export(.racusum_discretebeta_arl_sim)]]
int racusum_discretebeta_arl_sim(int r, double shape1, double shape2, NumericVector coeff, double h, double RA, int rs, double RQ) {
  double qn = 0, wt = 0, pt, pistar, logitp, s, logRA = log(RA);
  int rl = 0, y;
  do{                                             // In-/Out-of-Control ARL
    rl++;
    s = floor( R::rbeta(shape1, shape2) * rs );   // Sample and round s from beta distribution
    logitp = coeff[0] + s * coeff[1];             // Predicted probability
    pt = exp(logitp) / (1 + exp(logitp));         // True probability of death
    pistar = (RQ * pt) / (1 - pt + RQ * pt);      // Change based on Odds ratio change RQ
    y = (R::runif(0, 1) < pistar ? 1 : 0);        // Surgical outcome
    wt = -log(1 - pt + RA * pt) + y * logRA;      // LLR
    qn = fmax(0, qn + wt);
  } while (qn <= h);
  return rl;
}
