// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/vlad.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// optimal_k
double optimal_k(DataFrame pmix, double RA, bool yemp);
RcppExport SEXP _vlad_optimal_k(SEXP pmixSEXP, SEXP RASEXP, SEXP yempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type pmix(pmixSEXP);
    Rcpp::traits::input_parameter< double >::type RA(RASEXP);
    Rcpp::traits::input_parameter< bool >::type yemp(yempSEXP);
    rcpp_result_gen = Rcpp::wrap(optimal_k(pmix, RA, yemp));
    return rcpp_result_gen;
END_RCPP
}
// eocusum_arl_sim
int eocusum_arl_sim(int r, DataFrame pmix, double k, double h, double RQ, bool yemp, int side);
RcppExport SEXP _vlad_eocusum_arl_sim(SEXP rSEXP, SEXP pmixSEXP, SEXP kSEXP, SEXP hSEXP, SEXP RQSEXP, SEXP yempSEXP, SEXP sideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type pmix(pmixSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type RQ(RQSEXP);
    Rcpp::traits::input_parameter< bool >::type yemp(yempSEXP);
    Rcpp::traits::input_parameter< int >::type side(sideSEXP);
    rcpp_result_gen = Rcpp::wrap(eocusum_arl_sim(r, pmix, k, h, RQ, yemp, side));
    return rcpp_result_gen;
END_RCPP
}
// eocusum_ad_sim
int eocusum_ad_sim(int r, DataFrame pmix, double k, double h, double RQ, int side, int type, int m);
RcppExport SEXP _vlad_eocusum_ad_sim(SEXP rSEXP, SEXP pmixSEXP, SEXP kSEXP, SEXP hSEXP, SEXP RQSEXP, SEXP sideSEXP, SEXP typeSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type pmix(pmixSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type RQ(RQSEXP);
    Rcpp::traits::input_parameter< int >::type side(sideSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(eocusum_ad_sim(r, pmix, k, h, RQ, side, type, m));
    return rcpp_result_gen;
END_RCPP
}
// racusum_arl_mc
double racusum_arl_mc(NumericMatrix pmix, double RA, double RQ, double h, double scaling, int rounding, int method);
RcppExport SEXP _vlad_racusum_arl_mc(SEXP pmixSEXP, SEXP RASEXP, SEXP RQSEXP, SEXP hSEXP, SEXP scalingSEXP, SEXP roundingSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type pmix(pmixSEXP);
    Rcpp::traits::input_parameter< double >::type RA(RASEXP);
    Rcpp::traits::input_parameter< double >::type RQ(RQSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type scaling(scalingSEXP);
    Rcpp::traits::input_parameter< int >::type rounding(roundingSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(racusum_arl_mc(pmix, RA, RQ, h, scaling, rounding, method));
    return rcpp_result_gen;
END_RCPP
}
// racusum_crit_mc
double racusum_crit_mc(NumericMatrix pmix, double L0, double RA, double R, double scaling, int rounding, int method, int jmax, bool verbose);
RcppExport SEXP _vlad_racusum_crit_mc(SEXP pmixSEXP, SEXP L0SEXP, SEXP RASEXP, SEXP RSEXP, SEXP scalingSEXP, SEXP roundingSEXP, SEXP methodSEXP, SEXP jmaxSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type pmix(pmixSEXP);
    Rcpp::traits::input_parameter< double >::type L0(L0SEXP);
    Rcpp::traits::input_parameter< double >::type RA(RASEXP);
    Rcpp::traits::input_parameter< double >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type scaling(scalingSEXP);
    Rcpp::traits::input_parameter< int >::type rounding(roundingSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type jmax(jmaxSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(racusum_crit_mc(pmix, L0, RA, R, scaling, rounding, method, jmax, verbose));
    return rcpp_result_gen;
END_RCPP
}
// llr_score
double llr_score(DataFrame df, NumericVector coeff, double R0, double RA, bool yemp);
RcppExport SEXP _vlad_llr_score(SEXP dfSEXP, SEXP coeffSEXP, SEXP R0SEXP, SEXP RASEXP, SEXP yempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coeff(coeffSEXP);
    Rcpp::traits::input_parameter< double >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< double >::type RA(RASEXP);
    Rcpp::traits::input_parameter< bool >::type yemp(yempSEXP);
    rcpp_result_gen = Rcpp::wrap(llr_score(df, coeff, R0, RA, yemp));
    return rcpp_result_gen;
END_RCPP
}
// bcusum_arl_sim
int bcusum_arl_sim(int r, double h, DataFrame df, double R0, double RA);
RcppExport SEXP _vlad_bcusum_arl_sim(SEXP rSEXP, SEXP hSEXP, SEXP dfSEXP, SEXP R0SEXP, SEXP RASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< double >::type RA(RASEXP);
    rcpp_result_gen = Rcpp::wrap(bcusum_arl_sim(r, h, df, R0, RA));
    return rcpp_result_gen;
END_RCPP
}
// racusum_ad_sim
int racusum_ad_sim(int r, DataFrame pmix, double h, double RA, double RQ, int m, int type);
RcppExport SEXP _vlad_racusum_ad_sim(SEXP rSEXP, SEXP pmixSEXP, SEXP hSEXP, SEXP RASEXP, SEXP RQSEXP, SEXP mSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type pmix(pmixSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type RA(RASEXP);
    Rcpp::traits::input_parameter< double >::type RQ(RQSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(racusum_ad_sim(r, pmix, h, RA, RQ, m, type));
    return rcpp_result_gen;
END_RCPP
}
// racusum_arl_sim
int racusum_arl_sim(int r, DataFrame pmix, double h, double RA, double RQ, bool yemp);
RcppExport SEXP _vlad_racusum_arl_sim(SEXP rSEXP, SEXP pmixSEXP, SEXP hSEXP, SEXP RASEXP, SEXP RQSEXP, SEXP yempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type pmix(pmixSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type RA(RASEXP);
    Rcpp::traits::input_parameter< double >::type RQ(RQSEXP);
    Rcpp::traits::input_parameter< bool >::type yemp(yempSEXP);
    rcpp_result_gen = Rcpp::wrap(racusum_arl_sim(r, pmix, h, RA, RQ, yemp));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_vlad_optimal_k", (DL_FUNC) &_vlad_optimal_k, 3},
    {"_vlad_eocusum_arl_sim", (DL_FUNC) &_vlad_eocusum_arl_sim, 7},
    {"_vlad_eocusum_ad_sim", (DL_FUNC) &_vlad_eocusum_ad_sim, 8},
    {"_vlad_racusum_arl_mc", (DL_FUNC) &_vlad_racusum_arl_mc, 7},
    {"_vlad_racusum_crit_mc", (DL_FUNC) &_vlad_racusum_crit_mc, 9},
    {"_vlad_llr_score", (DL_FUNC) &_vlad_llr_score, 5},
    {"_vlad_bcusum_arl_sim", (DL_FUNC) &_vlad_bcusum_arl_sim, 5},
    {"_vlad_racusum_ad_sim", (DL_FUNC) &_vlad_racusum_ad_sim, 7},
    {"_vlad_racusum_arl_sim", (DL_FUNC) &_vlad_racusum_arl_sim, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_vlad(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
