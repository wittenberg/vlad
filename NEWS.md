# vlad 0.2.1
## Changes
* Use checkmate for error messages
* Calculate risk of death ahead of simulations to improve computation time
* Include various functions for distributions to compute the patient mix

## Improvements
* Improve simulation algorithms
* Improve testing

# vlad 0.2.0
## Changes
* Renamed functions:

original name         | new name
:---------------------|:-------------------
`cusum_arl_h_sim()`   | `bcusum_crit_sim()` 
`cusum_arl_sim()`     | `bcusum_arl_sim()`
`eocusum_arl_h_sim()` | `eocusum_crit_sim()`
`eocusum_adoc_sim()`  | `eocusum_ad_sim()`
`racusum_arl_h_sim()` | `racusum_crit_sim()`
`racusum_adoc_sim()`  | `racusum_ad_sim()`

## Improvements
* Include Markov Chain approximation for risk-adjusted CUSUM chart based on log-likelihood ratio statistic in function racusum_arl_mc() and racusum_crit_mc()
* Add option of calculation of optimal k based on mean outcome in function optimal_k()
* Add function for calculation of CUSUM scores for log-likelihood ratio CUSUM llr_cusum_scores()
and E-O CUSUM eo_cusum_scores()
* Extend unit tests for functions

# vlad 0.1.0
## Improvements
* Preparation for CRAN release
