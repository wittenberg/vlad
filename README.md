
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vlad

<!-- [![Build Status](https://travis-ci.org/wittenberg/vlad.svg)](https://travis-ci.org/wittenberg/vlad)  -->

<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/vlad)](http://cran.r-project.org/package=vlad)  -->

<!-- [![Coverage Status](https://codecov.io/gh/wittenberg/vlad/graph/badge.svg)](https://codecov.io/github/wittenberg/vlad?branch=master)  -->

<!-- [![Downloads](https://cranlogs.r-pkg.org/badges/vlad)](https://CRAN.R-project.org/package=vlad) -->

<!-- [![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/vlad?color=orange)](https://CRAN.R-project.org/package=vlad) -->

An R-package which contains functions to set up risk-adjusted quality
control charts in health care.

## Main features

  - Risk-adjusted CUSUM chart based on E-O by Wittenberg et al. (2018)
  - Risk-adjusted CUSUM chart based on log-likelihood ratio statistic by
    Steiner et al. (2000)
  - Algorithms are implemented using Rcpp and RcppArmadillo
  - High performance with parallel computation

## Installation

You can install the released version of **vlad** from
[CRAN](https://cran.r-project.org/web/packages/vlad/index.html) with:

``` r
install.packages("vlad")
```

And the development version from
[GitHub](https://github.com/wittenberg/vlad) with:

``` r
# install.packages("devtools")
devtools::install_github("wittenberg/vlad")
```

## Example

Load libraries:

``` r
library("vlad")
library("dplyr")
library("tidyr")
library("ggplot2")
```

Subset the dataset `cardiacsurgery` into Phase I (first two years) and
Phase II (five years) and estimate a risk model based on `phaseI`.

``` r
data("cardiacsurgery", package = "spcadjust")
cardiacsurgery <- cardiacsurgery %>% rename(s = Parsonnet) %>%
  mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
        phase = factor(ifelse(date < 2*365, "I", "II")))
head(cardiacsurgery)
#>   date time status  s surgeon y phase
#> 1    1   90      0 15       7 0     I
#> 2    1   90      0  9       3 0     I
#> 3    2   90      0  2       5 0     I
#> 4    3   90      0  8       7 0     I
#> 5    3   90      0  7       1 0     I
#> 6    3   90      0 40       1 0     I
phaseI <- filter(cardiacsurgery, phase == "I") %>% select(s, y)
coeff <- round(coef(glm(y ~ s, data = phaseI, family = "binomial")), 3)
print(coeff)
#> (Intercept)           s 
#>       -3.79        0.08
```

### Create VLADs for seven surgeons

By using the estimated risk model coefficients `coeff`, for each pair of
Parsonnet score `s` and operation outcome values `y`, the difference
between expected and observed outcome is calculated with the function
`calceo()`. Thereafter, differences are cummulated to create the VLAD.
This is done for all seven surgeons of the `cardiacsurgery` dataset.
Results are saved to the object `vlads7`.

``` r
vlads7 <- lapply(1:7, function(j){
  Si <- filter(cardiacsurgery, surgeon == j)
  EO <- sapply(seq_along(Si$s), function(i) calceo(df = Si[i, c("s", "y")], coeff = coeff))
  select(Si, surgeon, phase) %>%  mutate(n = 1:length(EO), cEO = cumsum(EO))
}) 
```

Create Variable life-adjusted Displays for each surgeon from the object
`vlads7`.

``` r
vlads7 %>% 
  bind_rows() %>%  
  gather(key = "Surgeon", value = value, c(-n, -surgeon, -phase)) %>%
  ggplot(aes(x = n, y = value, colour = phase, group = Surgeon)) +
    geom_hline(yintercept = 0, colour = "darkgreen", linetype = "dashed") +
    geom_line(size = 1.1) + facet_wrap( ~ surgeon, ncol = 2, scales = "free") +
    labs(x="Patient number n", y="CUSUM E-O") + theme_classic() +
    scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL))
```

<img src="man/figures/README-VLADS1-7-1.png" width="100%" style="display: block; margin: auto;" />

### Create a VLAD for surgeon 2

``` r
S2 <- filter(cardiacsurgery, surgeon == 2) %>% select(phase, s, y)
S2I <- subset(S2, c(phase == "I"))
S2II <- subset(S2, c(phase == "II"))
coeff <- coef(glm(y ~ s, data = S2I, family = "binomial"))
EO <- sapply(1:nrow(S2), function(i) calceo(df = S2[i, c("s", "y")], coeff = coeff))

df1 <- select(S2, phase) %>% mutate(n = row_number(), cEO = cumsum(EO))
df2 <- gather(df1, variable, value, c(-n, -phase))

p1 <- ggplot(df2, aes(x = n, y = value, colour = phase)) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_line() + geom_point() + 
  labs(x = "Patient number", y = "CUSUM E-O") + theme_classic() +
  scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL))
p1
```

<img src="man/figures/README-vladS2-1.png" width="100%" style="display: block; margin: auto;" />

### Compute thresholds of a risk-adjusted CUSUM chart for surgeon 2

Upper and lower control limits of the risk-adjusted CUSUM chart based on
log-likelihood ratio statistic can be computed with the function
`racusum_arl_h_sim()`. The implemention uses parallel simulation and a
multi-stage search procedure.

``` r
# set a random number generator for parallel computations
RNGkind("L'Ecuyer-CMRG")
# number of simulation runs
m <- 10^4
# assign cores
nc <- parallel::detectCores()
# verbose calculation 
UCL_sim <- racusum_crit_sim(L0 = 740, df = S2I[, c("s", "y")], coeff = coeff, m = m, RA = 2, nc = nc, verbose = TRUE)
#> h = 1    ARL = 73.748 
#> h = 2    ARL = 371.8713 
#> h = 3    ARL = 1272.2543 
#> h = 2.9  ARL = 1139.8029 
#> h = 2.8  ARL = 1014.7675 
#> h = 2.7  ARL = 894.4185 
#> h = 2.6  ARL = 794.149 
#> h = 2.5  ARL = 699.988 
#> h = 2.51     ARL = 711.0032 
#> h = 2.52     ARL = 718.1516 
#> h = 2.53     ARL = 728.3317 
#> h = 2.54     ARL = 737.2428 
#> h = 2.55     ARL = 747.9701 
#> h = 2.549    ARL = 747.2548 
#> h = 2.548    ARL = 746.1113 
#> h = 2.547    ARL = 744.8156 
#> h = 2.546    ARL = 744.311 
#> h = 2.545    ARL = 742.484 
#> h = 2.544    ARL = 741.7208 
#> h = 2.543    ARL = 740.6374 
#> h = 2.542    ARL = 738.9439 
#> h = 2.5421   ARL = 739.1099 
#> h = 2.5422   ARL = 739.4318 
#> h = 2.5423   ARL = 739.7967 
#> h = 2.5424   ARL = 739.7967 
#> h = 2.5425   ARL = 739.953 
#> h = 2.5426   ARL = 739.953 
#> h = 2.5427   ARL = 740.1047
# quite calculation
LCL_sim <- racusum_crit_sim(L0 = 740, df = S2I[, c("s", "y")], coeff = coeff, m = m, RA = 1/2, nc = nc, verbose = FALSE)
round(cbind(UCL_sim, LCL_sim), 3)
#>      UCL_sim LCL_sim
#> [1,]   2.543   2.278
```

### References

Wittenberg et al. (2018). [A simple signaling rule for variable
life‐adjusted display derived from an equivalent risk‐adjusted CUSUM
chart](https://doi.org/10.1002/sim.7647)

Steiner et al. (2000). [Monitoring surgical performance using
risk-adjusted cumulative sum
charts](https://doi.org/10.1093/biostatistics/1.4.441)

### Authors

Philipp Wittenberg and Sven Knoth

### License

GPL (\>= 2)
