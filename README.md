
<!--  -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# vlad

<!-- [![Build Status](https://travis-ci.org/wittenberg/vlad.svg)](https://travis-ci.org/wittenberg/vlad) -->

<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/vlad)](http://cran.r-project.org/package=vlad) -->

<!-- [![Coverage Status](https://codecov.io/gh/wittenberg/vlad/graph/badge.svg)](https://codecov.io/github/wittenberg/vlad?branch=master) -->

<!-- [![Downloads](https://cranlogs.r-pkg.org/badges/vlad)](https://CRAN.R-project.org/package=vlad) -->

<!-- [![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/vlad?color=orange)](https://CRAN.R-project.org/package=vlad) -->

An R-package which contains functions to set up risk-adjusted quality
control charts in health care.

## Main features

  - Risk-adjusted CUSUM chart control limit calculations based on fast
    and accurate Markov chain approximations (Knoth, Wittenberg, and Gan
    [2019](#ref-Knoth.etal_2019))
  - Risk-adjusted CUSUM chart based on E-O (Wittenberg, Gan, and Knoth
    [2018](#ref-Wittenberg.etal_2018))
  - Risk-adjusted CUSUM chart based on log-likelihood ratio statistic
    (Steiner et al. [2000](#ref-Steiner.etal_2000))
  - Algorithms are implemented using Rcpp, RcppArmadillo and parallel
    computation

## Installation

You can install the released version of **vlad** from
[CRAN](https://CRAN.R-project.org/package=vlad) with:

``` r
install.packages("vlad")
```

You can install the development version from
[GitHub](https://github.com/wittenberg/vlad) with:

``` r
# install.packages("devtools")
devtools::install_github("wittenberg/vlad")
```

The package compiles C++ source code during installation, therefore you
need the appropriate compilers:

On *Windows* you need
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) available from
CRAN.

On *macOS* you need the Clang 6.x compiler and the GNU Fortran compiler
from [macOS tools](https://cran.r-project.org/bin/macosx/tools/). Having
installed the compilers, you need to open a terminal and start R via
‘PATH=/usr/local/clang6/bin:$PATH R’. You can then install the package
via *devtools::install\_github(“wittenberg/vlad”)*

## Examples

[Supplemental material: Risk-adjusted CUSUM charts under model
error](https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fsim.8104&file=SIM_8104-Supp-0001-Supplemental-material-SIM-18-0571.R)

<!-- Load libraries: -->

<!-- ```{r load libraries, message=FALSE} -->

<!-- library("vlad") -->

<!-- library("dplyr") -->

<!-- library("tidyr") -->

<!-- library("ggplot2") -->

<!-- ``` -->

<!-- Subset the dataset `cardiacsurgery` into Phase I (first two years) and Phase II (five years) and estimate a risk model based on `phaseI`. -->

<!-- ```{r estimate risk model} -->

<!-- data("cardiacsurgery", package = "spcadjust") -->

<!-- cardiacsurgery <- cardiacsurgery %>% rename(s = Parsonnet) %>% -->

<!--   mutate(y = ifelse(status == 1 & time <= 30, 1, 0), -->

<!--         phase = factor(ifelse(date < 2*365, "I", "II"))) -->

<!-- head(cardiacsurgery) -->

<!-- phaseI <- filter(cardiacsurgery, phase == "I") %>% select(s, y) -->

<!-- coeff <- round(coef(glm(y ~ s, data = phaseI, family = "binomial")), 3) -->

<!-- print(coeff) -->

<!-- ``` -->

<!-- ### Create VLADs for seven surgeons -->

<!-- By using the estimated risk model coefficients `coeff`, for each pair of Parsonnet score `s` and operation outcome values `y`, the difference between expected and observed outcome is calculated with the function `calceo()`. -->

<!-- Thereafter, differences are cummulated to create the VLAD. This is done for all seven surgeons of the `cardiacsurgery` dataset. Results are saved to the object `vlads7`. -->

<!-- ```{r vlads7} -->

<!-- vlads7 <- lapply(1:7, function(j){ -->

<!--   Si <- filter(cardiacsurgery, surgeon == j) -->

<!--   EO <- sapply(seq_along(Si$s), function(i) calceo(df = Si[i, c("s", "y")], coeff = coeff)) -->

<!--   select(Si, surgeon, phase) %>%  mutate(n = 1:length(EO), cEO = cumsum(EO)) -->

<!-- })  -->

<!-- ``` -->

<!-- Create Variable life-adjusted Displays for each surgeon from the object `vlads7`. -->

<!-- ```{r VLADS1-7, fig.align='center', fig.width=8, fig.height=10} -->

<!-- vlads7 %>%  -->

<!--   bind_rows() %>%   -->

<!--   gather(key = "Surgeon", value = value, c(-n, -surgeon, -phase)) %>% -->

<!--   ggplot(aes(x = n, y = value, colour = phase, group = Surgeon)) + -->

<!--     geom_hline(yintercept = 0, colour = "darkgreen", linetype = "dashed") + -->

<!--     geom_line(size = 1.1) + facet_wrap( ~ surgeon, ncol = 2, scales = "free") + -->

<!--     labs(x="Patient number n", y="CUSUM E-O") + theme_classic() + -->

<!--     scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) + -->

<!--     scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) -->

<!-- ``` -->

<!-- ### Create a VLAD for surgeon 2 -->

<!-- ```{r vladS2, fig.align='center'} -->

<!-- S2 <- filter(cardiacsurgery, surgeon == 2) %>% select(phase, s, y) -->

<!-- S2I <- subset(S2, c(phase == "I")) -->

<!-- S2II <- subset(S2, c(phase == "II")) -->

<!-- coeff <- coef(glm(y ~ s, data = S2I, family = "binomial")) -->

<!-- EO <- sapply(1:nrow(S2), function(i) calceo(df = S2[i, c("s", "y")], coeff = coeff)) -->

<!-- df1 <- select(S2, phase) %>% mutate(n = row_number(), cEO = cumsum(EO)) -->

<!-- df2 <- gather(df1, variable, value, c(-n, -phase)) -->

<!-- p1 <- ggplot(df2, aes(x = n, y = value, colour = phase)) + -->

<!--   geom_hline(yintercept = 0, linetype = "dashed") + geom_line() + geom_point() +  -->

<!--   labs(x = "Patient number", y = "CUSUM E-O") + theme_classic() + -->

<!--   scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) + -->

<!--   scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) -->

<!-- p1 -->

<!-- ``` -->

<!-- ### Compute thresholds of a risk-adjusted CUSUM chart for surgeon 2  -->

<!-- Upper and lower control limits of the risk-adjusted CUSUM chart based on log-likelihood ratio statistic can be computed with the function `racusum_arl_h_sim()`. The implemention uses parallel simulation and a multi-stage search procedure.   -->

<!-- ```{r} -->

<!-- # set a random number generator for parallel computations -->

<!-- RNGkind("L'Ecuyer-CMRG") -->

<!-- # number of simulation runs -->

<!-- m <- 10^4 -->

<!-- # assign cores -->

<!-- nc <- parallel::detectCores() -->

<!-- # verbose calculation  -->

<!-- UCL_sim <- racusum_crit_sim(L0 = 740, df = S2I[, c("s", "y")], coeff = coeff, m = m, RA = 2, nc = nc, verbose = TRUE) -->

<!-- # quite calculation -->

<!-- LCL_sim <- racusum_crit_sim(L0 = 740, df = S2I[, c("s", "y")], coeff = coeff, m = m, RA = 1/2, nc = nc, verbose = FALSE) -->

<!-- round(cbind(UCL_sim, LCL_sim), 3) -->

<!-- ``` -->

### Authors

Philipp Wittenberg and Sven Knoth

### License

GPL (\>= 2)

### References

<div id="refs" class="references hanging-indent">

<div id="ref-Knoth.etal_2019">

Knoth, S, P Wittenberg, and FF Gan. 2019. “Risk-Adjusted CUSUM Charts
Under Model Error.” *Statistics in Medicine* 38 (12): 2206–18.
<https://doi.org/10.1002/sim.8104>.

</div>

<div id="ref-Steiner.etal_2000">

Steiner, SH, RJ Cook, VT Farewell, and T Treasure. 2000. “Monitoring
Surgical Performance Using Risk-Adjusted Cumulative Sum Charts.”
*Biostatistics* 1 (4): 441–52.
<https://doi.org/10.1093/biostatistics/1.4.441>.

</div>

<div id="ref-Wittenberg.etal_2018">

Wittenberg, P, FF Gan, and S Knoth. 2018. “A Simple Signaling Rule for
Variable Life-Adjusted Display Derived from an Equivalent Risk-Adjusted
CUSUM Chart.” *Statistics in Medicine* 37 (16): 2455–73.
<https://doi.org/10.1002/sim.7647>.

</div>

</div>
