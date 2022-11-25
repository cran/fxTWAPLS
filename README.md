
<!-- README.md is generated from README.Rmd. Please edit that file -->

## fxTWAPLS: An Improved Version of WA-PLS <img src="https://raw.githubusercontent.com/special-uor/fxTWAPLS/master/inst/images/logo.png" alt="logo" align="right" height=200px/>

<!-- badges: start -->
<!-- [![](https://img.shields.io/github/languages/code-size/special-uor/fxTWAPLS.svg)](https://github.com/special-uor/fxTWAPLS) -->

[![](https://img.shields.io/badge/devel%20version-0.1.2-yellow.svg)](https://github.com/special-uor/fxTWAPLS)
[![](https://www.r-pkg.org/badges/version/fxTWAPLS?color=black)](https://cran.r-project.org/package=fxTWAPLS)
[![R build
status](https://github.com/special-uor/fxTWAPLS/workflows/R-CMD-check/badge.svg)](https://github.com/special-uor/fxTWAPLS/actions)
[![](https://img.shields.io/badge/doi-10.1098/rspa.2020.0346-blue.svg)](https://doi.org/10.1098/rspa.2020.0346)
<!-- [![](https://app.codecov.io/gh/special-uor/fxTWAPLS/branch/master/graph/badge.svg?token=Q6SYL7AOGR)](https://app.codecov.io/gh/special-uor/fxTWAPLS) -->
<!-- [![R build status](https://github.com/special-uor/fxTWAPLS/workflows/R-CMD-check/badge.svg)](https://github.com/special-uor/fxTWAPLS/actions) -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/fxTWAPLS)](https://CRAN.R-project.org/package=fxTWAPLS) -->
<!-- badges: end -->

## Overview

The goal of this package is to provide an improved version of WA-PLS by
including the tolerances of taxa and the frequency of the sampled
climate variable. This package also provides a way of leave-out
cross-validation that removes both the test site and sites that are both
geographically close and climatically close for each cycle, to avoid the
risk of pseudo-replication.

## Installation

<!-- ### Create a Personal Access Token (PAT) for Github -->
<!-- This is needed to install packages from private repositories. Once configured, -->
<!-- there is no need to configure it again. -->

You can install the released version of fxTWAPLS from
[CRAN](https://cran.r-project.org/package=fxTWAPLS) with:

``` r
install.packages("fxTWAPLS")
```

And the development version from
[GitHub](https://github.com/special-uor/fxTWAPLS/) with:
<!-- You can install the development version from [GitHub](https://github.com/) with: -->

``` r
install.packages("remotes")
remotes::install_github("special-uor/fxTWAPLS", "dev")
```

## Publications

-   ***Latest:*** Liu, M., Shen, Y., González-Sampériz, P., Gil-Romera,
    G., ter Braak, C. J. F., Prentice, I. C., and Harrison, S. P.:
    Holocene climates of the Iberian Peninsula: pollen-based
    reconstructions of changes in the west-east gradient of temperature
    and moisture, Clim. Past Discuss. \[preprint\],
    <https://doi.org/10.5194/cp-2021-174>, in review, 2021.-
    [`fxTWAPLS v0.1.0`](https://github.com/special-uor/fxTWAPLS/releases/tag/v0.1.0/)

    ``` r
    install.packages("remotes")
    remotes::install_github("special-uor/fxTWAPLS@v0.1.0")
    ```

-   Liu Mengmeng, Prentice Iain Colin, ter Braak Cajo J. F., Harrison
    Sandy P.. 2020 An improved statistical approach for reconstructing
    past climates from biotic assemblages. *Proc. R. Soc. A.*
    **476**: 20200346. <https://doi.org/10.1098/rspa.2020.0346> -
    [`fxTWAPLS v0.0.2`](https://github.com/special-uor/fxTWAPLS/releases/tag/v0.0.2/)

    ``` r
    install.packages("remotes")
    remotes::install_github("special-uor/fxTWAPLS@v0.0.2")
    ```

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->

## Notes

The following functions can be executed in parallel:

-   [`cv.pr.w`](https://special-uor.github.io/fxTWAPLS/reference/cv.pr.w.html)
-   [`cv.w`](https://special-uor.github.io/fxTWAPLS/reference/cv.w.html)
-   [`get_distance`](https://special-uor.github.io/fxTWAPLS/reference/get_distance.html)
-   [`get_pseudo`](https://special-uor.github.io/fxTWAPLS/reference/get_pseudo.html)
-   [`sse.sample`](https://special-uor.github.io/fxTWAPLS/reference/sse.sample.html)

To do so, include the `cpus` parameter. For example:

``` r
cv_pr_tf_Tmin2 <- fxTWAPLS::cv.pr.w(
  taxa,
  modern_pollen$Tmin,
  nPLS = 5,
  fxTWAPLS::TWAPLS.w2,
  fxTWAPLS::TWAPLS.predict.w,
  pseudo_Tmin,
  usefx = TRUE,
  fx_method = "pspline",
  bin = 0.02,
  cpus = 2
)
```

Optionally, a progress bar can be displayed for long computations. Just
“pipe” the function call to `fxTWAPLS::pb()`.

``` r
`%>%` <- magrittr::`%>%`
cv_pr_tf_Tmin2 <- fxTWAPLS::cv.pr.w(
  taxa,
  modern_pollen$Tmin,
  nPLS = 5,
  fxTWAPLS::TWAPLS.w2,
  fxTWAPLS::TWAPLS.predict.w,
  pseudo_Tmin,
  usefx = TRUE,
  fx_method = "pspline",
  bin = 0.02,
  cpus = 2
) %>%
  fxTWAPLS::pb()
```

Alternatively, if you are not familiar with the “pipe” operator, you can
run the following code:

``` r
cv_pr_tf_Tmin2 <- fxTWAPLS::pb(
  fxTWAPLS::cv.pr.w(
    taxa,
    modern_pollen$Tmin,
    nPLS = 5,
    fxTWAPLS::TWAPLS.w2,
    fxTWAPLS::TWAPLS.predict.w,
    pseudo_Tmin,
    usefx = TRUE,
    fx_method = "pspline",
    bin = 0.02,
    cpus = 2
  )
)
```

## Example

#### Training

``` r
# Load modern data
modern_pollen <- read.csv("/path/to/modern_pollen.csv")

# Extract modern pollen taxa
taxaColMin <- which(colnames(modern_pollen) == "taxa0")
taxaColMax <- which(colnames(modern_pollen) == "taxaN")
taxa <- modern_pollen[, taxaColMin:taxaColMax]

# Set the binwidth to get the sampling frequency of the climate (fx),
# the fit is almost insenitive to binwidth when choosing pspline method.
bin <- 0.02

# Use fxTWAPLSv2 to train
fit_tf_Tmin2 <- fxTWAPLS::TWAPLS.w2(
  taxa,
  modern_pollen$Tmin,
  nPLS = 5,
  usefx = TRUE,
  fx_method = "pspline",
  bin = bin
)
```

#### Cross validation

``` r
# Set CPUS to run in parallel
CPUS <- 6

# Import pipe operator to use with the progress bar
`%>%` <- magrittr::`%>%`

# Get the location information of each sample
point <- modern_pollen[, c("Long", "Lat")]

# Get the distance between each point
dist <- fxTWAPLS::get_distance(point, cpus = CPUS)

# Get the pseudo sites (which are both geographically close and climatically
# close to the test site) which should be removed in cross validation
pseudo_Tmin <- fxTWAPLS::get_pseudo(
  dist,
  modern_pollen$Tmin,
  cpus = CPUS
)

# Leave-out cross validation
cv_pr_tf_Tmin2 <- fxTWAPLS::cv.pr.w(
  taxa,
  modern_pollen$Tmin,
  nPLS = 5,
  fxTWAPLS::TWAPLS.w2,
  fxTWAPLS::TWAPLS.predict.w,
  pseudo_Tmin,
  usefx = TRUE,
  fx_method = "pspline",
  bin = bin,
  cpus = CPUS,
  test_mode = FALSE
) %>%
  fxTWAPLS::pb()

# Random t test to the cross validation result
rand_pr_tf_Tmin2 <-
  fxTWAPLS::rand.t.test.w(cv_pr_tf_Tmin2, n.perm = 999)
```

#### Reconstruction

``` r
# Load fossil data
Holocene <- read.csv("/path/to/Holocene.csv")

# Extract fossil pollen taxa
taxaColMin <- which(colnames(Holocene) == "taxa0")
taxaColMax <- which(colnames(Holocene) == "taxaN")
core <- Holocene[, taxaColMin:taxaColMax]

# Choose nsig (the last significant number of components) based on the p-value
nsig <- 3

# Predict
fossil_tf_Tmin2 <- fxTWAPLS::TWAPLS.predict.w(fit_tf_Tmin2, core)

# Get the sample specific errors
sse_tf_Tmin2 <- fxTWAPLS::sse.sample(
  modern_taxa = taxa,
  modern_climate = modern_pollen$Tmin,
  fossil_taxa = core,
  trainfun = fxTWAPLS::TWAPLS.w2,
  predictfun = fxTWAPLS::TWAPLS.predict.w,
  nboot = nboot,
  nPLS = 5,
  nsig = nsig,
  usefx = TRUE,
  fx_method = "pspline",
  bin = bin,
  cpus = CPUS
) %>%
  fxTWAPLS::pb()
# Output
recon_result <-
  cbind.data.frame(
    recon_Tmin = fossil_tf_Tmin2[["fit"]][, nsig],
    sse_recon_Tmin = sse_tf_Tmin2
  )
```
