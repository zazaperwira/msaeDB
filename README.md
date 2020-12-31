
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msaeDB

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/zazaperwira/msaeDB.svg?branch=master)](https://travis-ci.com/zazaperwira/msaeDB)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/zazaperwira/msaeDB?branch=master&svg=true)](https://ci.appveyor.com/project/zazaperwira/msaeDB)
<!-- badges: end -->

The goal of msaeDB is to implement Benchmarking Method for Multivariate
Small Area Estimation under Fay Herriot Model.  
Multivariate Small Area Estimation (MSAE) is a development of Univariat
Small Area Estimation that considering the correlation among response
variables and borrowing the strength from auxiliary variables
effectiveness of a domain sample size, the multivariat model in this
package is based on Multivariate model 1 proposed by Roberto Benavent
and Domigo Morales (2015) \<DOI: 10.1016/j.csda.2015.07.013.\>.
Benchmarking in Small Area Estimation is a modification of Small Area
Estimation model to guarantees that the aggreagate weighted mean of the
county predictors equals the corresponding weighted mean of survey
estimates. Difference Benchmarking is the simplest but widely used by
multiplying EBLUP estimator by the common adjustment factor (J.N.K Rao
and Isabela Molina, 2013).

## Installation

You can install the released version of msaeDB from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("msaeDB")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(msaeDB)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
