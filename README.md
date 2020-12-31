
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msaeDB

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/zazaperwira/msaeDB.svg?branch=master)](https://travis-ci.com/zazaperwira/msaeDB)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/zazaperwira/msaeDB?branch=master&svg=true)](https://ci.appveyor.com/project/zazaperwira/msaeDB)
[![Codecov test
coverage](https://codecov.io/gh/zazaperwira/msaeDB/branch/master/graph/badge.svg)](https://codecov.io/gh/zazaperwira/msaeDB?branch=master)
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
multiplying EBLUP estimator by the common adjustment factors (J.N.K Rao
and Isabela Molina, 2013).

## Authors

Zaza Yuda Perwira, Azka Ubaidillah

## Maintainer

Zaza Yuda Perwira <221710086@stis.ac.id>

## Functions

  - `msaedb()` Produces EBLUPs, MSE, and Aggregarion of Multivariat SAE
    with Difference Benchmarking
  - `usaedb()` Produces EBLUPs, MSE, and Aggregarion of Univariate SAE
    with Difference Benchmarking
  - `msaefh()` Produces EBLUPs and MSE of Multivariat SAE
  - `usaefh()` Produces EBLUPs and MSE of Univariate SAE

## References

  - Benavent, Roberto & Morales, Domigo. (2016). Multivariate
    Fay-Herriot models for small area estimation. Computational
    Statistics and Data Analysis 94 2016 372-390. \<DOI:
    10.1016/j.csda.2015.07.013.\>
  - Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New
    York: John Wiley and Sons, Inc.
  - Steorts, Rebecca & Ghosh, Malay. (2013). On estimation of mean
    square Errors of Benchmarked Empirical Bayes Estimators. Article in
    Statistica Sinica April 2013. \<DOI: 10.5705/ss.2012.053.\>
  - Ubaidillah, Azka et al. (2019). Multivariate Fay-Herriot models for
    small area estimation with application to household consumption per
    capita expenditure in Indonesia. Journal of Applied Statistics.
    46:15. 2845-2861. \<DOI: 10.1080/02664763.2019.1615420\>
  - Permatasari, Novia. (2020). Pembangunan paket R pada model Fay
    Herriot multivariat untuk pendugaan area kecil (Bachelor Thesis).
    Jakarta: Polytechnic Statistics of STIS
