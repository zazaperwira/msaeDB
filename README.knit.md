---
output: github_document
---


<!-- README.md is generated from README.Rmd. Please edit that file -->



# msaeDB

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/zazaperwira/msaeDB.svg?branch=master)](https://travis-ci.com/zazaperwira/msaeDB)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/zazaperwira/msaeDB?branch=master&svg=true)](https://ci.appveyor.com/project/zazaperwira/msaeDB)
[![Codecov test coverage](https://codecov.io/gh/zazaperwira/msaeDB/branch/master/graph/badge.svg)](https://codecov.io/gh/zazaperwira/msaeDB?branch=master)
<!-- badges: end -->

The goal of msaeDB is to implement Benchmarking Method for Multivariate Small Area Estimation under Fay Herriot Model.  
    Multivariate Small Area Estimation (MSAE) is a development of Univariat Small Area Estimation that
    considering the correlation among response variables and borrowing the strength from auxiliary variables 
    effectiveness of a domain sample size, the multivariat model in this package is based on Multivariate
    model 1 proposed by Roberto Benavent and Domigo Morales (2015) <DOI: 10.1016/j.csda.2015.07.013.>.
    Benchmarking in Small Area Estimation is a modification of Small Area Estimation model to guarantees that the 
    aggreagate weighted mean of the county predictors equals the corresponding weighted mean of survey estimates.
    Difference Benchmarking is the simplest but widely used  by multiplying EBLUP estimator by the common adjustment 
    factors (J.N.K Rao and Isabela Molina, 2013).

## Authors

Zaza Yuda Perwira, Azka Ubaidillah

## Maintainer 
Zaza Yuda Perwira <221710086@stis.ac.id>

## Functions
* `msaedb()` Produces EBLUPs, MSE, and Aggregarion of Multivariat SAE with Difference Benchmarking
* `usaedb()` Produces EBLUPs, MSE, and Aggregarion of Univariate SAE with Difference Benchmarking
* `msaefh()` Produces EBLUPs and MSE of Multivariat SAE
* `usaefh()` Produces EBLUPs and MSE of Univariate SAE

