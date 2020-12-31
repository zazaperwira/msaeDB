#' msaeDB : Multivariate Small Area Estimation with Difference Benchmarking
#'
#' This Package is implementing Benchmarking Method for Multivariate Small Area Estimation under Fay-Herriot Model. Multivariate Small Area Estimation (MSAE) is a development of Univariate Small Area Estimation that considering the correlation among response variables and borrowing the strength from auxiliary variables effectiveness of a domain sample size, the multivariate model in this package is based on Multivariate model 1 proposed by Roberto Benavent and Domigo Morales (2016). Benchmarking in Small Area Estimation is a modification of Small Area Estimation model to guarantees that the aggregate weighted mean of the county predictors equals the corresponding weighted mean of survey estimates. Difference Benchmarking is the simplest but widely used  by multiplying EBLUP estimator by the common adjustment factor (J.N.K Rao and Isabela Molina, 2013).
#'
#' @section Author(s):
#' Zaza Yuda Perwira, Azka Ubaidillah
#'
#' \strong{Maintainer}: Zaza Yuda Perwira \email{221710086@@stis.ac.id}
#'
#' @section Functions:
#' \describe{
#'   \item{\code{link{msaedb}}}{Produces EBLUPs, MSE, and Aggregarion of Multivariat SAE with Difference Benchmarking}
#'   \item{\code{link{usaedb}}}{Produces EBLUPs, MSE, and Aggregarion of Univariate SAE with Difference Benchmarking}
#'   \item{\code{link{msaefh}}}{Produces EBLUPs and MSE of Multivariat SAE}
#'   \item{\code{link{usaefh}}}{Produces EBLUPs and MSE of Univariate SAE}
#'}
#' @section Reference:
#' \itemize{
#'   \item{Benavent, Roberto & Morales, Domigo. (2016). Multivariate Fay-Herriot models for small area estimation. Computational Statistics and Data Analysis 94 2016 372-390. DOI: 10.1016/j.csda.2015.07.013.}
#'   \item{Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New York: John Wiley and Sons, Inc.}
#'   \item{Steorts, Rebecca & Ghosh, Malay. (2013). On estimation of mean square Errors of Benchmarked Empirical Bayes Estimators. Article in Statistica Sinica April 2013. DOI: 10.5705/ss.2012.053.}
#'   \item{Permatasari, Novia. (2020). Pembangunan paket R pada model Fay Herriot multivariat untuk pendugaan area kecil (Bachelor Thesis). Jakarta: Polytechnic Statistics of STIS}
#'   \item{Ubaidillah, Azka et al. (2019). Multivariate Fay-Herriot models for small area estimation with application to household consumption per capita expenditure in Indonesia. Journal of Applied Statistics. 46:15. 2845-2861. DOI: 10.1080/02664763.2019.1615420.}
#' }
#'
#'
#' @docType package
#' @name msaeDB
#'
#' @import magic
#' @import MASS
#' @import stats

NULL
