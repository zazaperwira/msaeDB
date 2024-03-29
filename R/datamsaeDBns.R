#' @title Sample Data for Multivariate Small Area Estimation with Difference Benchmarking with clustering
#' @description Dataset to simulate difference benchmarking of Multivariate Fay Herriot model
#' This data is generated base on multivariate Fay Herriot model by these following steps:
#'   \enumerate{
#'     \item Generate explanatory variables \code{X1 and X2}. Take \eqn{\mu_{X1}}{\muX1} = \eqn{\mu_{X1}}{\muX1} = 10, \eqn{\sigma_{X11}}{\sigmaX11}=1, \eqn{\sigma_{X2}}{\sigmaX22}=2, and \eqn{\rho_{x}}{\rhox}= 1/2.
#'     \cr Sampling error \code{e} is generated with the following \eqn{\sigma_{e11}}{\sigmae11} = 0.15,  \eqn{\sigma_{e22}}{\sigmae22} = 0.25, \eqn{\sigma_{e33}}{\sigmae33} = 0.35, and \eqn{\rho_{e}}{\rhoe} = 1/2.
#'     \cr For random effect \code{u}, we set \eqn{\sigma_{u11}}{\sigmau11}= 0.2, \eqn{\sigma_{u22}}{\sigmau22}= 0.6, and \eqn{\sigma_{u33}}{\sigmau33}= 1.8.
#'     \cr For the weight we generate \code{w1 w2 w3} by set the \code{w1 ~ U(25,30)} , \code{w2 ~ U(25,30)}, \code{w3 ~ U(25,30)}
#'     \cr Calculate direct estimation \code{Y1 Y2 Y3} where \eqn{Y_{i}}{Yi} = \eqn{X * \beta + u_{i} + e_{i}}{X\beta+ui+ei}
#'     \item Then combine the direct estimations \code{Y1 Y2 Y3}, explanatory variables \code{X1 X2}, weights \code{w1 w2 w3}, and sampling varians covarians \code{v1 v12 v13 v2 v23 v3} in a data frame then named as datamsaeDB
#'   }
#'
#' @format A data frame with 30 rows and 18 variables:
#' \describe{
#'   \item{clY1}{cluster information of Y1}
#'   \item{clY2}{cluster information of Y2}
#'   \item{clY3}{cluster information of Y3}
#'   \item{nonsample}{A column with logical values, \code{TRUE} if the area is non-sampled}
#'   \item{Y1}{Direct Estimation of Y1}
#'   \item{Y2}{Direct Estimation of Y2}
#'   \item{Y3}{Direct Estimation of Y3}
#'   \item{X1}{Auxiliary variable of X1}
#'   \item{X2}{Auxiliary variable of X2}
#'   \item{w1}{Known proportion of units in small areas of Y1}
#'   \item{w2}{Known proportion of units in small areas of Y2}
#'   \item{w3}{Known proportion of units in small areas of Y3}
#'   \item{v1}{Sampling Variance of Y1}
#'   \item{v12}{Sampling Covariance of Y1 and Y2}
#'   \item{v13}{Sampling Covariance of Y1 and Y3}
#'   \item{v2}{Sampling Variance of Y2}
#'   \item{v23}{Sampling Covariance of Y2 and Y3}
#'   \item{v3}{Sampling Variance of Y3}
#' }
#'
"datamsaeDBns"
