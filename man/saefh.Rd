% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saefh.R
\name{saefh}
\alias{saefh}
\title{EBLUPs under Univariate Fay Herriot Model}
\usage{
saefh(formula, vardir, samevar = FALSE, MAXITER = 100, PRECISION = 1e-04, data)
}
\arguments{
\item{formula}{List of formula that describe the fitted model}

\item{vardir}{Sampling variances of direct estimations,if it is included in data frame so it is the vector with the name of sampling variances.if it is not, it is a data frame of sampling variance in order : \code{var1, cov12,.,cov1r,var2,cov23,.,cov2r,.,cov(r-1)(r),var(r)}}

\item{samevar}{Whether the variances of the data are same or not. Logical input with default \code{FALSE}}

\item{MAXITER}{Maximum number of iteration in Fisher-scoring algorithm with default \code{100}}

\item{PRECISION}{Limit of Fisher-scoring convergence tolerance with default \code{1e-4}}

\item{data}{The data frame}
}
\value{
This function returns a list of the following objects:
   \item{SAE_Eblup}{A dataframe with the values of the EBLUPs estimators}
   \item{MSE_Eblup}{A dataframe with the values of estimated mean square errors of EBLUPs estimators}
   \item{fit}{A list containing the following objects:}
     \itemize{
       \item method : The fitting method (this function is using "REML")
       \item convergence : The convergence result of Fisher-scoring algorithm (Logical Value)
       \item iterations : The number of Fisher-Scoring algorithm iterations
       \item estcoef : A dataframe with the estimated model coefficient, standard error,t statistics, p-values of the significance of each coefficient
       \item refvar : A dataframe with estimated random effect variances
       \item informationFisher : A matrix of information fisher from Fisher-scoring algorithm
     }
}
\description{
This function produces EBLUPs, MSE of Univariate SAE
}
\examples{

##load dataset
data(datamsaeDB)

#Compute Fitted model for Y1, Y2, and Y3
#Y1 ~ X1 + X2
#Y2 ~ X2
#Y3 ~ X1

##Using parameter 'data'
formula = list(f1 = Y1~X1+X2,
               f2 = Y2~X2,
               f3 = Y3~X1)
vardir = c("v1","v12","v13","v2","v23","v3")
#Note : in real data for univariate SAE, if you does not have the values of covariances,
#       set covariancse as zero in the dataframe

saeFH <- saefh(formula, vardir, data=datamsaeDB)

##Do not use parameter 'data'
formula = list(f1 = datamsaeDB$Y1~datamsaeDB$X1+datamsaeDB$X2,
               f2 = datamsaeDB$Y2~datamsaeDB$X2,
               f3 = datamsaeDB$Y3~datamsaeDB$X1)
vardir = datamsaeDB[,c("v1","v12","v13","v2","v23","v3")]
#Note : in real data for univariate SAE, if you does not have the values of covariances,
#       set covariancse as zero in the dataframe
saeFH_d <- saefh(formula, vardir)

saeFH$SAE_Eblup       #to see EBLUP Estimators
saeFH$MSE_Eblup        #to see estimated MSE of EBLUP estimators


}
