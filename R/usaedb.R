#' @title EBLUPs under Univariate Fay Herriot Model with Difference Benchmarking
#' @description This function produces EBLUPs, MSE, and Aggregarion of Univariat SAE with Difference Benchmarking
#' @param formula list of formula that describe the fitted model
#' @param vardir  Sampling variance of direct estimation, if the data is included in data frame so it is the vector with the name of sampling variances. if it is not, it is a data frame of sampling variance in order : \code{var1, cov12,.,cov1r,var2,cov23,.,cov2r.,.cov(r-1)(r),var(r)}
#' @param weight  Vector of proportion of units in small area
#' @param samevar Whether the variances of the data is same or not. Logical input with default \code{FALSE}
#' @param MAXITER Maximum number of iteration in Fisher-scoring algorithm with default \code{100}
#' @param PRECISION Limit of fisher-scoring convergence tolerance with default \code{1e-4}
#' @param data The data frame
#'
#' @return This function returns a list of the following objects:
#'    \item{USAE_Eblup}{a dataframe with the values of the EBLUPs estimators}
#'    \item{MSE_Eblup}{a dataframe with the values of estimated mean square errors of EBLUPs estimators}
#'    \item{fit}{a list containing the following objects:}
#'      \itemize{
#'        \item method : The fitting method (this function is using "REML")
#'        \item convergence : the convergence result of fisher scoring algorithm (Logical Value)
#'        \item iterations : The number of Fisher-Scoring algorithm iterations
#'        \item estcoef : a dataframe with the estimated model coefficient, standard error,t statistics, p-values of the significance of each coefficient
#'        \item refvar : a dataframe with estimated random effect variance
#'        \item informationFisher : a matrix of information fisher from Fisher-scoring algorithm
#'      }
#'
#'    \item{difference_benchmarking}{a list containing the following objects:}
#'      \itemize{
#'        \item Estimation : a dataframe with the value of Benchmarked EBLUPs estimators
#'        \item Aggregation : The Aggregation of benchmarked EBLUPs estimations, EBLUPs Estimations and direct estimations
#'        \item MSE_DB : a dataframe with the values of estimated mean square errors of benchmarked EBLUPs estimators
#'        \item g.4a : first component of g4 in difference benchmarking MSE estimation formula
#'        \item g.4b : second component of g4 in difference benchmarking MSE estimation formula
#'      }
#'
#' @examples
#' ##load dataset
#' data(datamsaeDB)
#'
#' #Compute Fitted model for Y1, Y2, and Y3
#' #Y1 ~ X1 + X2
#' #Y2 ~ X2
#' #Y3 ~ X1
#'
#' ##Using parameter 'data'
#' formula = list(f1 = Y1~X1+X2,
#'                f2 = Y2~X2,
#'                f3 = Y3~X1)
#' vardir = c("v1","v12","v13","v2","v23","v3")
#' #Note : in real data for univariate SAE, if you does not have the valuse of covariances,
#' #       set covariancse as zero in the dataframe
#'
#' weight = c("w1","w2","w3")
#' usaeDB <- usaedb(formula, vardir, weight, data=datamsaeDB)
#'
#' ##Do not use parameter 'data'
#' formula = list(f1 = datamsaeDB$Y1~datamsaeDB$X1+datamsaeDB$X2,
#'                f2 = datamsaeDB$Y2~datamsaeDB$X2,
#'                f3 = datamsaeDB$Y3~datamsaeDB$X1)
#' vardir = datamsaeDB[,c("v1","v12","v13","v2","v23","v3")]
#' #Note : in real data for univariate SAE, if you does not have the valuse of covariances,
#' #       set covariancse as zero in the dataframe
#' weight = datamsaeDB[,c("w1","w2","w3")]
#' usaeDB <- usaedb(formula, vardir, weight = weight)
#'
#' usaeDB$USAE_Eblup       #to see EBLUP Estimators
#' usaeDB$MSE_Eblup        #to see estimated MSE of EBLUP estimators
#' usaeDB$difference_benchmarking$Estimation   #to see Benchmarked EBLUP Estimators
#' usaeDB$difference_benchmarking$MSE_DB       #to see estimated MSE of Benchmarked EBLUP Estimators
#' usaeDB$difference_benchmarking$Aggregation  #to see the aggregation of, benchmarking
#'
#' @export usaedb
#'
usaedb <- function (formula, vardir, weight, samevar = FALSE, MAXITER = 100, PRECISION = 1e-04,
                    data) {

  result = list(USAE_Eblup = NA, MSE_Eblup = NA,
                fit = list(method = NA, convergence = NA, iterations = NA,
                           estcoef = NA, refvar = NA, informationFisher = NA),
                difference_benchmarking = list(Estimation = NA,
                                               Aggregation = NA,
                                               MSE_DB = NA,
                                               g4.a = NA,
                                               g4.b = NA))
  r <- length(formula)
  RIn_function <- function(vardir, n,r){
    if (r == 1){
      RR <- diag(vardir)
    } else {
      it <- 0
      it2 <- 0
      RIn <- list()
      rmat2 <- matrix(0,n,n)
      for (i in 1:r){
        for (j in 1:r){
          it <- it +1
          if (i>j){
            RIn[[it+it2]] <- rmat2
            it <- it-1
            it2 <- it2+1
          }else { RIn[[it+it2]] <- diag(vardir[,it])
          }
        }
      }
      RIN <- list()
      for ( i in 1:r){
        if (i == 1){
          RIN[[i]] <- RIn[[i]]
          for (j in 1:(r-1)){
            RIN[[i]] <- cbind(RIN[[i]],RIn[[j+1]])
          }
        } else {
          RIN[[i]] <- RIn[[i*r-r+1]]
          for (j in 1:(r-1)){
            RIN[[i]] <- cbind(RIN[[i]],RIn[[(i*r-r+1)+j]])
          }
        }
      }
      RR <- do.call(rbind,RIN)
      RR <- (t(RR)+RR)*(matrix(1,n*r,n*r)-diag(0.5,n*r))
    }
    RIn<- return(RR)
  }
  if (!missing(data)) {
    formuladata <- formula
    for(i in 1:r) {formuladata[[i]] <- model.frame(formula[[i]], na.action = na.omit, data)}

    y.vec <- unlist(lapply(formuladata, function(x){x[1][1]}))


    x.matrix <- formula
    for(i in 1:r) {x.matrix[[i]] <- model.matrix(formula[[i]], na.action = na.omit, data)}
    x.matrix = Reduce(adiag,x.matrix)
    w.matrix = as.matrix(data[,weight])

    n = length(y.vec)/r

    if (any(is.na(data[, weight])))
      stop("Object weight contains NA values.")
    if (!all(weight %in% names(data)))
      stop("Object weight is not appropiate with data")
    if (length(weight) != r)
      stop("Length of weight is not appropiate, the length must be ",r)
    if (any(is.na(data[, vardir])))
      stop("Object vardir contains NA values.")
    if (!all(vardir %in% names(data)))
      stop("Object vardir is not appropiate with data")
    if (length(vardir) != sum(1:r))
      stop("Length of vardir is not appropiate, the length must be ", sum(1:r))

    RIn = RIn_function(data[, vardir],n,r)*diag(n*r)

  } else {

    formuladata <- formula
    for(i in 1:r) {formuladata[[i]] <- model.frame(formula[[i]], na.action = na.omit)}
    y.vec <- unlist(lapply(formuladata, function(x){x[1][1]}))

    x.matrix <- formula
    for(i in 1:r) {x.matrix[[i]] <- model.matrix(formula[[i]], na.action = na.omit)}
    x.matrix = Reduce(adiag,x.matrix)

    w.matrix = as.matrix(weight)

    n = length(y.vec)/r

    if ((dim(vardir)[2] != sum(1:r)) && (dim(vardir)[1] != n)) {
      stop("Object vardir is not appropiate with data, it must be ",n," x ",r," matrix")
    }
    if (any(is.na(vardir)))
      stop("Object vardir contains NA values.")
    if ((dim(weight)[2] != r) && (dim(vardir)[1] != n)) {
      stop("Object weight is not appropiate with data, it must be ",n," x ",r," matrix")
    }
    if (any(is.na(weight)))
      stop("Object weight contains NA values.")
    RIn = RIn_function(vardir,n, r)*diag(n*r)

  }

  for (i in 1:r) {
    if (attr(attributes(formuladata[[i]])$terms, "response") ==
        1)
      textformula = paste(formula[[i]][2], formula[[i]][1],
                          formula[[i]][3])
    else textformula = paste(formula[[i]][1], formula[[i]][2])
    if (length(na.action(formuladata[[i]])) > 0) {
      stop("Argument formula= ", textformula, " contains NA values.")
    }
  }

  varnames_Y <- lapply(formula, function(x) {x[[2]]})
  In = diag(n)
  Ir = diag(r)
  d.sigma <- lapply(formula, function(x){x=matrix(0,r,r)})
  for (i in 1:r) {d.sigma[[i]][i, i] = 1}
  d.SIGMA <- lapply(d.sigma, function(x){kronecker(x,In)})

  convergence = TRUE
  if (samevar) {
    Varu <- median(diag(RIn))
    k    <- 0
    diff <- rep(PRECISION + 1, r)
    while (any(diff > PRECISION) & (k < MAXITER)) {
      Varu1<- Varu
      G    <- kronecker(Varu, Ir)
      GIn  <- kronecker(G, In)
      SIGMA<-  (GIn + RIn)
      SIGMA_inv <- solve(SIGMA)
      Xt_Si<- t(SIGMA_inv %*% x.matrix)
      Q    <- solve(Xt_Si %*% x.matrix,tol = 1e-30)
      P    <- SIGMA_inv - t(Xt_Si) %*% Q %*% Xt_Si
      Py   <- P %*% y.vec
      s    <- (-0.5) * sum(diag(P)) + 0.5 * (t(Py) %*% Py)
      F    <- 0.5 * sum(diag(P %*% P))
      Varu <- Varu1 + solve(F) %*% s
      diff <- abs((Varu - Varu1)/Varu1)
      k <- k + 1
    }
    Varu = as.vector((rep(max(Varu,0), r)))
    names(Varu) = varnames_Y
    if (k >= MAXITER && diff >= PRECISION) {
      convergence = FALSE
    }

    GIn   <- kronecker(diag(Varu), In)
    SIGMA <-  (GIn + RIn)
    SIGMA_inv <- solve(SIGMA)
    Xt_Si <- t(SIGMA_inv %*% x.matrix)
    Q     <- solve(Xt_Si %*% x.matrix)
    P     <- SIGMA_inv - t(Xt_Si) %*% Q %*% Xt_Si
    Py    <- P %*% y.vec
    beta.REML <- Q %*% Xt_Si %*% y.vec
    resid <- y.vec - x.matrix %*% beta.REML
    USAE_Eblup<- data.frame(matrix(x.matrix %*% beta.REML +GIn %*% SIGMA_inv %*% resid,n, r))
    colnames(USAE_Eblup) = varnames_Y
    std.err.beta   <- sqrt(diag(Q))
    tvalue         <- beta.REML/std.err.beta
    pvalue         <- 2 * pnorm(abs(tvalue), lower.tail = FALSE)
    coef           <- cbind(beta.REML, std.err.beta, tvalue, pvalue)
    colnames(coef) = c("beta", "std.error", "t.statistics","p.value")
    Bi  <- RIn%*%solve(SIGMA)
    m   <- dim(x.matrix)[1]
    p   <- dim(x.matrix)[2]
    I   <- diag(m)
    g1d <- diag(Bi%*%GIn)
    g2d <- diag(Bi%*%x.matrix%*%Q%*%t(x.matrix)%*%t(Bi))
    dg  <- SIGMA_inv - (I-Bi) %*% SIGMA_inv
    g3d <- diag(dg %*% SIGMA %*% t(dg))/F

    W   <- diag(as.vector(w.matrix))
    g4.a<- matrix(0,m,m)
    for (i in 1:r){
      g4.a <- g4.a + d.SIGMA[[i]]*sum(diag(d.SIGMA[[i]]%*%W%*%Bi%*%SIGMA%*%t(W)%*%t(Bi)))
    }
    g4.a   <- diag(g4.a)

    g4.b   <- matrix(0,r,r)
    for (k in 0:(r-1)) {
      for (i in (k*30+1):((k+1)*n)) {
        xdi<- matrix(x.matrix[i, ], nrow = 1, ncol = p)
        for (j in (k*30+1):((k+1)*n)) {
          xdj <- matrix(x.matrix[j, ], nrow = 1, ncol = p)
          g4.b<- g4.b + d.sigma[[k+1]]*as.numeric(W[i,i]*W[j,j]*Bi[i,i]*Bi[j,j]*as.numeric(xdi %*% Q %*% t(xdj)))
        }
      }
    }
    g4.b   <- diag(kronecker(g4.b,In))
    g4d    <- g4.a - g4.b


    g4.a  <- as.data.frame(matrix(g4.a,n,r))
    g4.a  <- apply(g4.a, 2, median)
    names(g4.a) <- varnames_Y
    g4.b  <- as.data.frame(matrix(g4.b,n,r))
    g4.b  <- apply(g4.b, 2, median)
    names(g4.b) <- varnames_Y

    MSE_Eblup <- g1d + g2d + 2 * g3d
    MSE_DB    <- g1d + g2d + 2 * g3d + g4d

    MSE_Eblup <- data.frame(matrix(MSE_Eblup, n, r))
    MSE_DB    <- data.frame(matrix(MSE_DB, n, r))
    names(MSE_Eblup) = varnames_Y
    names(MSE_DB)    = varnames_Y

  } else {

    Varu <- apply(matrix(diag(RIn), n, r), 2, median)
    k    <- 0
    diff <- rep(PRECISION + 1, r)
    while (any(diff > rep(PRECISION, r)) & (k < MAXITER)) {
      Varu1 <- Varu
      if (r == 1) {
        G <- Varu1
      }
      else {
        G <- diag(as.vector(Varu1))
      }
      GIn   <- kronecker(G, In)
      SIGMA <- GIn + RIn
      SIGMA_inv <- solve(SIGMA)
      Xt_Si <- t(SIGMA_inv %*% x.matrix)
      Q     <- solve(Xt_Si %*% x.matrix)
      P     <- SIGMA_inv - t(Xt_Si) %*% Q %*% Xt_Si
      Py    <- P %*% y.vec
      s     <- vector()
      for (i in 1:r){s[i] <- (-0.5) * sum(diag(P %*% d.SIGMA[[i]])) + 0.5 * (t(Py) %*% d.SIGMA[[i]] %*% Py)}
      F     <- matrix(NA,r,r)
      for (i in 1:r){
        for (j in 1:r){
          F[j,i] <-  0.5*sum(diag(P %*% d.SIGMA[[i]] %*% P %*% d.SIGMA[[j]]))
        }
      }

      Varu  <- Varu1 + solve(F) %*% s
      diff  <- abs((Varu - Varu1)/Varu1)
      k     <- k + 1
    }
    Varu    <- as.vector(sapply(Varu, max,0))
    names(Varu) = varnames_Y
    if (k >= MAXITER && diff >= PRECISION) {
      convergence = FALSE
    }
    if (r == 1) {
      G <- Varu
    }else {
      G <- diag(as.vector(Varu))
    }
    GIn   <- kronecker(G, In)
    SIGMA <- GIn + RIn
    SIGMA_inv <- solve(SIGMA)
    Xt_Si <- t(SIGMA_inv %*% x.matrix)
    Q     <- solve(Xt_Si %*% x.matrix)
    P     <- SIGMA_inv - t(Xt_Si) %*% Q %*% Xt_Si
    Py    <- P %*% y.vec
    beta.REML <- Q %*% Xt_Si %*% y.vec
    resid <- y.vec - x.matrix %*% beta.REML
    USAE_Eblup <- data.frame(matrix(x.matrix %*% beta.REML + GIn %*% SIGMA_inv %*% resid,n, r))
    colnames(USAE_Eblup) = varnames_Y
    std.err.beta <- sqrt(diag(Q))
    tvalue       <- beta.REML/std.err.beta
    pvalue       <- 2 * pnorm(abs(tvalue), lower.tail = FALSE)
    coef         <- cbind(beta.REML, std.err.beta, tvalue, pvalue)
    colnames(coef)= c("beta", "std.error", "t.statistics","p.value")
    F_inv <- solve(F)
    Bi    <- RIn%*%solve(SIGMA)
    m     <- dim(x.matrix)[1]
    p     <- dim(x.matrix)[2]
    I     <- diag(m)

    g1d   <- diag(Bi%*%GIn)
    g2d   <- diag(Bi%*%x.matrix%*%Q%*%t(x.matrix)%*%t(Bi))

    dg    <- lapply(d.SIGMA, function(x) x %*% SIGMA_inv - GIn %*%
                      SIGMA_inv %*% x %*% SIGMA_inv)
    g3d = list()
    for (i in 1:r) {
      for (j in 1:r) {
        g3d[[(i - 1) * r + j]] = F_inv[i, j] * (dg[[i]] %*%
                                                  SIGMA %*% t(dg[[j]]))
      }
    }
    g3d   <- diag(Reduce("+", g3d))

    W     <- diag(as.vector(w.matrix))
    g4.a  <- matrix(0,m,m)
    for (i in 1:r){
      g4.a<- g4.a + d.SIGMA[[i]]*sum(diag(d.SIGMA[[i]]%*%W%*%Bi%*%SIGMA%*%t(W)%*%t(Bi)))
    }
    g4.a  <- diag(g4.a)
    g4.b  <- matrix(0,r,r)
    for (k in 0:(r-1)) {
      for (i in (k*30+1):((k+1)*n)) {
        xdi <- matrix(x.matrix[i, ], nrow = 1, ncol = p)
        for (j in (k*30+1):((k+1)*n)) {
          xdj <- matrix(x.matrix[j, ], nrow = 1, ncol = p)
          g4.b<- g4.b + d.sigma[[k+1]]*as.numeric(W[i,i]*W[j,j]*Bi[i,i]*Bi[j,j]*as.numeric(xdi %*% Q %*% t(xdj)))
        }
      }
    }
    g4.b  <- diag(kronecker(g4.b,In))
    g4d   <- g4.a - g4.b

    g4.a  <- as.data.frame(matrix(g4.a,n,r))
    g4.a  <- apply(g4.a, 2, median)
    names(g4.a) <- varnames_Y
    g4.b  <- as.data.frame(matrix(g4.b,n,r))
    g4.b  <- apply(g4.b, 2, median)
    names(g4.b) <- varnames_Y

    MSE_Eblup <- g1d + g2d + 2 * g3d
    MSE_DB    <- g1d + g2d + 2 * g3d + g4d

    MSE_Eblup <- data.frame(matrix(MSE_Eblup, n, r))
    MSE_DB    <- data.frame(matrix(MSE_DB, n, r))
    names(MSE_Eblup) = varnames_Y
    names(MSE_DB)    = varnames_Y
  }

  y.direct  <- matrix(y.vec, n,r)
  colnames(y.direct) = varnames_Y
  W <- as.matrix(w.matrix)

  USAE_DB <- matrix(NA,n,r)
  for (i in 1:r) {
    USAE_DB[,i] <- USAE_Eblup[,i] + (sum(y.direct[,i]*W[,i])-sum(USAE_Eblup[,i]*W[,i]))
  }
  colnames(USAE_DB) <- varnames_Y
  USAE_DB <- as.data.frame(USAE_DB)

  Aggregation_Direct <- colSums(as.matrix(y.direct)*(W))
  Aggregation_DB     <- colSums(as.matrix(USAE_DB)*(W))
  Aggregation_EBLUP  <- colSums(as.matrix(USAE_Eblup)*(W))

  Aggregation <- matrix(unlist(rbind(Aggregation_Direct,Aggregation_DB,Aggregation_EBLUP)),3,r)
  rownames(Aggregation) <- c("Aggregation_Direct","Aggregation_DB","Aggregation_EBLUP")
  colnames(Aggregation) <- varnames_Y


  result$USAE_Eblup = signif(USAE_Eblup, digits = 5)
  result$MSE_Eblup = signif(MSE_Eblup, digits = 5)
  result$fit$method = "REML"
  result$fit$convergence = convergence
  result$fit$iterations = k
  result$fit$estcoef = signif(coef, digits = 5)
  result$fit$refvar = signif(data.frame(t(Varu)), digits = 5)
  result$fit$informationFisher = signif(F, digits = 5)
  result$difference_benchmarking$Estimation = USAE_DB
  result$difference_benchmarking$Aggregation = Aggregation
  result$difference_benchmarking$MSE_DB = MSE_DB
  result$difference_benchmarking$g4.a = g4.a
  result$difference_benchmarking$g4.b = g4.b
  return(result)
}
