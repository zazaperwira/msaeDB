#' @title EBLUPs under Multivariate Fay Herriot Model with Difference Benchmarking for non-sampled area
#' @description This function produces EBLUPs, MSE, and aggregation of Multivariate SAE with Difference Benchmarking for non-sampled area
#' @param formula List of formula that describe the fitted model
#' @param vardir Sampling variances of direct estimations included in data frame as the vector with the name of sampling variances in order : \code{var1, cov12,.,cov1r,var2,cov23,.,cov2r,.,cov(r-1)(r),var(r)}
#' @param weight  Known proportion of units in small areas, where \eqn{\sum_{d=1}^{D}}{sum from d=1 to D of} \eqn{W_{rd}}{Wrd} = 1 . \code{d = 1 ... D} is the number of small areas, and \code{r = 1 ... R} is the number of response variables
#' @param cluster cluster information
#' @param nonsample A column with logical values, \code{TRUE} if the area is non-sampled
#' @param samevar Whether the variances of the data are same or not. Logical input with default \code{FALSE}
#' @param MAXITER Maximum number of iteration in Fisher-scoring algorithm with default \code{100}
#' @param PRECISION Limit of Fisher-scoring convergence tolerance with default \code{1e-4}
#' @param data The data frame
#'
#' @return This function returns a list of the following objects:
#'    \item{MSAE_Eblup_sampled}{A dataframe with the values of the EBLUPs estimators for sampled areas}
#'    \item{MSAE_Eblup_all}{A dataframe with the values of the EBLUPs estimators for all areas}
#'    \item{MSE_Eblup_sampled}{A dataframe with the values of estimated mean square errors of EBLUPs estimators for sampled areas}
#'    \item{MSE_Eblup_all}{A dataframe with the values of estimated mean square errors of EBLUPs estimators for all areas}
#'    \item{randomEffect_sampled}{a dataframe with the values of the random effect estimators for sampled areas}
#'    \item{randomEffect_all}{a dataframe with the values of the random effect estimators for all areas}
#'    \item{Rmatrix_sampled}{a block diagonal matrix composed of sampling errors for sampled areas}
#'    \item{fit}{A list containing the following objects:}
#'      \itemize{
#'        \item method : The fitting method (this function is using "REML")
#'        \item convergence : The convergence result of Fisher-scoring algorithm (Logical Value)
#'        \item iterations : The number of Fisher-Scoring algorithm iterations
#'        \item estcoef : A dataframe with the estimated model coefficient, standard error,t statistics, p-values of the significance of each coefficient
#'        \item refvar : A dataframe with estimated random effect variances
#'        \item informationFisher : A matrix of information fisher from Fisher-scoring algorithm
#'      }
#'
#'    \item{difference_benchmarking}{a list containing the following objects:}
#'      \itemize{
#'        \item Estimation_sampled : A dataframe with the values of benchmarked EBLUPs estimators for sampled areas
#'        \item Estimation_all : A dataframe with the values of benchmarked EBLUPs estimators for all areas
#'        \item Aggregation_sampled : The aggregation of benchmarked EBLUPs estimators, EBLUPs estimators and direct estimations for sampled areas
#'        \item Aggregation_all : The aggregation of benchmarked EBLUPs estimators, EBLUPs estimators and direct estimations for all areas
#'        \item MSE_DB_sampled : A dataframe with the values of estimated mean square errors of benchmarked EBLUPs estimators for sampled areas
#'        \item MSE_DB_all : A dataframe with the values of estimated mean square errors of benchmarked EBLUPs estimators for all areas
#'        \item g.4a : First component of g4 in difference benchmarking MSE estimation formula
#'        \item g.4b : Second component of g4 in difference benchmarking MSE estimation formula
#'      }
#'
#' @examples
#' ##load dataset
#' data(datamsaeDBns)
#'
#' #Compute Fitted model for Y1, Y2, and Y3
#' #Y1 ~ X1 + X2
#' #Y2 ~ X1 + X2
#' #Y3 ~ X1 + X2
#'
#' ##Using parameter 'data'
#' formula = list(f1 = Y1~X1+X2,
#'                f2 = Y2~X1+X2,
#'                f3 = Y3~X1+X2)
#' vardir = c("v1","v12","v13","v2","v23","v3")
#' weight = c("w1","w2","w3")
#' cluster = c("clY1","clY2","clY3")
#' nonsample = "nonsample"
#' msaeDBns <- msaedbns(formula,vardir, weight,cluster, nonsample, data=datamsaeDBns)
#'
#' @export msaedbns
#'
msaedbns <- function (formula, vardir, weight,cluster,nonsample, samevar = FALSE, MAXITER = 100, PRECISION = 1e-04,
                      data) {

  result = list(MSAE_Eblup_sampled = NA,MSAE_Eblup_all = NA, MSE_Eblup_sampled = NA, MSE_Eblup_all = NA,   randomEffect_sampled = NA, randomEffect_all = NA, Rmatrix_sampled = NA,
                fit = list(method = NA, convergence = NA, iterations = NA,
                           estcoef = NA, refvar = NA, informationFisher = NA),
                difference_benchmarking = list(Estimation_sampled = NA,
                                               Estimation_all = NA,
                                               Aggregation_sampled = NA,
                                               Aggregation_all = NA,
                                               MSE_DB_sampled = NA,
                                               MSE_DB_all = NA,
                                               g4.a = NA,
                                               g4.b = NA))

  if (!(TRUE %in% data[,"nonsample"])){
    stop("this msaedbns() function is used for at least 1 observation has no sample, check your 'nonsample' variable or use msaedb() instead ")

  }
  if (length(formula)<=1){
    stop("this msaedbns() function is used for at least 2 response variables, numbers of your response variables is ",length(formula),". use saedbns() function instead")
  }
  r <- length(formula)
  RIn_function <- function(vardir, n,r){
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
        }else { RIn[[it+it2]] <- diag(vardir[,it])}
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
    RIn<- return(RR)
  }
  index <- cbind(rep(1:dim(data)[1]))
  data <- cbind(index,data)
  data_sampled <- data[which(data$nonsample == FALSE), ]

  formuladata <- formula
  for(i in 1:r) {formuladata[[i]] <- model.frame(formula[[i]], na.action = na.omit, data_sampled )}

  y.vec <- unlist(lapply(formuladata, function(x){x[1][1]}))

  x.matrix <- formula
  for(i in 1:r) {x.matrix[[i]] <- model.matrix(formula[[i]], na.action = na.omit, data_sampled)}
  x.matrix = Reduce(adiag,x.matrix)
  w.matrix_temp = as.matrix(data_sampled[,weight])
  w.matrix <- matrix(0,dim(w.matrix_temp)[1],r)
  for (i in 1:dim(w.matrix)[1]) {
    for (j in 1:r){
      w.matrix[i,j] <- w.matrix_temp[i,j]/sum(w.matrix_temp[,j])
    }
  }

  n = length(y.vec)/r
  dfnonsample <- data[which(data$nonsample == TRUE), ]
  anns <- dfnonsample$index #area number non sample
  vardirns   <- dfnonsample[,vardir]
  x.matrixns <- formula
  for(i in 1:r) {x.matrixns[[i]] <- model.matrix(formula[[i]], na.action = na.omit, dfnonsample)}
  x.matrixns = Reduce(adiag,x.matrixns)

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

  RIn = RIn_function(data_sampled[, vardir],n,r)


  for (i in 1:r) {
    if (attr(attributes(formuladata[[i]])$terms, "response") == 1)
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
    MSAE_Eblup<- data.frame(matrix(x.matrix %*% beta.REML +GIn %*% SIGMA_inv %*% resid,n, r))
    colnames(MSAE_Eblup) = varnames_Y
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

    g4.b  <- matrix(0,r,r)
    for (l in 0:(r-1)) {
      for (i in ((l*n)+1):((l+1)*n)) {
        xdi <- matrix(x.matrix[i, ], nrow = 1, ncol = p)
        for (j in ((l*n)+1):((l+1)*n)) {
          xdj <- matrix(x.matrix[j, ], nrow = 1, ncol = p)
          g4.b<- g4.b + d.sigma[[l+1]]*as.numeric(W[i,i]*W[j,j]*Bi[i,i]*Bi[j,j]*as.numeric(xdi %*% Q %*% t(xdj)))
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
      G <- diag(as.vector(Varu1))
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
    G <- diag(as.vector(Varu))
    GIn   <- kronecker(G, In)
    SIGMA <- GIn + RIn
    SIGMA_inv <- solve(SIGMA)
    Xt_Si <- t(SIGMA_inv %*% x.matrix)
    Q     <- solve(Xt_Si %*% x.matrix)
    P     <- SIGMA_inv - t(Xt_Si) %*% Q %*% Xt_Si
    Py    <- P %*% y.vec
    beta.REML <- Q %*% Xt_Si %*% y.vec
    resid <- y.vec - x.matrix %*% beta.REML
    MSAE_Eblup <- data.frame(matrix(x.matrix %*% beta.REML +GIn %*% SIGMA_inv %*% resid,n, r))
    colnames(MSAE_Eblup) = varnames_Y
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
        g3d[[(i - 1) * r + j]] = F_inv[i, j] * (dg[[i]] %*% SIGMA %*% t(dg[[j]]))
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
    for (l in 0:(r-1)) {
      for (i in ((l*n)+1):((l+1)*n)) {
        xdi <- matrix(x.matrix[i, ], nrow = 1, ncol = p)
        for (j in ((l*n)+1):((l+1)*n)) {
          xdj <- matrix(x.matrix[j, ], nrow = 1, ncol = p)
          g4.b<- g4.b + d.sigma[[l+1]]*as.numeric(W[i,i]*W[j,j]*Bi[i,i]*Bi[j,j]*as.numeric(xdi %*% Q %*% t(xdj)))
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

  randomEffect <- GIn%*%SIGMA_inv%*%resid
  randomEffect <- as.data.frame(matrix(randomEffect, n, r))
  names(randomEffect) <- varnames_Y

  XBns   <- matrix(x.matrixns%*%beta.REML,length(anns),r)
  XBns   <- cbind(XBns,dfnonsample[,cluster])
  estrEff <- cbind(randomEffect,data_sampled[,cluster],nonsample=data_sampled$nonsample,data_sampled$index)
  avrEffc <- matrix(0,max(data[,cluster]),r)

  for (i in 1:r){
    for (j in 1:max(data[,cluster])){
      if (!(j %in% estrEff[,i+r])){
        newline <- rep(NA,dim(estrEff)[2])
        newline[1:r] <-0
        newline[i+r] <- j
        estrEff <- rbind(estrEff,newline)
      }
    }
  }

  for (i in 1:r){
    avrEffc[,i] <- sapply(split(estrEff[,i], estrEff[,i+r]), mean)
  }
  colnames(avrEffc) <- varnames_Y
  avrEffc <- cbind(cluster = c(1:max(data[,cluster])), avrEffc)

  EBLUPCluster <- matrix(0,length(anns),r)

  for (i in 1:length(anns)){
    for (j in 1:r){
      EBLUPCluster[i,j] <- XBns[i,j]+avrEffc[XBns[i,j+r],j+1]
    }
  }
  colnames(EBLUPCluster) <- varnames_Y
  EBLUPCluster <- data.frame(an = anns,EBLUPCluster )

  totalArea <- dim(MSAE_Eblup)[1]+dim(EBLUPCluster)[1]
  idx <- cbind(index = rep(1:totalArea))
  MSAE_Eblup_temp <- cbind(index=data_sampled$index,MSAE_Eblup)

  MSAE_Eblup_all <- merge(x=as.matrix(idx), y=MSAE_Eblup_temp, by = "index", all.x=TRUE)
  MSAE_Eblup_all <- MSAE_Eblup_all[,-1]

  for (i in 1:totalArea){
    for (j in 1:length(anns)){
      if (i==anns[j]){
        MSAE_Eblup_all[i,] <- EBLUPCluster[which(EBLUPCluster$an == i), ][,-1]
      }
    }
  }

  ######### Aggregation Sampled area #########
  y.vec  <- matrix(y.vec, n,r)
  colnames(y.vec) = varnames_Y
  W <- as.matrix(w.matrix)

  alfa <- ginv(t(W))%*%(t(W)%*%y.vec-t(W)%*%as.matrix(MSAE_Eblup))

  MSAE_DB <- as.data.frame(MSAE_Eblup + alfa)
  colnames(MSAE_DB) <- varnames_Y

  Aggregation_Direct <- colSums(as.matrix(y.vec)*(W))
  Aggregation_DB     <- colSums(as.matrix(MSAE_DB)*(W))
  Aggregation_EBLUP  <- colSums(as.matrix(MSAE_Eblup)*(W))

  Aggregation <- matrix(unlist(rbind(Aggregation_Direct,Aggregation_DB,Aggregation_EBLUP)),3,r)
  rownames(Aggregation) <- c("Aggregation_Direct","Aggregation_DB","Aggregation_EBLUP")
  colnames(Aggregation) <- varnames_Y

  #Benchmarking All sample
  formuladata_all <- formula
  for(i in 1:r) {formuladata_all[[i]] <- model.frame(formula[[i]], na.action = na.omit, data )}

  y.vec_all <- unlist(lapply(formuladata_all, function(x){x[1][1]}))
  y.vecns <- y.vec_all
  #Make sure y of non sample is equal to 0
  for (i in 1:r ){
    for (j in 1:(length(anns))){
      y.vecns[anns[j]+((i-1)*n)]<-0
    }
  }
  y.direct_all  <- matrix(y.vecns, totalArea,r)
  colnames(y.direct_all) = varnames_Y
  W_all = as.matrix(data[,weight])

  alfa_all <- ginv(t(W_all))%*%(t(W_all)%*%y.direct_all-t(W_all)%*%as.matrix(MSAE_Eblup_all))

  MSAE_DB_all <- as.data.frame(MSAE_Eblup_all + alfa_all)
  colnames(MSAE_DB_all) <- varnames_Y

  Aggregation_Direct <- colSums(as.matrix(y.direct_all)*(W_all))
  Aggregation_DB     <- colSums(as.matrix(MSAE_DB_all)*(W_all))
  Aggregation_EBLUP  <- colSums(as.matrix(MSAE_Eblup_all)*(W_all))

  Aggregation_all <- matrix(unlist(rbind(Aggregation_Direct,Aggregation_DB,Aggregation_EBLUP)),3,r)
  rownames(Aggregation_all) <- c("Aggregation_Direct","Aggregation_DB","Aggregation_EBLUP")
  colnames(Aggregation_all) <- varnames_Y

  v <- matrix(diag(RIn),n,r)
  colnames(v) <- varnames_Y
  SampledV <- cbind(v,data_sampled[,cluster],an = data_sampled$index,nonsample=data_sampled$nonsample)
  avVardir <- matrix(0,max(data[,cluster]),r)
  for (i in 1:r){
    for (j in 1:max(data[,cluster])){
      if (!(j %in% SampledV[,i+r])){
        newline <- rep(NA,dim(SampledV)[2])
        newline[1:r] <-0
        newline[i+r] <- j
        SampledV <- rbind(SampledV,newline)
      }
    }
  }

  for (i in 1:r){
    avVardir[,i] <- sapply(split(SampledV[,i], SampledV[,r+i]), mean)
  }

  colnames(avVardir) <- varnames_Y
  avVardir <- cbind(cluster = c(1:max(data[,cluster])), avVardir)

  vardirns <- cbind(dfnonsample[, vardir][,1:r],dfnonsample[,cluster],an =dfnonsample$index)
  vardir.ns <- matrix(0,length(anns),r)

  for (i in 1:length(anns)){
    for (j in 1:r){
      vardir.ns[i,j] <- avVardir[vardirns[i,j+r],j+1]
    }
  }

  RInns <- diag(as.vector(vardir.ns))
  GInns <- kronecker(diag(Varu),diag(length(anns)))
  Bins  <- RInns%*%solve(GInns+RInns)
  g1dns <- Bins%*%GInns
  g1dns <- matrix(diag(g1dns),length(anns),r)

  SIGMAns <- GInns + RInns
  SIGMA_invns <- solve(SIGMAns)
  Xt_Sins <- t(SIGMA_invns %*% x.matrixns)
  Qns     <- ginv(Xt_Sins%*%x.matrixns)
  xqxns <- (x.matrixns)%*%Qns%*%t(x.matrixns)
  g2dns   <- (Bins%*%xqxns%*%t(Bins))
  g2dns   <- matrix(diag(g2dns),length(anns),r)

  g3eblup <- matrix(g3d,n,r)
  colnames(g3eblup) <- varnames_Y
  Sampledg3d <- cbind(g3eblup,data_sampled[,cluster],an = data_sampled$index ,nonsample=data_sampled$nonsample)

  avg3d <- matrix(0,max(data[,cluster]),r)

  for (i in 1:r){
    for (j in 1:max(data[,cluster])){
      if (!(j %in% Sampledg3d[,i+r])){
        newline <- rep(NA,dim(Sampledg3d)[2])
        newline[1:r] <-0
        newline[i+r] <- j
        Sampledg3d <- rbind(Sampledg3d,newline)
      }
    }
  }

  for (i in 1:r){
    avg3d[,i] <- sapply(split(Sampledg3d[,i], Sampledg3d[,r+i]), mean)
  }
  colnames(avg3d) <- varnames_Y
  avg3d <- cbind(cluster = c(1:max(data[,cluster])), avg3d)

  g3dns <- matrix(0,length(anns),r)

  for (i in 1:length(anns)){
    for (j in 1:r){
      g3dns[i,j] <-avg3d[vardirns[i,j+r],j+1]
    }
  }

  MSE_Eblupns = g1dns+g2dns+2*g3dns
  MSE_Eblupns <- data.frame(an=anns,MSE_Eblupns)

  MSE_Eblup_temp <- cbind(index=data_sampled$index,MSE_Eblup)
  MSE_Eblup_all <- merge(x=as.matrix(idx), y=MSE_Eblup_temp, by = "index", all.x=TRUE)
  MSE_Eblup_all <- MSE_Eblup_all[,-1]

  for (i in 1:totalArea){
    for (j in 1:length(anns)){
      if (i==anns[j]){
        MSE_Eblup_all[i,] <- MSE_Eblupns[which(MSE_Eblupns$an == i), ][,-1]
      }
    }
  }

  g4dns  <- colMeans(matrix(g4d,n,r))
  g4dns <- matrix(NA,totalArea,r)
  for (i in 1:totalArea){
    g4dns[i,] <- colMeans(matrix(g4d,n,r))
  }
  MSE_DB_all <-MSE_Eblup_all+g4dns

  randomEffectns <- matrix(0,length(anns),r)
  randomEffectns <- cbind(randomEffectns,dfnonsample[,cluster],an =dfnonsample$index)

  for (i in 1:length(anns)){
    for (j in 1:r){
      randomEffectns[i,j] <- randomEffectns[i,j]+avrEffc[XBns[i,j+r],j+1]
    }
  }

  randomEffectns <- randomEffectns[,1:r]
  randomEffectns <- data.frame(an=anns,randomEffectns)

  randomEffect_temp <- cbind(index=data_sampled$index,randomEffect)
  randomEffect_all <- merge(x=as.matrix(idx), y=randomEffect_temp, by = "index", all.x=TRUE)
  randomEffect_all <- randomEffect_all[,-1]

  for (i in 1:totalArea){
    for (j in 1:length(anns)){
      if (i==anns[j]){
        randomEffect_all[i,] <- randomEffectns[which(randomEffectns$an == i), ][,-1]
      }
    }
  }
  result$MSAE_Eblup_sampled = MSAE_Eblup
  result$MSAE_Eblup_all = MSAE_Eblup_all
  result$MSE_Eblup_sampled = MSE_Eblup
  result$MSE_Eblup_all = MSE_Eblup_all
  result$randomEffect_sampled = signif(randomEffect, digits = 5)
  result$randomEffect_all = signif(randomEffect_all, digits = 5)
  result$Rmatrix_sampled = signif(RIn, digits = 5)
  result$fit$method = "REML"
  result$fit$convergence = convergence
  result$fit$iterations = k
  result$fit$estcoef = signif(coef, digits = 5)
  result$fit$refvar = signif(data.frame(t(Varu)), digits = 5)
  result$fit$informationFisher = signif(F, digits = 5)
  result$difference_benchmarking$Estimation_sampled = MSAE_DB
  result$difference_benchmarking$Estimation_all = MSAE_DB_all
  result$difference_benchmarking$Aggregation_sampled = Aggregation
  result$difference_benchmarking$Aggregation_all = Aggregation_all
  result$difference_benchmarking$MSE_DB_sampled = MSE_DB
  result$difference_benchmarking$MSE_DB_all = MSE_DB_all
  result$difference_benchmarking$g4.a = g4.a
  result$difference_benchmarking$g4.b = g4.b
  return(result)
}
