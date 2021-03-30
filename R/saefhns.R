#' @title EBLUPs under Univariate Fay Herriot Model for non-sampled area
#' @description This function produces EBLUPs, MSE, and aggregation of Univariate SAE for non-sampled area
#' @param formula List of formula that describe the fitted model
#' @param vardir  Sampling variances of direct estimations included in data frame as the vector with the name of sampling variances in order : \code{var1, cov12,.,cov1r,var2,cov23,.,cov2r,.,cov(r-1)(r),var(r)}
#' @param cluster cluster information
#' @param nonsample A column with logical values, \code{TRUE} if the area is non-sampled
#' @param samevar Whether the variances of the data are same or not. Logical input with default \code{FALSE}
#' @param MAXITER Maximum number of iteration in Fisher-scoring algorithm with default \code{100}
#' @param PRECISION Limit of Fisher-scoring convergence tolerance with default \code{1e-4}
#' @param data The data frame
#'
#' @return This function returns a list of the following objects:
#'    \item{SAE_Eblup_sampled}{A dataframe with the values of the EBLUPs estimators for sampled areas}
#'    \item{SAE_Eblup_all}{A dataframe with the values of the EBLUPs estimators for all areas}
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
#' cluster = c("clY1","clY2","clY3")
#' nonsample = "nonsample"
#' saeFHns <- saefhns(formula,vardir,cluster, nonsample, data=datamsaeDBns)
#'
#' @export saefhns
#'
saefhns <- function (formula, vardir,cluster,nonsample, samevar = FALSE, MAXITER = 100, PRECISION = 1e-04,
                     data) {

  result = list(SAE_Eblup_sampled = NA,SAE_Eblup_all = NA, MSE_Eblup_sampled = NA, MSE_Eblup_all = NA,   randomEffect_sampled = NA, randomEffect_all = NA, Rmatrix_sampled = NA,
                fit = list(method = NA, convergence = NA, iterations = NA,
                           estcoef = NA, refvar = NA, informationFisher = NA)
                )

  if (!(TRUE %in% data[,"nonsample"])){
    stop("this saefhns() function is used for at least 1 observation has no sample, check your 'nonsample' variable or use saefh() instead ")

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
  n = length(y.vec)/r
  dfnonsample <- data[which(data$nonsample == TRUE), ]
  anns <- dfnonsample$index #area number non sample
  vardirns   <- dfnonsample[,vardir]
  x.matrixns <- formula
  for(i in 1:r) {x.matrixns[[i]] <- model.matrix(formula[[i]], na.action = na.omit, dfnonsample)}
  x.matrixns = Reduce(adiag,x.matrixns)

  if (any(is.na(data[, vardir])))
    stop("Object vardir contains NA values.")
  if (!all(vardir %in% names(data)))
    stop("Object vardir is not appropiate with data")
  if (length(vardir) != sum(1:r))
    stop("Length of vardir is not appropiate, the length must be ", sum(1:r))

  if(r==1){
    RIn = diag(as.vector(data_sampled[, vardir]))
  } else{
    RIn = RIn_function(data_sampled[, vardir],n,r)*diag(n*r)
  }


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
    SAE_Eblup<- data.frame(matrix(x.matrix %*% beta.REML +GIn %*% SIGMA_inv %*% resid,n, r))
    colnames(SAE_Eblup) = varnames_Y
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


    MSE_Eblup <- g1d + g2d + 2 * g3d

    MSE_Eblup <- data.frame(matrix(MSE_Eblup, n, r))
    names(MSE_Eblup) = varnames_Y

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
    }
    else {
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
    SAE_Eblup <- data.frame(matrix(x.matrix %*% beta.REML +GIn %*% SIGMA_inv %*% resid,n, r))
    colnames(SAE_Eblup) = varnames_Y
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


    MSE_Eblup <- g1d + g2d + 2 * g3d
    MSE_Eblup <- data.frame(matrix(MSE_Eblup, n, r))
    names(MSE_Eblup) = varnames_Y

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

  totalArea <- dim(SAE_Eblup)[1]+dim(EBLUPCluster)[1]
  idx <- cbind(index = rep(1:totalArea))
  SAE_Eblup_temp <- cbind(index=data_sampled$index,SAE_Eblup)

  SAE_Eblup_all <- merge(x=as.matrix(idx), y=SAE_Eblup_temp, by = "index", all.x=TRUE)
  SAE_Eblup_all <- SAE_Eblup_all[,-1]

  for (i in 1:totalArea){
    for (j in 1:length(anns)){
      if (i==anns[j]){
        if(r==1){
          SAE_Eblup_all[i] <- EBLUPCluster[which(EBLUPCluster$an == i), ][,-1]
        }else{
          SAE_Eblup_all[i,] <- EBLUPCluster[which(EBLUPCluster$an == i), ][,-1]
        }
      }
    }
  }

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

  if (r==1){
    vardirns <- cbind(dfnonsample[, vardir],dfnonsample[,cluster],an =dfnonsample$index)
  } else {
    vardirns <- cbind(dfnonsample[, vardir][,1:r],dfnonsample[,cluster],an =dfnonsample$index)
  }

  vardir.ns <- matrix(0,length(anns),r)

  for (i in 1:length(anns)){
    for (j in 1:r){
      vardir.ns[i,j] <- avVardir[vardirns[i,j+r],j+1]
    }
  }

  RInns <- diag(as.vector(vardir.ns))
  if (r==1){
    GInns <- kronecker(Varu,diag(length(anns)))
  } else {
    GInns <- kronecker(diag(Varu),diag(length(anns)))
  }
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
        if (r==1){
          MSE_Eblup_all[i] <- MSE_Eblupns[which(MSE_Eblupns$an == i), ][,-1]
        }else{
          MSE_Eblup_all[i,] <- MSE_Eblupns[which(MSE_Eblupns$an == i), ][,-1]
        }
      }
    }
  }


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
        if (r==1){
          randomEffect_all[i] <- randomEffectns[which(randomEffectns$an == i), ][,-1]
        } else {
          randomEffect_all[i,] <- randomEffectns[which(randomEffectns$an == i), ][,-1]
        }
      }
    }
  }
  result$SAE_Eblup_sampled = SAE_Eblup
  result$SAE_Eblup_all = SAE_Eblup_all
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
  return(result)
}
