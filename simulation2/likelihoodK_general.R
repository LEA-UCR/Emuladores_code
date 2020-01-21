likelihoodBanerjee <- function(nu,phi,beta0,beta1,sigma2,taue,model,type){
  C <- cExpMat(knotsMRA[[1]],hh,type,phi,sigma2,nu)
  Cstar <- cExpMat(knotsMRA[[1]],knotsMRA[[1]],type,phi,sigma2,nu)
  Sigma <- t(C) %*% chol2inv(chol(Cstar)) %*% C 
  Y <- hh$Y
  if(model == "SVC"){
    X <- as.vector(scale(hh$X))
  } else {X <- hh$X}
  XX <- diag(X)
  muhat <- beta0+beta1*X
  Sigmainv <- chol2inv(chol(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))
  m2logv <- log(det(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))+
    t(Y-muhat)%*%Sigmainv%*%(Y-muhat)
  return(m2logv)
}

Blockmatrix <- function(iP){
  ndim <- length(iP)
  imatrix <- matrix(1,nrow = ndim,ncol = ndim)
  for(i in 1:ndim){
    for(j in 1:i){
      if(iP[i]!=iP[j]) imatrix[i,j] <- imatrix[j,i] <- 0
    }
  }
  return(imatrix)
}

ExtractBlocks <- function(mat, plot.graph = FALSE) {
  stopifnot(nrow(mat) == ncol(mat))
  x <- mat
  diag(x) <- 1
  edges <- as.matrix(summary(x)[c("i", "j")])
  library(igraph)
  g <- graph.edgelist(edges, directed = FALSE)
  if (plot.graph) plot(g)
  groups <- unique(Map(sort, neighborhood(g, nrow(mat))))
  sub.Mat <- Map(`[`, list(as.matrix(mat)), groups, groups, drop = FALSE)
  sub.mat <- Map(as.matrix, sub.Mat)
  return(sub.mat)
}

likelihoodFSA_Block <- function(nu,phi,beta0,beta1,sigma2,taue,model,type){
  Y <- hh$Y
  if(model == "SVC"){
    X <- as.vector(scale(hh$X))
  } else {X <- hh$X}
  XX <- diag(X)
  muhat <- beta0+beta1*X
  C <- cExpMat(knotsMRA[[1]],hh,type,phi,sigma2,nu)
  Cstar <- cExpMat(knotsMRA[[1]],knotsMRA[[1]],type,phi,sigma2,nu)
  Sigmaw <- t(C) %*% chol2inv(chol(Cstar)) %*% C 
  Sigma <- cExpMat(hh,hh,type,phi,sigma2,nu)
  Kappa <- Blockmatrix(hh$iK2)
  Sigmae <- (Sigma-Sigmaw)*Kappa
  XSigmae <- t(XX)%*%Sigmae%*%XX+sigma2*diag(dim(Sigmae)[1])
  library(Matrix)
  XSigmae <- as(XSigmae ,'dgCMatrix')
  blocks <- unique(ExtractBlocks(XSigmae))
  blocksinv <- purrr::map(blocks,~chol2inv(chol(.)))
  blocksdet <- purrr::map_dbl(blocks,~det(.))
  XSigmaeinv <- bdiag(blocksinv)
  SigmaYinv <- XSigmaeinv-XSigmaeinv%*%t(XX)%*%t(C)%*%
    solve(Cstar+C%*%XX%*%XSigmaeinv%*%XX%*%t(C))%*%C%*%XX%*%XSigmaeinv
  SigmaYdet <- det(Cstar+C%*%XX%*%XSigmaeinv%*%XX%*%t(C))*(det(Cstar))^(-1)*
    prod(blocksdet)
  m2logv <- as.numeric(log(SigmaYdet)+t(Y-muhat)%*%SigmaYinv%*%(Y-muhat))
  return(m2logv)
}

likelihoodGaussian  <- function(nu,phi,beta0,beta1,sigma2,taue,model,type){
  Sigma <- cExpMat(hh,hh,type,phi,sigma2,nu)
  Y <- hh$Y
  if(model == "SVC"){
    X <- as.vector(scale(hh$X))
  } else {X <- hh$X}
  XX <- diag(X)
  muhat <- beta0+beta1*X
  Sigmainv <- chol2inv(chol(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))
  m2logv <- log(det(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))+
    t(Y-muhat)%*%Sigmainv%*%(Y-muhat)
  return(m2logv)
}

# likelihood <- function(nu=nu,phi,beta0=beta0,beta1=beta1,
#                        taub,taue,model=model,type=type){
#   #Sigma <- cExpMat(hh,hh,type,phi,sigma2,nu)
#   Y <- hh$Y
#   if(model == "SVC"){
#     X <- scale(hh$X)
#   } else {X <- hh$X}
#   XX <- diag(X)
#   muhat <- beta0+beta1*X
#   singlelikelihoods = dnorm(Y, mean = muhat, sd = sqrt(1/taue), log = T)
#   m2logv = sum(singlelikelihoods)
#   return(m2logv)
# }

