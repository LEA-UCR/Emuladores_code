likelihoodGaussian  <- function(nu,phi,beta0,beta1,sigma2,taue,model,type){
  Sigma <- cExpMat(hh,hh,type,phi,sigma2,nu, 
                   acau=acau, bcau=bcau)
  Y <- hh$Y
  if(model == "SVC"){
    X <- as.vector(scale(hh$X))
    XX <- diag(hh$X)
    muhat <- beta0+beta1*X
  } else {
    X <- hh$X
    XX <- diag(length(X))
    muhat <- beta0+beta1*X
  }
  Sigmainv <- chol2inv(chol(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))
  m2logv <- log(det(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))+
    t(Y-muhat)%*%Sigmainv%*%(Y-muhat)
  return(m2logv)
}

W_maker <- function(nCov){
  Wlist <- list()
  WSlist <- list()
  # Initialize matrices
  for(m in 0:nn){
    Wlist[[m+1]] <- list()
    WSlist[[m+1]] <- list()
    for(l in 0:m){
      Wlist[[m+1]][[l+1]] <- list()
    }
  }
  
  for(m in 0:nn){
    show(m)
    for(l in 0:m){
      if(l==0){
        k <- 0
        Wlist[[m+1]][[l+1]][[k+1]] <- Matrix(cExpMat_mult(Qlist[[m+1]],
                                                   Qlist[[l+1]],A,nCov,
                                                   types,range=phi,nu=nu, 
                                                   alpha=acau, beta=bcau))
        WSlist[[m+1]][[k+1]] <- Matrix(cExpMat_mult(Qlist[[nn+2]],
                                                    Qlist[[m+1]],A,nCov,
                                                    types,range=phi,nu=nu, 
                                                    alpha=acau, beta=bcau))
      }else{
        Wlist[[m+1]][[l+1]][[1]] <- Matrix(cExpMat_mult(Qlist[[m+1]],
                                                          Qlist[[l+1]],A,nCov,
                                                          types,range=phi,nu=nu, 
                                                          alpha=acau, beta=bcau))
        for(k in 0:(l-1)){
          Wlist[[m+1]][[l+1]][[k+2]] <- (Wlist[[m+1]][[l+1]][[k+1]]-Wlist[[m+1]][[k+1]][[k+1]]%*%solve(Wlist[[k+1]][[k+1]][[k+1]])%*%t(Wlist[[l+1]][[k+1]][[k+1]]))*Tm[[m+1]][[l+1]][[k+2]]
          WSlist[[m+1]][[k+2]] <- (WSlist[[m+1]][[k+1]]-WSlist[[k+1]][[k+1]]%*%solve(Wlist[[k+1]][[k+1]][[k+1]])%*%t(Wlist[[m+1]][[k+1]][[k+1]]))*TSm[[m+1]][[k+2]]
        }
      }
    }
  }
  Lambda_list <- purrr::map(1:(nn+1),~return(Wlist[[.]][[.]][[.]]))
  B_list <- purrr::map(1:(nn+1),~return(WSlist[[.]][[.]]))
  matrices_r <- list(Lambda=Lambda_list,B=B_list)
  return(matrices_r)
}

likelihoodMRA <- function(nu,phi,beta,A,nCov,taue,model,type,nMRA,Y,X,XR){
  XR <- as.matrix(XR)
  XX <- bdiag(purrr::map(1:dim(XR)[1],~t(XR[.,])))
  muhat <- X%*%beta
  
  Y <- as.matrix(Y) 
  muhat <- as.matrix(muhat)
  matrices <- W_maker(nCov)
  Lambda_m <- matrices$Lambda
  B_m <- matrices$B
  
  quad_SigmaY <- 0
  
  browser()
  for(k in nMRA:0){
    show(k)
    if(k==nMRA){
      Sigma_w <- B_m[[k+1]]%*%solve(Lambda_m[[k+1]])%*%t(B_m[[k+1]]) 
      XSigmae <- XX%*%Sigma_w%*%t(XX)+(1/taue) * diag(dim(XX)[1])
      SigmaYinv <- solve(XSigmae)
      logSigmaYdet <- log(det(XSigmae))
    }else{
      SigmaYinv_old <- SigmaYinv
      SigmaYinv <- SigmaYinv_old - SigmaYinv_old %*% 
        XX %*% B_m[[k+1]] %*% 
        solve(Lambda_m[[k+1]]+t(B_m[[k+1]])%*%
                t(XX) %*% SigmaYinv_old %*%
                XX %*% B_m[[k+1]])%*%t(B_m[[k+1]]) %*% t(XX) %*% SigmaYinv_old
      logSigmaYdet <- logSigmaYdet+log(det(Lambda_m[[k+1]]+
                                             t(B_m[[k+1]])%*%
                                             t(XX) %*% SigmaYinv_old %*%
                                             XX %*% B_m[[k+1]]))-
        log(det(Lambda_m[[k+1]]))
      }
  }
  m2logv <- logSigmaYdet+t(Y-muhat)%*%SigmaYinv%*%(Y-muhat)
  return(m2logv)
}
