likelihoodBanerjee <- function(nu,phi,beta0,beta1,sigma2,taue,model,type){
  C <- cExpMat(knotsMRA[[1]],hh,type,phi,variance=sigma2,nu)
  Cstar <- cExpMat(knotsMRA[[1]],knotsMRA[[1]],type,phi,variance=sigma2,nu)
  Sigma <- t(C) %*% chol2inv(chol(Cstar)) %*% C 
  Y <- hh$Y
  if(model == "SVC"){
    X <- as.vector(scale(hh$X))
    XX <- diag(X)
    muhat <- beta0+beta1*X
  } else {X <- hh$X
  XX <- diag(X)
  muhat <- beta0+beta1*X
  }
  Sigmainv <- chol2inv(chol(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))
  m2logv <- log(det(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))+
    t(Y-muhat)%*%Sigmainv%*%(Y-muhat)
  return(as.numeric(m2logv))
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
    XX <- diag(hh$X)
    muhat <- beta0+beta1*X
  } else {
    X <- hh$X
    XX <- diag(X)
    muhat <- beta0+beta1*X
  }
  
  C <- cExpMat(knotsMRA[[1]],hh,type,phi,variance=sigma2,nu)
  Cstar <- cExpMat(knotsMRA[[1]],knotsMRA[[1]],type,phi,variance=sigma2,nu)
  Sigmaw <- t(C) %*% chol2inv(chol(Cstar)) %*% C 
  Sigma <- cExpMat(hh,hh,type,phi,variance=sigma2,nu)
  Kappa <- Blockmatrix(hh$iK2)
  Sigmae <- (Sigma-Sigmaw)*Kappa
  XSigmae <- t(XX)%*%Sigmae%*%XX+(1/taue)*diag(dim(Sigmae)[1])
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
  Sigma <- cExpMat(hh,hh,type,phi,variance=sigma2,nu)
  Y <- hh$Y
  if(model == "SVC"){
    X <- as.vector(scale(hh$X))
    XX <- diag(hh$X)
    muhat <- beta0+beta1*X
  } else {
    X <- hh$X
    XX <- diag(X)
    muhat <- beta0+beta1*X
  }
  
  Sigmainv <- chol2inv(chol(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))
  m2logv <- log(det(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))+
    t(Y-muhat)%*%Sigmainv%*%(Y-muhat)
  return(as.numeric(m2logv))
}

WQXYmakerMatern <- function(sigma2){
  Qlist <- list()
  Wlist <- list()
  Xlist <- list()
  Ylist <- list()
  for(m in 0:(nn-1)){
    M <- m+1
    Qlist[[M]] <- list()
    Xlist[[M]] <- list()
    Ylist[[M]] <- list()
    Wlist[[M]] <- list()
    for(l in 1:M){
      Wlist[[M]][[l]] <- list()
    }
    for(jm in 1:(dim(indicesW[[M]])[1])){
      Qlist[[M]][[jm]] <- knotsMRA[[M]] %>% 
        filter(.data[[paste0('iP',M)]]==jm)
      Xlist[[M]][[jm]] <- diag((hh %>% 
            filter(.data[[paste0('iP',M)]]==jm) %>%
            dplyr::select(X))$X)
      Ylist[[M]][[jm]] <- (hh %>%
            filter(.data[[paste0('iP',M)]]==jm) %>%
            dplyr::select(Y))$Y
      indicesjerarq <- Qlist[[M]][[jm]] %>% 
        dplyr::select(starts_with('iP'))%>%
        st_drop_geometry()
      for(l in 1:M){
        #        show(l)
        jl <- as.numeric(indicesjerarq %>% 
             dplyr::select(.data[[paste0('iP',l)]]) %>% 
               unique())
        factorW <- 0
        if(l!=1){
          factorW <- 0
          for(k in 1:(l-1)){
            jk <- as.numeric(indicesjerarq %>% 
               dplyr::select(.data[[paste0('iP',k)]]) %>% 
                 unique())  
  #diag(Wlist[[k]][[k]][[jk]])<-diag(Wlist[[k]][[k]][[jk]])+
            # rep(sigma2,dim(Wlist[[k]][[k]][[jk]])[1])
            factorW <- factorW + Wlist[[M]][[k]][[jm]]%*%
              chol2inv(chol(Wlist[[k]][[k]][[jk]]))%*%
              t(Wlist[[l]][[k]][[jl]])
          }
        }
        Wlist[[M]][[l]][[jm]] <- cExpMat(Qlist[[M]][[jm]],
                                         Qlist[[l]][[jl]],
                                         type,phi,variance=sigma2,nu)-factorW
    rownames(Wlist[[M]][[l]][[jm]]) <- as.character(Qlist[[M]][[jm]]$indice_m)
    colnames(Wlist[[M]][[l]][[jm]]) <- as.character(Qlist[[l]][[jl]]$indice_m)
        
        #Wlist[[M]][[l]][[jm]] <- corrMaternduo(Qlist[[M]][[jm]],
        #                                       Qlist[[l]][[jl]],
        #                                       kappa,
        #                                       sigma2)-factorW
        #        Wlist[[M]][[l]][[jm]] <- corrMaternduo_fields(Qlist[[M]][[jm]],
        #                                               Qlist[[l]][[jl]],sigma2)-factorW
        
        #image.plot(Wlist[[M]][[l]][[jm]],legend.lab = paste(M,l,sep = '-'))
      }
    }
  }
  matrices_r <- list(W=Wlist,X=Xlist,Y=Ylist)
  return(matrices_r)
}

likelihoodMRA <- function(nu,phi,beta0,beta1,sigma2,taue,model,
                          type, MRA_num){
  sigma2 <- 1/taub
  matrices <- WQXYmakerMatern(sigma2)
  Y <- hh$Y
  if(model == "SVC"){
    X <- as.vector(scale(hh$X))
    XX <- diag(hh$X)
    muhat <- beta0+beta1*X
  } else {
    X <- hh$X
    XX <- diag(X)
    muhat <- beta0+beta1*X
  }
  
  Wmat <- matrices$W
  
  MRA.decompose <- function(j){
    C <- matrix(0,nrow = max(hh$indice_m),
                ncol = max(knotsMRA[[j]]$indice_m))
    rownames(C) <- as.character(1:max(hh$indice_m))
    colnames(C) <- as.character(1:max(knotsMRA[[j]]$indice_m))
    for(i in 1:length(Wmat[[4]][[j]])){
      C[rownames(Wmat[[4]][[j]][[i]]),
        colnames(Wmat[[4]][[j]][[i]])] <- Wmat[[4]][[j]][[i]]
    }
    
    Cstar <- matrix(0,nrow = max(knotsMRA[[j]]$indice_m),
                    ncol = max(knotsMRA[[j]]$indice_m))
    rownames(Cstar) <- as.character(1:max(knotsMRA[[j]]$indice_m))
    colnames(Cstar) <- rownames(Cstar)
    for(i in 1:length(Wmat[[j]][[j]])){
      Cstar[rownames(Wmat[[j]][[j]][[i]]),
            colnames(Wmat[[j]][[j]][[i]])] <- Wmat[[j]][[j]][[i]]
    }
    
    Cstar <- as(Cstar ,'dgCMatrix')
    blocks_Cstar <- unique(ExtractBlocks(Cstar))
    blocks_Cstarinv <- purrr::map(blocks_Cstar,~chol2inv(chol(.)))
    Cstarinv <- bdiag(blocks_Cstarinv)
    Sigmaw <- C %*% Cstarinv %*% t(C)
    return(list(Sigmaw,C,Cstar))
  }
  
  for(k in MRA_num:1){
    matrices_MRA <- MRA.decompose(k)
    Sigmaw <- matrices_MRA[[1]]
    C <- matrices_MRA[[2]] 
    Cstar <- matrices_MRA[[3]]
    if(k==MRA_num){
      SigmaB <- Sigmaw
      XSigmae <- XX %*% Sigmaw %*% XX + (1/taue) * 
        diag(dim(Sigmaw)[1])
      blocks_XSigmae <- unique(ExtractBlocks(XSigmae))
      blocks_inv <- purrr::map(blocks_XSigmae,~chol2inv(chol(.)))
      SigmaYinv <- XSigmae_inv <- bdiag(blocks_inv)
      blocks_det <- purrr::map_dbl(blocks_XSigmae,~det(.))
      logSigmaYdet <- sum(log(blocks_det))
    }else{
      SigmaB <- SigmaB+Sigmaw
      SigmaYinv <- SigmaYinv-SigmaYinv%*%t(XX)%*%C%*%
        solve(Cstar+t(C)%*%XX%*%SigmaYinv%*%
                t(XX)%*%C)%*%t(C)%*%XX%*%SigmaYinv
      logSigmaYdet <- log(det(Cstar+t(C)%*%XX%*%SigmaYinv%*%
                                t(XX)%*%C))-log(det(Cstar))+
        logSigmaYdet
    }
  }
  m2logv <- logSigmaYdet+t(Y-muhat)%*%SigmaYinv%*%(Y-muhat)
  return(as.numeric(m2logv))
  #return(list(m2logv,SigmaB))
}
