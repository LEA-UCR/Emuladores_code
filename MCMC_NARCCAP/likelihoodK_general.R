likelihoodBanerjee <- function(nu,phi,beta0,beta1,sigma2,taue,model,type){
  C <- cExpMat(knotsMRA[[1]],hh,type,phi,sigma2,nu, 
               acau=acau, bcau=bcau)
  Cstar <- cExpMat(knotsMRA[[1]],knotsMRA[[1]],type,phi,sigma2,nu, 
                   acau=acau, bcau=bcau)
  Sigma <- t(C) %*% chol2inv(chol(Cstar)) %*% C 
  Y <- hh$Y
  if(model == "SVC"){
    X <- as.vector(scale(hh$X))
    XX <- diag(X)
    muhat <- beta0+beta1*X
  } else {X <- hh$X
  XX <- diag(length(X))
  muhat <- beta0+beta1*X
  }
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

ExtractBlocks_pre <- function(mat, plot.graph = FALSE) {
  stopifnot(nrow(mat) == ncol(mat))
  x <- mat
  diag(x) <- 1
  edges <- as.matrix(summary(x)[c("i", "j")])
  library(igraph)
  browser()
  g <- graph.edgelist(edges, directed = FALSE)
  if (plot.graph) plot(g)
  groups <- unique(Map(sort, neighborhood(g, nrow(mat))))
  sub.Mat <- Map(`[`, list(as.matrix(mat)), groups, groups, drop = FALSE)
  sub.mat <- Map(as.matrix, sub.Mat)
  return(unique(sub.mat))
}

ExtractBlocks <- function(mat, plot.graph = FALSE) {
  stopifnot(nrow(mat) == ncol(mat))
  x <- mat
  diag(x) <- 1
  edges <- as.matrix(summary(x)[c("i", "j")])
  library(igraph)
  g <- graph.edgelist(edges, directed = FALSE)
  if (plot.graph) plot(g)
  groups <- unique(Map(as.matrix,(Map(sort, ego(g, 1)))))
  sub.Mat <- Map(`[`, list(as.matrix(mat)), groups, groups, drop = FALSE)
  sub.mat <- Map(as.matrix, sub.Mat)
  return(unique(sub.mat))
}

likelihoodFSA_Block <- function(nu,phi,beta0,beta1,sigma2,taue,model,type){
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

  C <- cExpMat(knotsMRA[[1]],hh,type,phi,sigma2,nu, 
               acau=acau, bcau=bcau)
  Cstar <- cExpMat(knotsMRA[[1]],knotsMRA[[1]],type,phi,sigma2,nu,
                   acau=acau, bcau=bcau)
  Sigmaw <- t(C) %*% chol2inv(chol(Cstar)) %*% C 
  Sigma <- cExpMat(hh,hh,type,phi,sigma2,nu, 
                   acau=acau, bcau=bcau)
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

WQXYmakerMatern <- function(){
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
      #Xlist[[M]][[jm]] <- diag((hh %>% 
      #                            filter(.data[[paste0('iP',M)]]==jm) %>%
      #                            dplyr::select(X))$X)
      #Ylist[[M]][[jm]] <- (hh %>%
      #                       filter(.data[[paste0('iP',M)]]==jm) %>%
      #                       dplyr::select(Y))$Y
      indicesjerarq <- Qlist[[M]][[jm]] %>% 
        dplyr::select(starts_with('iP'))%>%
        st_drop_geometry()
      for(l in 1:M){
        #        show(l)
        jl <- as.numeric(indicesjerarq %>% 
                           dplyr::select(.data[[paste0('iP',l)]]) %>% unique())
        factorW <- 0
        if(l!=1){
          factorW <- 0
          for(k in 1:(l-1)){
            jk <- as.numeric(indicesjerarq %>% 
                               dplyr::select(.data[[paste0('iP',k)]]) %>% unique())  
            #diag(Wlist[[k]][[k]][[jk]])<-diag(Wlist[[k]][[k]][[jk]])+
            # rep(sigma2,dim(Wlist[[k]][[k]][[jk]])[1])
            factorW <- factorW + Wlist[[M]][[k]][[jm]]%*%
              chol2inv(chol(Wlist[[k]][[k]][[jk]]))%*%
              t(Wlist[[l]][[k]][[jl]])
          }
        }
        Wlist[[M]][[l]][[jm]] <- cExpMat(Qlist[[M]][[jm]],
                                         Qlist[[l]][[jl]],
                                         type,phi,sigma2,nu, 
                                         acau=acau, bcau=bcau)-factorW
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
  #matrices_r <- list(W=Wlist,X=Xlist,Y=Ylist)
  matrices_r <- list(W=Wlist)
  return(matrices_r)
}

likelihoodMRA <- function(nu,phi,beta0,beta1,sigma2,taue,model,type,nMRA){
  source('likelihoodK_general.R')
  library(Matrix)
  library(plot.matrix)
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
  
  rownames(XX) <- colnames(XX) <- as.character(1:max(hh$indice_m))
  Y <- as.matrix(Y)
  muhat <- as.matrix(muhat)
  rownames(Y) <- rownames(muhat) <- rownames(XX)
  
  matrices <- WQXYmakerMatern()
  Wmat <- matrices$W

  MRA.decompose <- function(j){
    
    blocks_C <- list()
    blocks_Cstar <- list()
    blocks_Cstarinv <- list()
    blocks_Sigmaw <- list()
    blocks_XSigmae <- list()
    blocks_XX <-list()
    blocks_Y <-list()
    blocks_muhat <-list()
    indices_blocks <- list()
#    Sigmaw <- matrix(0,nrow = max(hh$indice_m),ncol = max(hh$indice_m))
#    rownames(Sigmaw) <- colnames(Sigmaw) <- as.character(1:max(hh$indice_m))
#    C <- matrix(0,nrow = max(hh$indice_m),ncol = max(knotsMRA[[j]]$indice_m))
#    rownames(C) <- as.character(1:max(hh$indice_m))
#    colnames(C) <- as.character(1:max(knotsMRA[[j]]$indice_m))
#    Cstar <- matrix(0,nrow = max(knotsMRA[[j]]$indice_m),
#                     ncol = max(knotsMRA[[j]]$indice_m))
#    rownames(Cstar) <- as.character(1:max(knotsMRA[[j]]$indice_m))
#    colnames(Cstar) <- rownames(Cstar)

    for(i in 1:length(Wmat[[j]][[j]])){
      blocks_C[[i]] <- Wmat[[nMRA]][[j]][[i]]
      blocks_Cstar[[i]] <- Wmat[[j]][[j]][[i]]
      blocks_Cstarinv[[i]] <- chol2inv(chol(blocks_Cstar[[i]]))
      blocks_Sigmaw[[i]] <- blocks_C[[i]] %*% blocks_Cstarinv[[i]] %*% 
        t(blocks_C[[i]])
      indices_blocks[[i]] <- rownames(Wmat[[nMRA]][[j]][[i]])
      # Sigmaw[rownames(Wmat[[nMRA]][[j]][[i]]),
      #      rownames(Wmat[[nMRA]][[j]][[i]])] <- blocks_Sigmaw[[i]]
      # C[rownames(Wmat[[nMRA]][[j]][[i]]),
      #      colnames(Wmat[[nMRA]][[j]][[i]])] <- blocks_C[[i]]
      # Cstar[rownames(Wmat[[j]][[j]][[i]]),
      #      colnames(Wmat[[j]][[j]][[i]])] <- blocks_Cstar[[i]]
      blocks_XX[[i]] <- XX[indices_blocks[[i]],indices_blocks[[i]]]
      blocks_Y[[i]] <- Y[indices_blocks[[i]],]
      blocks_muhat[[i]] <- muhat[indices_blocks[[i]],]
      if(j==nMRA){
        blocks_XSigmae[[i]] <- blocks_XX[[i]] %*% blocks_Sigmaw[[i]] %*% 
          blocks_XX[[i]] + (1/taue) * diag(dim(blocks_Sigmaw[[i]])[1])
      }
    }
    
    #return(list(Sigmaw,blocks_C,blocks_Cstar,blocks_XSigmae,blocks_XX,blocks_Y,blocks_muhat,indices_blocks))
    return(list(blocks_C,blocks_Cstar,blocks_XSigmae,
                blocks_XX,blocks_Y,blocks_muhat,indices_blocks))
  }
  
  quad_SigmaY <- 0
  
  for(k in nMRA:1){
    #show(k)
    matrices_MRA <- MRA.decompose(k)
    #Sigmaw <- matrices_MRA[[1]]
    blocks_C <- matrices_MRA[[1]] 
    blocks_Cstar <- matrices_MRA[[2]]
    blocks_XSigmae <- matrices_MRA[[3]]
    blocks_XX <- matrices_MRA[[4]]
    blocks_Y <- matrices_MRA[[5]]
    blocks_muhat <- matrices_MRA[[6]]
    indices_blocks <- matrices_MRA[[7]]
    if(k==nMRA){
      #SigmaB <- Sigmaw
      blocks_SigmaYinv <- purrr::map(blocks_XSigmae,~chol2inv(chol(.)))
      SigmaYinv <- matrix(0,nrow = max(hh$indice_m),ncol = max(hh$indice_m))
      rownames(SigmaYinv) <- colnames(SigmaYinv) <- as.character(1:max(hh$indice_m))

      for(i in 1:length(blocks_SigmaYinv)){
        SigmaYinv[indices_blocks[[i]],indices_blocks[[i]]] <- blocks_SigmaYinv[[i]]
        quad_SigmaY <- quad_SigmaY+t(blocks_Y[[i]]-blocks_muhat[[i]])%*%
          blocks_SigmaYinv[[i]]%*%(blocks_Y[[i]]-blocks_muhat[[i]]) 
    #    show(quad_SigmaY)
      }
      blocks_det <- purrr::map_dbl(blocks_XSigmae,~det(.))
      logSigmaYdet <- sum(log(blocks_det))
    }else{
      #SigmaB <- SigmaB+Sigmaw
      SigmaYinv_old <- SigmaYinv
      SigmaYinv <- matrix(0,nrow = max(hh$indice_m),ncol = max(hh$indice_m))
      rownames(SigmaYinv) <- colnames(SigmaYinv) <- as.character(1:max(hh$indice_m))
      for(i in 1:length(indices_blocks)){
        blocks_SigmaYinv_old <- SigmaYinv_old[indices_blocks[[i]],indices_blocks[[i]]]
        blocks_SigmaYinv <- blocks_SigmaYinv_old - blocks_SigmaYinv_old %*% 
          blocks_XX[[i]] %*% blocks_C[[i]] %*% 
          solve(blocks_Cstar[[i]]+t(blocks_C[[i]])%*%
                  blocks_XX[[i]]%*%blocks_SigmaYinv_old%*%
                  blocks_XX[[i]]%*%blocks_C[[i]])%*%t(blocks_C[[i]])%*%
          blocks_XX[[i]]%*%blocks_SigmaYinv_old
        logSigmaYdet <- logSigmaYdet+log(det(blocks_Cstar[[i]]+
        t(blocks_C[[i]])%*%blocks_XX[[i]]%*%blocks_SigmaYinv%*%
          blocks_XX[[i]]%*%blocks_C[[i]]))-
          log(det(blocks_Cstar[[i]]))
        quad_SigmaY <- quad_SigmaY+
          t(blocks_Y[[i]]-blocks_muhat[[i]])%*%blocks_SigmaYinv%*%
          (blocks_Y[[i]]-blocks_muhat[[i]])
        SigmaYinv[indices_blocks[[i]],indices_blocks[[i]]] <- blocks_SigmaYinv
      }
      
      #SigmaYinv <- SigmaYinv-SigmaYinv%*%t(XX)%*%C%*%
      #  solve(Cstar+t(C)%*%XX%*%SigmaYinv%*%t(XX)%*%C)%*%t(C)%*%XX%*%SigmaYinv
      #logSigmaYdet <- log(det(Cstar+t(C)%*%XX%*%SigmaYinv%*%t(XX)%*%C))-log(det(Cstar))+
      #  logSigmaYdet
    }
  }
  #m2logv <- logSigmaYdet+t(Y-muhat)%*%SigmaYinv%*%(Y-muhat)
  m2logv <- logSigmaYdet+quad_SigmaY
  return(m2logv)
  #return(list(m2logv,SigmaB))
}
