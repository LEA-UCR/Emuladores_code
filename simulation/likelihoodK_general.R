WQXYmaker <- function(){
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
      #      show(paste(m,jm,sep = '-'))
      Qlist[[M]][[jm]] <- knotsMRA[[M]] %>% 
        filter(.data[[paste0('iP',M)]]==jm)
      Xlist[[M]][[jm]] <- diag((hh %>% 
                                  filter(.data[[paste0('iP',M)]]==jm) %>%
                                  dplyr::select(Xcov))$Xcov)
      Ylist[[M]][[jm]] <- (hh %>%
                             filter(.data[[paste0('iP',M)]]==jm) %>%
                             dplyr::select(Yresp))$Yresp
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
                                         type,
                                         range,
                                         sigma2,nu)-factorW
        #        Wlist[[M]][[l]][[jm]] <- corrMaternduo_fields(Qlist[[M]][[jm]],
        #                                               Qlist[[l]][[jl]],sigma2)-factorW
        
        #image.plot(Wlist[[M]][[l]][[jm]],legend.lab = paste(M,l,sep = '-'))
      }
    }
  }
  matrices_r <- list(W=Wlist,X=Xlist,Y=Ylist)
  return(matrices_r)
}


likelihoodKatzfuss <- function(nu,range,sigma2,taue,beta,type){
  
  matrices_r <- WQXYmaker()
  Wlist <- matrices_r$W
  Xlist <- matrices_r$X
  Ylist <- matrices_r$Y
  
  dj <- list()
  uj <- list()
  Atilde <- list()
  omegatilde <- list()
  A <- list()
  omega <- list()
  
  
  for(m in (nn-1):0){
    M <- m+1
    dj[[M]] <- list()
    uj[[M]] <- list()
    Atilde[[M]] <- list()
    omegatilde[[M]] <- list()
    A[[M]] <- list()
    omega[[M]] <- list()
    if(M==nn){
      for(jm in 1:(dim(indicesW[[M]])[1])){
        Sigma_m <- Xlist[[M]][[jm]]%*%Wlist[[M]][[M]][[jm]]%*%Xlist[[M]][[jm]]+(1/taue)*diag(dim(Xlist[[M]][[jm]])[1])
        Sigma_m_inv <- chol2inv(chol(Sigma_m))
        dj[[M]][[jm]] <- log(det(Sigma_m))
        uj[[M]][[jm]] <- t(Ylist[[M]][[jm]]-Xlist[[M]][[jm]]%*%rep(beta,dim(Xlist[[M]][[jm]])[1]))%*% Sigma_m_inv %*% (Ylist[[M]][[jm]]-Xlist[[M]][[jm]]%*%rep(beta,dim(Xlist[[M]][[jm]])[1]))
        Atilde[[M]][[jm]] <- list()
        omegatilde[[M]][[jm]] <- list()
        for(l in 1:M){
          Atilde[[M]][[jm]][[l]] <- list()
          omegatilde[[M]][[jm]][[l]] <- t(Xlist[[M]][[jm]]%*%Wlist[[M]][[l]][[jm]])%*%Sigma_m_inv%*%(Ylist[[M]][[jm]]-Xlist[[M]][[jm]]%*%rep(beta,dim(Xlist[[M]][[jm]])[1]))
          for(k in l:M){
            Atilde[[M]][[jm]][[l]][[k]] <- t(Xlist[[M]][[jm]]%*%Wlist[[M]][[k]][[jm]])%*%Sigma_m_inv%*%Xlist[[M]][[jm]]%*%Wlist[[M]][[l]][[jm]]
            #            omegatilde[[M]][[jm]][[l]][[k]] <- t(Xlist[[M]][[jm]]%*%Wlist[[M]][[k]][[jm]])%*%Sigma_m_inv%*%Ylist[[M]][[jm]]
          }
        }
      }
    }else{
      for(jm in 1:(dim(indicesW[[M]])[1])){
        A[[M]][[jm]] <- list()
        omega[[M]][[jm]] <- vector(mode = 'list',length = M)
        jmp1 <- indicesW[[M+1]] %>% left_join(indicesW[[M]]) %>%
          filter(.data[[paste0('iP',M)]]==jm) %>%
          dplyr::select(.data[[paste0('iP',M+1)]])
        for(l in 1:M){
          A[[M]][[jm]][[l]] <- vector(mode = 'list',length = M)
          omega[[M]][[jm]][[l]] <- 0
          for(jp in unlist(jmp1)){
            omega[[M]][[jm]][[l]] <- omega[[M]][[jm]][[l]]+omegatilde[[M+1]][[jp]][[l]]
          }
          for(k in l:M){
            A[[M]][[jm]][[l]][[k]] <- 0
            for(jp in unlist(jmp1)){
              A[[M]][[jm]][[l]][[k]] <- A[[M]][[jm]][[l]][[k]]+Atilde[[M+1]][[jp]][[l]][[k]]
            }
          }
        }
        Ktildeinv <- Wlist[[M]][[M]][[jm]]+A[[M]][[jm]][[M]][[M]]
        Ktilde <- chol2inv(chol(Ktildeinv))
        dj[[M]][[jm]] <- log(det(Ktildeinv))-log(det(Wlist[[M]][[M]][[jm]]))
        uj[[M]][[jm]] <- -t(omega[[M]][[jm]][[M]])%*%Ktilde%*%omega[[M]][[jm]][[M]]
        dtemp <- utemp <- 0
        for(jp in unlist(jmp1)){
          dtemp <- dtemp+dj[[M+1]][[jp]]
          utemp <- utemp+uj[[M+1]][[jp]]
        }
        dj[[M]][[jm]] <- dj[[M]][[jm]] + dtemp
        uj[[M]][[jm]] <- uj[[M]][[jm]] + utemp
        Atilde[[M]][[jm]] <- list()
        omegatilde[[M]][[jm]] <- vector(mode = 'list',length = M)
        for(l in 1:M){
          Atilde[[M]][[jm]][[l]] <- vector(mode = 'list',length = M)
          omegatilde[[M]][[jm]][[l]] <- omega[[M]][[jm]][[l]] - A[[M]][[jm]][[k]][[M]]%*%Ktilde%*%omega[[M]][[jm]][[M]]
          for(k in l:M){
            Atilde[[M]][[jm]][[l]][[k]] <- A[[M]][[jm]][[l]][[k]]-A[[M]][[jm]][[k]][[M]] %*% Ktilde %*% A[[M]][[jm]][[l]][[M]]
          }
        }
      }
    }
  }
    m2logv <- sum(unlist(dj))+sum(unlist(uj))
  return(m2logv)
}

likelihoodGaussian <- function(nu,range,sigma2,taue,betas,type){
  Sigma <- cExpMat(hh,hh,type,range,sigma2,nu)
  Y <- hh$Yresp
  X <- hh$Xcov
  XX <- diag(X)
  Sigmainv <- chol2inv(chol(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))
  m2logv <- log(det(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))+
    t(Y-betas*X)%*%Sigmainv%*%(Y-betas*X)
  return(m2logv)
}

likelihoodBanerjee <- function(nu,range,sigma2,taue,betas,type){
  #Sigma <- cExpMat(hh,hh,type,range,sigma2,nu)
  C <- cExpMat(knotsMRA[[1]],hh,type,range,sigma2,nu)
  Cstar <- cExpMat(knotsMRA[[1]],knotsMRA[[1]],type,range,sigma2,nu)
  Sigma <- t(C) %*% chol2inv(chol(Cstar)) %*% C 
  Y <- hh$Yresp
  X <- hh$Xcov
  XX <- diag(X)
  Sigmainv <- chol2inv(chol(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))
  m2logv <- log(det(XX%*%Sigma%*%t(XX)+(1/taue)*diag(dim(Sigma)[1])))+
    t(Y-betas*X)%*%Sigmainv%*%(Y-betas*X)
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


likelihoodFSA_Block <- function(nu,range,sigma2,taue,betas,type){
  Y <- hh$Yresp
  X <- hh$Xcov
  XX <- diag(X)
  
  C <- cExpMat(knotsMRA[[1]],hh,type,range,sigma2,nu)
  Cstar <- cExpMat(knotsMRA[[1]],knotsMRA[[1]],type,range,sigma2,nu)
  Sigmaw <- t(C) %*% chol2inv(chol(Cstar)) %*% C 
  Sigma <- cExpMat(hh,hh,type,range,sigma2,nu)
  Kappa <- Blockmatrix(hh$iK2)
  Sigmae <- (Sigma-Sigmaw)*Kappa
  XSigmae <- t(XX)%*%Sigmae%*%XX+sigma2*diag(dim(Sigmae)[1])
  library(Matrix)
  XSigmae <- as(XSigmae ,'dgCMatrix')
  blocks <- unique(ExtractBlocks(XSigmae))
  blocksinv <- purrr::map(blocks,~chol2inv(chol(.)))
  blocksdet <- purrr::map_dbl(blocks,~det(.))
  XSigmaeinv <- bdiag(blocksinv)
  SigmaYinv <- XSigmaeinv-XSigmaeinv%*%t(XX)%*%t(C)%*%solve(Cstar+C%*%XX%*%XSigmaeinv%*%XX%*%t(C))%*%C%*%XX%*%XSigmaeinv
  SigmaYdet <- det(Cstar+C%*%XX%*%XSigmaeinv%*%XX%*%t(C))*(det(Cstar))^(-1)*prod(blocksdet)
  return(as.numeric(log(SigmaYdet)+t(Y-betas*X)%*%SigmaYinv%*%(Y-betas*X)))
}
  