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
      #      show(paste(m,jm,sep = '-'))
      Qlist[[M]][[jm]] <- knotsMRA[[M]] %>% 
        filter(.data[[paste0('iP',M)]]==jm)
      Xlist[[M]][[jm]] <- diag((regionalpoints %>% 
                                  filter(.data[[paste0('iP',M)]]==jm) %>%
                                  select(Xcov))$Xcov)
      Ylist[[M]][[jm]] <- (regionalpoints %>%
                             filter(.data[[paste0('iP',M)]]==jm) %>%
                             select(Yresp))$Yresp
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
        Wlist[[M]][[l]][[jm]] <- corrMaternduo(Qlist[[M]][[jm]],
                                               Qlist[[l]][[jl]],
                                               kappa,
                                               sigma2)-factorW
        #        Wlist[[M]][[l]][[jm]] <- corrMaternduo_fields(Qlist[[M]][[jm]],
        #                                               Qlist[[l]][[jl]],sigma2)-factorW
        
        #image.plot(Wlist[[M]][[l]][[jm]],legend.lab = paste(M,l,sep = '-'))
      }
    }
  }
  matrices_r <- list(W=Wlist,X=Xlist,Y=Ylist)
  return(matrices_r)
}


likelihoodKatzfuss <- function(kappa,sigma2,taue){
  
  matrices_r <- WQXYmakerMatern()
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
        uj[[M]][[jm]] <- t(Ylist[[M]][[jm]])%*% Sigma_m_inv %*% Ylist[[M]][[jm]]
        Atilde[[M]][[jm]] <- list()
        omegatilde[[M]][[jm]] <- list()
        for(l in 1:M){
          Atilde[[M]][[jm]][[l]] <- list()
          omegatilde[[M]][[jm]][[l]] <- t(Xlist[[M]][[jm]]%*%Wlist[[M]][[l]][[jm]])%*%Sigma_m_inv%*%Ylist[[M]][[jm]]
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
          select(.data[[paste0('iP',M+1)]])
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