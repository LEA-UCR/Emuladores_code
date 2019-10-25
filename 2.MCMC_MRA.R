library(pdist)
source("1.MRA_resolution.R")

corrMaternduo <- function(points_sf1,points_sf2,kappa, variance, nu=1) {
  coords1 <- st_coordinates(points_sf1$geometry)
  coords2 <- st_coordinates(points_sf2$geometry)
  m <- as.matrix(pdist(coords1,coords2))
  m <- variance*exp((1-nu)*log(2) + nu*log(kappa*m)-
                      lgamma(nu))*besselK(m*kappa, nu)
  m[is.nan(m)] <- variance
  #diag(m) <- variance
  return(m)
}

tau2 <- 1

WQmaker <- function(){
  Qlist <- list()
  Wlist <- list()
  for(m in 0:(nn-1)){
    M <- m+1
    Qlist[[M]] <- list()
    Wlist[[M]] <- list()
    for(l in 1:M){
      Wlist[[M]][[l]] <- list()
    }
    for(jm in 1:(dim(indicesW[[M]])[1])){
      show(paste(m,jm,sep = '-'))
      Qlist[[M]][[jm]] <- knotsMRA[[M]] %>% filter(.data[[paste0('iP',M)]]==jm)
      indicesjerarq <- Qlist[[M]][[jm]] %>% dplyr::select(starts_with('iP'))%>%
        st_drop_geometry()
      for(l in 1:M){
        jl <- as.numeric(indicesjerarq %>% dplyr::select(.data[[paste0('iP',l)]]) %>% unique())
        factorW <- 0
        if(l!=1){
          factorW <- 0
          for(k in 1:(l-1)){
            jk <- as.numeric(indicesjerarq %>% dplyr::select(.data[[paste0('iP',k)]]) %>% unique())     
            factorW <- factorW + Wlist[[M]][[k]][[jm]]%*%solve(Wlist[[k]][[k]][[jk]])%*%t(Wlist[[l]][[k]][[jl]])
          }
        }
        Wlist[[M]][[l]][[jm]] <- corrMaternduo(Qlist[[M]][[jm]],Qlist[[l]][[jl]],kappa,sigma2)-factorW  
      }
    }
  }
  return(Wlist)
}

Wlist <- WQmaker()
