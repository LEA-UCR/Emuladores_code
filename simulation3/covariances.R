library(pdist)
library(fields)
library(RandomFields)

corrMaternduo <- function(points_sf1,points_sf2,kappa, variance, nu=1) {
  coords1 <- st_coordinates(points_sf1$geometry)
  coords2 <- st_coordinates(points_sf2$geometry)
  m <- ifelse(identical(coords1,coords2)==TRUE,list(as.matrix(dist(coords1,diag = T,upper = T))),
              list(as.matrix(pdist(coords1,coords2))))
  m <- variance*exp((1-nu)*log(2) + nu*log(kappa*m[[1]])-
                      lgamma(nu))*besselK(m[[1]]*kappa, nu)
  m[is.nan(m)] <- variance
  #diag(m) <- variance
  return(m)
}

cExpMat <- function(points_sf1,points_sf2,type,range=1,
                    variance,nu=1,acau,bcau){
  coords1 <- st_coordinates(points_sf1$geometry)
  coords2 <- st_coordinates(points_sf2$geometry)
  if(type=='Exponential'){
    m <- stationary.cov(coords1,coords2,Covariance = type,
                        Distance = 'rdist.earth',
                        theta = range,phi=variance)
  }
  if(type=='Matern'){
    m <- stationary.cov(coords1,coords2,Covariance = type,
                        Distance = 'rdist.earth',
                        theta = range,phi=variance,nu=nu)
  }
  if(type=='Cauchy'){
    m <- rdist.earth(coords1,coords2)
    m <- (1+abs(m)^acau)^(-bcau/acau) #Gneiting & Schlather, 2004
  }
  return(m)
}

covKolmogorovHurst <- function(points_sf1,points_sf2,H,variance){
  coords1 <- st_coordinates(points_sf1$geometry)
  coords2 <- st_coordinates(points_sf2$geometry)
  m <- ifelse(identical(coords1,coords2)==TRUE,
              list(as.matrix(dist(coords1,diag = T,upper = T))),
              list(as.matrix(pdist(coords1,coords2))))
  m <- variance*m^(4*H-4)
}