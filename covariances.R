library(pdist)
library(fields)

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


corrMaternduo_fields <- function(points_sf1,points_sf2,variance) {
  coords1 <- st_coordinates(points_sf1$geometry)
  coords2 <- st_coordinates(points_sf2$geometry)
  m <- matrix(0,nrow = dim(coords1)[1],ncol = dim(coords2)[1])
  #if(identical(coords1,coords2)){
  #  diag(m) <- variance
  #}
  m <- variance*stationary.cov(coords1,coords2,Covariance = 'Matern',
                               Distance = 'rdist.earth')
  return(m)
}

covKolmogorovHurst <- function(points_sf1,points_sf2,H,variance){
  coords1 <- st_coordinates(points_sf1$geometry)
  coords2 <- st_coordinates(points_sf2$geometry)
  m <- ifelse(identical(coords1,coords2)==TRUE,list(as.matrix(dist(coords1,diag = T,upper = T))),
              list(as.matrix(pdist(coords1,coords2))))
  m <- variance*m^(4*H-4)
}