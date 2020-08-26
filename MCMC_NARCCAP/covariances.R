library(pdist)
library(fields)
library(RandomFields)
library(furrr)

corrMaternduo <- function(points_sf1,points_sf2,
                          kappa, variance, nu=1) {
  coords1 <- st_coordinates(points_sf1$geometry)
  coords2 <- st_coordinates(points_sf2$geometry)
  m <- ifelse(identical(coords1,coords2)==TRUE,
              list(as.matrix(dist(coords1,diag = T,upper = T))),
              list(as.matrix(pdist(coords1,coords2))))
  m <- variance*exp((1-nu)*log(2) + nu*log(kappa*m[[1]])-
                      lgamma(nu))*besselK(m[[1]]*kappa, nu)
  m[is.nan(m)] <- variance
  #diag(m) <- variance
  return(m)
}

# Covariance function between sf objects.
cExpMat <- function(points_sf1,points_sf2,type,range=1,
                    variance=1,nu=1,alpha=1,beta=1){
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
    m <- (1+abs(m)^alpha)^(-beta/alpha) #Gneiting & Schlather, 2004
  }
  return(m)
}

#Augmented Covariance function (Co-regionalization model)
cExpMat_mult <- function(points_sf1,points_sf2,A,nCov=1,typesCov,range,
                         nu,alpha,beta){
  matricesCov <- NULL
  Fun_aug <- function(i){
    matricesCov[[i]] <- cExpMat(points_sf1,points_sf2,typesCov[i],range,
                                variance=1,nu,alpha,beta)
    return(as.matrix(kronecker(matricesCov[[i]],
           Matrix(A[,i]%*%t(A[,i]),sparse = T))))
  }
  #plan(multiprocess)
  #matrices_aug <- future_map(1:nCov,~Fun_aug(.))
  matrices_aug <- purrr::map(1:nCov,~Fun_aug(.))
  
  return(Reduce("+",matrices_aug))
}