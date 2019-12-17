library(fields)
library(dplyr)
library(sf)
library(rgeos)
library(plotly)

set.seed(1)

#Grid size and definition
nlat <- 30
nlon <- 30
gridbase <- expand.grid(lon=seq(230,300,length.out = nlon),
                        lat=seq(30,60,length.out = nlat))

gridlist <- list(x=seq(230,300,length.out = nlon),
                 y=seq(30,60,length.out = nlat))


n <- dim(gridbase)[1] #Number of spatial points
k <- 1 #Observations through time 

#Non-spatial parameters
beta0 <- 0
beta1 <- 1

#Random field generator
rExpMat <- function(n,coords,type,range,variance,nu=1){
  if(type=='Exponential'){
    m <- stationary.cov(coords,Covariance = type,
                        Distance = 'rdist.earth',
                        theta = range,phi=variance)
  }
  if(type=='Matern'){
    m <- stationary.cov(coords,Covariance = type,
                        Distance = 'rdist.earth',
                        theta = range,phi=variance,nu=nu)
  }
  return(drop(crossprod(chol(m),
                        matrix(rnorm(nrow(coords)*n), ncol=n))))
}

#Number of covariates
ncov <- 1

#Covariance structure parameters
nu <- 1.5
range <- 4
sigma2 <- 1
type='Exponential'

#Spatial parameter generation
beta1s <- rExpMat(ncov,gridbase,type = type,
                  range = range,variance = sigma2)

#Nugget effect generation
taue <- 20 ##Precision parameter for error
error <- rnorm(n*k, 0, sqrt(1/taue)) ### error in the observation

#Covariates and dependent variable 
X <- runif(n*k)
y <- beta0+(beta1+beta1s)*X+error ##Simulate the observations

dataset<-tibble(Yresp=y,Xcov=X) 


#Spatial structure 
hh <- SpatialPointsDataFrame(coords = gridbase,data = dataset)
proj4string(hh) <- '+proj=longlat +datum=WGS84'

save(dataset, hh, file="dataset_general.Rdata")
