library(fields)
library(dplyr)
library(sf)
library(rgeos)
library(plotly)
library(RandomFields)

# model <- "SVC"; type <- "Exponential"; i<-1

function.to.gen.data <- function(model, type, i){
set.seed(1*i) 
### ############################ ###
###  REAL VALUES FOR PARAMETERS  ###
### ############################ ###
#Non-spatial parameters
beta0 <- 0 # intercept
beta1 <- 2 # beta 1 (without spatial structure)
#Number of covariates
ncov <- 1
taue <- 1 # Precision parameter for error
#Covariance structure parameters
nu <- 1  # (fixed)
range <- 0.9 # range parameter
taub <- 1 # 
sigma2 = 1/taub
### ############################ ###
###   Grid size and definition   ###
### ############################ ###
nlat <- 100
nlon <- 100
gridbase <- expand.grid(lon=seq(0,1,length.out = nlon),
                        lat=seq(0,1,length.out = nlat))
gridlist <- list(x=seq(0,1,length.out = nlon),
                 y=seq(0,1,length.out = nlat))
n <- dim(gridbase)[1] #Number of spatial points
k <- 1 #Observations through time 
#Random field generator
rExpMat <- function(n,coords,type,range,variance,nu=1.5){
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
### ############################ ###
### Spatial parameter generation ###
### ############################ ###
beta1s <- rExpMat(ncov,gridbase,type = type,
                 range = range,variance = sigma2)
#model_betas <- RMexp(var = sigma2,scale = range)
#beta1s <- RFsimulate(model_betas,x = gridbase$lon,y = gridbase$lat)
#beta1s <- beta1s$variable1
#Nugget effect generation
error <- rnorm(n*k, 0, sqrt(1/taue)) ### nugget
### ############################ ###
###  Cov and dependent variable  ###
### ############################ ###
  if(model=='SVI'){
    X <- rep(1,n*k)
    # Model SVI
    y <- beta0+(beta1)*X + (beta1s)*X + error 
  }
  if(model=='SVC'){
    X <- rgamma(n*k,2,1)
    # Model beta(s)
    y <- beta0+(beta1)*scale(X) + (beta1s)*X + error 
  }
  dataset<-tibble(Y=as.vector(y),X=X) 
### ############################ ###
###        Save data sets        ###
### ############################ ###
#Spatial structure 
hh <- SpatialPointsDataFrame(coords = gridbase,data = dataset)
proj4string(hh) <- '+proj=longlat +datum=WGS84'
save(dataset, hh, file=paste0("sim_data/dataset",
                              model,type,i,".Rdata"))
}

nsimulations <- 1

lapply(1:nsimulations,function(i)function.to.gen.data("SVC","Exponential",i))
#lapply(1:nsimulations,function(i)function.to.gen.data("SVI","Exponential",i))
#lapply(1:nsimulations,function(i)function.to.gen.data("SVC","Matern",i))
#lapply(1:nsimulations,function(i)function.to.gen.data("SVI","Matern",i))
