library(fields)
library(dplyr)
library(sf)
library(rgeos)
library(plotly)

set.seed(1)

load("datos/resolucion/coordinates.Rdata")
load("datos/resolucion/domainfinal.Rdata")

## Cut a window:
# 230 a 300 y de 30 a 60

global<-final[final$lonval>240 &final$lonval<290
              & final$latval>30 & final$latval<60,c("lonval","latval")]
regional<-t(aa)[t(aa)[,1]>30&t(aa)[,1]<60
                &t(aa)[,2]>240&t(aa)[,2]<290,2:1]
plot(global, cex=0.01)
points(regional,cex=0.01,col="blue")

N <- dim(global)[1] #Number of spatial points global
n <- dim(regional)[1] #Number of spatial points regional
k <- 1 #Observations through time 

beta0 <- 0
beta1 <- 1
kappa <- 1.5
sigma2 <- 1/4
ncov <- 1
rho <- c(0.7, 0.5) 
rMatern <- function(n, coords, kappa, variance, nu=1) {
  m <- as.matrix(dist(coords))
  m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
  diag(m) <- 1
  return(drop(crossprod(chol(variance*m),
                        matrix(rnorm(nrow(coords)*n), ncol=n))))
}
beta1s <- rMatern(ncov, regional, kappa, sigma2)

## simulate the covariate values (equal values for same quadrant)
hh<-SpatialPointsDataFrame(global, data.frame(ID=1:dim(global)[1],
                                              y=runif(N*k)))

regionalpoints <- SpatialPoints(regional)
globalpoints <- SpatialPoints(global)
proj4string(regionalpoints) <- '+proj=longlat +datum=WGS84'
proj4string(globalpoints) <- '+proj=longlat +datum=WGS84'
proj4string(hh) <- '+proj=longlat +datum=WGS84'


dev.off()
plot(global, cex=0.01)
points(regional,cex=0.01,col="blue")
hh2<-gBuffer(hh, width=0.75, quadsegs=1, 
             capStyle='SQUARE', byid=T)

plot(hh2, add=TRUE)
points(regionalpoints, cex=0.1)
points(globalpoints, cex=0.5, col="blue")

## match regional with global
aa <- regionalpoints%over%geometry(hh2)
coo<-cbind(regional,global[unlist(aa),], unlist(aa))
plot(coo[,1:2], cex=0.01,col=coo[,5])
points(coo[,3:4],cex=0.01)
names(coo)<-c("lonreg", "latreg", "longlo", "latglo", "IDglo")

taue <- 20 ##Precision parameter for error
error <- rnorm(n*k, 0, sqrt(1/taue)) ### error in the observation

hh4<-hh[coo[,5],2]$y

y <- beta0+(beta1+beta1s)*hh4+error ##Simulate the observations



