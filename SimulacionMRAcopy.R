library(fields)
library(dplyr)
library(sf)
library(rgeos)
library(plotly)

set.seed(1)

load("/../home/Emuladores/datos/resolucion/coordinates.Rdata")
load("/../home/Emuladores/datos/resolucion/domainfinal.Rdata")

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
hh2<-gBuffer(hh, width=0.7, quadsegs=1, capStyle='SQUARE', byid=T)
plot(hh2, add=TRUE,col="red")
points(regional, cex=0.1)
points(global, cex=0.5, col="blue")

## match regional with global
aa <- over(SpatialPoints(regional), geometry(hh2), 
                                 returnList = TRUE)
coo<-cbind(regional,global[unlist(aa),], unlist(aa))
plot(coo[,1:2], cex=0.01,col=coo[,5])
points(coo[,3:4],cex=0.01)
names(coo)<-c("lonreg", "latreg", "longlo", "latglo", "IDglo")

taue <- 20 ##Precision parameter for error
error <- rnorm(n*k, 0, sqrt(1/taue)) ### error in the observation

hh4<-hh2[coo[,5],2]$y

y <- beta0+(beta1+beta1s)*hh4+error ##Simulate the observations

### y is the response variable with regional resolution
### hh4 are     the xs with global resolution

dataset<-tibble(coo,resp=y,covariate=hh4)

ggplot(data = dataset, mapping = aes(x = coo$lonreg,
                                     y = coo$latreg,
                                     color=resp)) + 
  geom_tile(size=2, alpha=0.4) + theme_bw() +
  scale_color_continuous(low="blue", high="yellow") +
  ylab("latitude") + xlab("longitude") 

ggplot(data = dataset, mapping = aes(x = coo$lonreg,
                                     y = coo$latreg,
                                     color=covariate)) + 
  geom_tile(size=2, alpha=0.4) + theme_bw() +
  scale_color_continuous(low="blue", high="yellow") +
  ylab("latitude") + xlab("longitude") 

plot(dataset$covariate,dataset$resp, ylab="Response", xlab="Covariate")

### Use MRA to analyze this dataset:

source('MRAfunctions.R')
library(ltsa) # for fast simulation of truth
library(fields) # needed for rdist function

### covariance function
theta=NA
cov.fun=function(locs1,locs2,theta=NA) {
  if( (length(locs1)==0) || (length(locs2)==0) ){ 0 } else { 
    x=rdist(locs1,locs2)/.05
    exp(-x) + (abs(x)<1e-8)*.05
    (1 + sqrt(3)*x)*exp(-sqrt(3)*x)  # Matern with smoothness 1.5
  }
}

# partition 2D !!!! (simulation2D.ji)

# cov.mats=vector("list",M+1)
#for(m in 0:M){
#  cov.mats[[m+1]]=MRA.illus(NA,cov.fun,part.illus$data,part.illus$knots,part.illus$indices,M.cov.plot=m)[[2]]
#}
#u=10; xgrid=(1:(u*n.all))/(n.all*u); part.illus=partition.1D(J,M,r,domain,xgrid,1:(n.all*u))
#bfs=MRA.illus(NA,cov.fun,part.illus$data,part.illus$knots,part.illus$indices)[[1]]


