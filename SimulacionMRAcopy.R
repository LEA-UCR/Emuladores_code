library(fields)
library(dplyr)
library(sf)
library(rgeos)
library(plotly)

set.seed(1)

load("../resolucion/coordinates.Rdata")
load("../resolucion/domainfinal.Rdata")

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
proj4string(hh) <- '+proj=longlat +datum=WGS84'
regionalpoints <- SpatialPoints(regional)
globalpoints <- SpatialPoints(global)
proj4string(regionalpoints) <- '+proj=longlat +datum=WGS84'
proj4string(globalpoints) <- '+proj=longlat +datum=WGS84'

hh2<-gBuffer(hh, width=0.7, quadsegs=1, capStyle='SQUARE', byid=T)
plot(hh2, add=TRUE,col="red")
points(regionalpoints, cex=0.1)
points(globalpoints, cex=0.5, col="blue")

## match regional with global
aa <- over(regionalpoints, geometry(hh2), 
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
#plot(dataset$covariate,dataset$resp, ylab="Response", xlab="Covariate")

### Use MRA to analyze this dataset:
#source('MRAfunctions.R')
library(sf)
library(raster)

# create a raster with the points and a data table:

bordes <- bbox(hh2)
crsglobal <- CRS('+proj=longlat +datum=WGS84')

### Aquí pueden cambiarse las particiones dependiendo del problema.
partitions <- list(1,c(1:4),seq(1,2*2*2*11), 
                   seq(1,2*2*3*2*11), 
                   seq(1,2*2*3*3*2*11))
length(partitions) # how many levels
nc <- list(1, 2, 2*2,  2*2*3, 2*2*3*3)
nr <- list(1, 2, 2*11, 2*11,  2*11)

indicesglob  <- list()
indicesreg   <- list()
cellloc.glob <- list()
cellloc.reg  <- list()

# ¿Cuántas veces queremos partir el dominio y cuál es el borde?
nn <- length(nc)
temp <- raster(xmn=bordes[1],ymn=bordes[2],xmx=bordes[3],
               ymx=bordes[4],val=1,
               crs=crsglobal,ncols=1,nrows=1)

for(i in 1:nn){
globraster <- raster(xmn=bordes[1],ymn=bordes[2],xmx=bordes[3],
                       ymx=bordes[4],val=partitions[[i]],
                       crs=crsglobal,ncols=nc[[i]],nrows=nr[[i]])
  
indicesglob[[i]]    <- extract(globraster,globalpoints, cellnumbers=TRUE)[,1]
cellloc.glob[[i]]   <- rowColFromCell(globraster,indicesglob[[i]])
indicesreg[[i]]     <- extract(globraster,regionalpoints, cellnumbers=TRUE)[,1]
cellloc.reg[[i]]    <- rowColFromCell(globraster,indicesreg[[i]])

plot(globraster)
plot(rasterToPolygons(globraster),add=T)
}

# table for regional indices:
indicesreg   <- data.frame(matrix(unlist(indicesreg), ncol=nn,
                        nrow=length(indicesreg[[1]])))
names(indicesreg) <- as.character(unlist(lapply(1:nn,
                                  function(i)paste0("iP",i))))

all <- cbind(dataset, indicesreg)

str(all)

cbind(all$coo$IDglo, all$iP5)


regionalpoints <- SpatialPointsDataFrame(cbind(all$coo$lonreg, 
                                               all$coo$latreg), 
                           as.data.frame(cbind(
                             longlo= all$coo$longlo,
                             latglo= all$coo$latglo,
                             IDglo = all$coo$IDglo,
                             Yresp = all$resp,
                             Xcov  = all$covariate,
                             iP1   = all$iP1,
                             iP2   = all$iP2,
                             iP3   = all$iP3,
                             iP4   = all$iP4,
                             iP5   = all$iP5)),
                  proj4string = CRS('+proj=longlat +datum=WGS84'),
                  bbox = bordes)

regionalpoints = st_as_sf(regionalpoints)
plot(regionalpoints)
plot(regionalpoints %>% filter(iP3==1) )

## Aquí es donde se puede hacer la aplicación para cada nivel y partición

library(purrr)
regionalpoints %>%
  split(.$iP4) %>%
  map(summary) 




