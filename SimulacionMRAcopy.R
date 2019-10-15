library(fields)
library(dplyr)
library(sf)
library(rgeos)
library(plotly)
library(lwgeom)
library(pdist)
library(stringr)

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
nc_n <- c(1,2,2,3,3) #c(2,2,3,3)
nr_n <- c(1,2,11,1,1) #c(2,11,1,1)



indicesglob  <- list()
indicesreg   <- list()
cellloc.glob <- list()
cellloc.reg  <- list()
indicesglobK  <- list()
indicesregK   <- list()

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

indicesglobtemp <- as.data.frame(cellloc.glob[[i]]) %>% 
  mutate(rown = row%%nr_n[i],coln=col%%nc_n[i]) %>% 
  mutate(rown=ifelse(rown==0,nr_n[i],rown),
         coln=ifelse(coln==0,nc_n[i],coln))

indicesregtemp <- as.data.frame(cellloc.reg[[i]]) %>% 
  mutate(rown = row%%nr_n[i],coln=col%%nc_n[i]) %>% 
  mutate(rown=ifelse(rown==0,nr_n[i],rown),
         coln=ifelse(coln==0,nc_n[i],coln))

indexmatrix <- as.data.frame(expand.grid(1:nr_n[i],1:nc_n[i]))
indexmatrix <- indexmatrix %>% dplyr::select(rown=Var1,coln=Var2)%>%
  mutate(celln=1:(nr_n[i]*nc_n[i]))

indicesglobK[[i]] <- as.numeric((indicesglobtemp %>% left_join(indexmatrix,by = c('rown','coln')) %>%
  dplyr::select(celln))$celln)
indicesregK[[i]] <- as.numeric((indicesregtemp %>% left_join(indexmatrix,by = c('rown','coln')) %>%
  dplyr::select(celln))$celln)


#plot(globraster)
#plot(rasterToPolygons(globraster),add=T)
}



# table for regional indices:
indicesregK   <- data.frame(matrix(unlist(indicesregK), ncol=nn,
                        nrow=length(indicesregK[[1]])))
names(indicesregK) <- as.character(unlist(lapply(1:nn,
                                  function(i)paste0("iP",i))))

all <- cbind(dataset, indicesregK)

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

# table for global indices: (or eventually knots)
datasetglob <- dataset$coo %>% dplyr::select(IDglo,latglo,longlo) %>% 
  distinct(IDglo,latglo,longlo)

indicesglobK   <- data.frame(matrix(unlist(indicesglobK), ncol=nn,
                                  nrow=length(indicesglobK[[1]])))
names(indicesglobK) <- as.character(unlist(lapply(1:nn,
                                                function(i)paste0("iP",i))))

allglob <- cbind(datasetglob, indicesglobK)


globalpoints <- SpatialPointsDataFrame(cbind(allglob$longlo, 
                                             allglob$latglo), 
                                         as.data.frame(cbind(
                                           iP1   = allglob$iP1,
                                           iP2   = allglob$iP2,
                                           iP3   = allglob$iP3,
                                           iP4   = allglob$iP4,
                                           iP5   = allglob$iP5)),
                                         proj4string = CRS('+proj=longlat +datum=WGS84'),
                                         bbox = bordes)

globalpoints = st_as_sf(globalpoints)
plot(globalpoints)
plot(globalpoints %>% filter(iP3==1) )

## Aquí es donde se puede hacer la aplicación para cada nivel y partición

library(purrr)
regionalpoints %>%
  split(.$iP4) %>%
  map(summary) 


##Ejercicio calculo de W según formula (6) de Katzfuss

Q5a <- globalpoints %>% filter(iP4==1)
Q5b <- globalpoints %>% filter(iP1==1)

corrMatern <- function(points_sf,kappa, variance, nu=1) {
  coords <- st_coordinates(points_sf$geometry)
  m <- as.matrix(dist(coords))
  m <- variance*exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
  m[is.nan(m)] <- variance
  #diag(m) <- variance
  return(m)
}
corrMatern(Q5a,kappa,sigma2)

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
MM <- corrMaternduo(Q5a,Q5b,kappa,sigma2)


Wmaker <- function(globalpoints,M,nc,nr){
  Wlist <- list()
  for(m in 0:M){
    Qlist <- list()
    Jm <- nc[[m+1]]*nr[[m+1]]
    for(j in 1:Jm){
     Qlist[[j]] <- globalpoints %>% filter(iP1==1)  
    }
  }
}



globalpoints <- globalpoints %>% mutate(iG1=as.numeric(paste0(iP1)),iG2=as.numeric(paste0(iP1,iP2)),
                               iG3=as.numeric(paste0(iP1,iP2,iP3)),
                               iG4=as.numeric(paste0(iP1,iP2,iP3,iP4)),
                               iG4=as.numeric(paste0(iP1,iP2,iP3,iP4,iP5))) 



indicesW <- list(sort(unique(globalpoints$iG1)),
                 sort(unique(globalpoints$iG2)),
                 sort(unique(globalpoints$iG3)),
                 sort(unique(globalpoints$iG4)),
                 sort(unique(globalpoints$iG5)))




Qlist <- list()
Wlist <- list()
#m=0 j0=1 l=0
Q0 <- globalpoints %>% filter(iG1==1)
Qlist[[1]] <- list()
Qlist[[1]][[1]] <- Q0
Wlist[[1]] <- list()
W0 <- corrMatern(Q0,kappa,sigma2)
Wlist[[1]][[1]] <- W0

#m=1 j1=1,..,4 l=0,1
Qlist[[2]] <- list()
Wlist[[2]] <- list()
Wlist[[2]][[1]] <- list()
Wlist[[2]][[2]] <- list()
for(j1 in 1:length(indicesW[[2]])){
  Qlist[[2]][[j1]] <- globalpoints %>% filter(iG2==indicesW[[2]][j1])
  Wlist[[2]][[1]][[j1]] <- corrMaternduo(Qlist[[2]][[j1]],Qlist[[1]][[1]],kappa,sigma2)
  Wlist[[2]][[2]][[j1]] <- corrMaternduo(Qlist[[2]][[j1]],Qlist[[2]][[j1]],kappa,sigma2)-
    Wlist[[2]][[1]][[j1]]%*%solve(Wlist[[1]][[1]])%*%t(Wlist[[2]][[1]][[j1]])
}

#m=2 j2=1,..,11 l=0,1,2
Qlist[[3]] <- list()
Wlist[[3]] <- list()
Wlist[[3]][[1]] <- list()
Wlist[[3]][[2]] <- list()
Wlist[[3]][[3]] <- list()
for(j2 in 1:length(indicesW[[3]])){
  Qlist[[3]][[j2]] <- globalpoints %>% filter(iG3==indicesW[[3]][j2])
  Wlist[[3]][[1]][[j2]] <- corrMaternduo(Qlist[[3]][[j2]],Qlist[[1]][[1]],kappa,sigma2)
  #Qlist[[3]][[j2]]$iG2[1]
  #Wlist[[3]][[2]][[j2]] <- corrMaternduo(Qlist[[3]][[j2]],Qlist[[2]][[j1]],kappa,sigma2)-
  #  Wlist[[2]][[1]][[j2]]%*%solve(Wlist[[1]][[1]])%*%t(Wlist[[2]][[1]][[j2]])
}



Q01_1 <- globalpoints %>% filter(iP1==1,iP2==1,iP3==1)
Q01_2 <- globalpoints %>% filter(iP1==1,iP2==1,iP3==2)
Q01_3 <- globalpoints %>% filter(iP1==1,iP2==1,iP3==3)
Q01_4 <- globalpoints %>% filter(iP1==1,iP2==1,iP3==4)
Q01_5 <- globalpoints %>% filter(iP1==1,iP2==1,iP3==5)


W01_1_0 <- corrMaternduo(Q01_1,Q0,kappa,sigma2)
W01_2_0 <- corrMaternduo(Q01_2,Q0,kappa,sigma2)
W01_3_0 <- corrMaternduo(Q01_3,Q0,kappa,sigma2)
W01_4_0 <- corrMaternduo(Q01_4,Q0,kappa,sigma2)
W01_5_0 <- corrMaternduo(Q01_5,Q0,kappa,sigma2)

W01_1_1 <- corrMaternduo(Q01_1,Q01,kappa,sigma2)-W01_1_0%*%solve(W0)%*%t(W01_0)
W01_2_1 <- corrMaternduo(Q01_2,Q01,kappa,sigma2)-W01_2_0%*%solve(W0)%*%t(W01_0)
W01_3_1 <- corrMaternduo(Q01_3,Q01,kappa,sigma2)-W01_3_0%*%solve(W0)%*%t(W01_0)
W01_4_1 <- corrMaternduo(Q01_4,Q01,kappa,sigma2)-W01_4_0%*%solve(W0)%*%t(W01_0)
W01_5_1 <- corrMaternduo(Q01_5,Q01,kappa,sigma2)-W01_5_0%*%solve(W0)%*%t(W01_0)

W01_1_2 <- corrMaternduo(Q01_1,Q01_1,kappa,sigma2)-W01_1_0%*%solve(W0)%*%t(W01_1_0)-W01_1_1%*%solve(W01_1)%*%t(W01_1_1)
W01_2_2 <- corrMaternduo(Q01_2,Q01_2,kappa,sigma2)-W01_2_0%*%solve(W0)%*%t(W01_2_0)-W01_2_1%*%solve(W01_1)%*%t(W01_2_1)
W01_3_2 <- corrMaternduo(Q01_3,Q01_3,kappa,sigma2)-W01_3_0%*%solve(W0)%*%t(W01_3_0)-W01_3_1%*%solve(W01_1)%*%t(W01_3_1)
W01_4_2 <- corrMaternduo(Q01_4,Q01_4,kappa,sigma2)-W01_4_0%*%solve(W0)%*%t(W01_4_0)-W01_4_1%*%solve(W01_1)%*%t(W01_4_1)
W01_5_2 <- corrMaternduo(Q01_5,Q01_5,kappa,sigma2)-W01_5_0%*%solve(W0)%*%t(W01_5_0)-W01_5_1%*%solve(W01_1)%*%t(W01_5_1)
