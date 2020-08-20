library(fields)
library(dplyr)
library(sf)
library(rgeos)
library(plotly)
library(lwgeom)
library(pdist)
library(stringr)
library(tidyr)
library(rlang)

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

hh2<-gBuffer(hh, width=0.71, quadsegs=1, capStyle='SQUARE', byid=T)
plot(hh2, add=TRUE,col="red")
points(regionalpoints, cex=0.1)
points(globalpoints, cex=0.5, col="blue")

## match regional with global
aa <- over(regionalpoints, geometry(hh2))
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


corrMatern <- function(points_sf,kappa, variance, nu=1) {
  coords <- st_coordinates(points_sf$geometry)
  m <- as.matrix(dist(coords))
  m <- variance*exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
  m[is.nan(m)] <- variance
  #diag(m) <- variance
  return(m)
}
#corrMatern(Q5a,kappa,sigma2)

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


indicesW <- list()
indicesW[[1]] <- globalpoints %>% expand(iP1) %>% arrange(iP1)%>%
  mutate(indexg1 = 1:n())
indicesW[[2]] <- globalpoints %>% expand(iP1,iP2) %>% 
  arrange(iP1,iP2) %>% mutate(indexg2 = 1:n())
indicesW[[3]] <- globalpoints %>% expand(iP1,iP2,iP3) %>%
  arrange(iP1,iP2,iP3) %>% mutate(indexg3 = 1:n())
indicesW[[4]] <- globalpoints %>% expand(iP1,iP2,iP3,iP4) %>%
  arrange(iP1,iP2,iP3,iP4) %>% mutate(indexg4 = 1:n())
indicesW[[5]] <- globalpoints %>% expand(iP1,iP2,iP3,iP4,iP5) %>%
  arrange(iP1,iP2,iP3,iP4,iP5) %>% mutate(indexg5 = 1:n())

tablaindicesW <- indicesW[[5]] %>% left_join(indicesW[[4]])%>%
  left_join(indicesW[[3]]) %>% left_join(indicesW[[2]]) %>%
  left_join(indicesW[[1]])


globalpoints <- globalpoints %>% left_join(tablaindicesW) 
regionalpoints <- regionalpoints %>% left_join(tablaindicesW) 


##Matrices W according to Katzfuss, 2017 (need to be optimized!)
##Indexing: [[m]][[l]][[partition number]]


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
      Qlist[[M]][[jm]] <- globalpoints %>% filter(.data[[paste0('indexg',M)]]==jm)
      indicesjerarq <- Qlist[[M]][[jm]] %>% dplyr::select(starts_with('indexg'))%>%
        st_drop_geometry()
      for(l in 1:M){
        jl <- as.numeric(indicesjerarq %>% dplyr::select(.data[[paste0('indexg',l)]]) %>% unique())
        factorW <- 0
        if(l!=1){
          factorW <- 0
          for(k in 1:(l-1)){
            jk <- as.numeric(indicesjerarq %>% dplyr::select(.data[[paste0('indexg',k)]]) %>% unique())     
            factorW <- factorW + Wlist[[M]][[k]][[jm]]%*%solve(Wlist[[k]][[k]][[jk]])%*%t(Wlist[[l]][[k]][[jl]])
          }
        }  
        Wlist[[M]][[l]][[jm]] <- corrMaternduo(Qlist[[M]][[jm]],Qlist[[l]][[jl]],kappa,sigma2)-factorW  
      }
    }
  }
  return(Wlist)
}


#save(Qlist,Wlist,file = 'WQmatrices.RData')
#load('WQmatrices.RData')

###Likelihood calculation (fixed parameters)
