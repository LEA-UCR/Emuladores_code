library(raster) 
library(tidyverse)
library(sf)
library(stars)
library(gtools)
library(rlang)

gen_resolution <- function(dataset_file,nCov_v){
  load(datasetfile) # data loading
  hh<- hh[coordinates(hh)[,"lon"]>=240 &
          coordinates(hh)[,"lon"]<=265 &
          coordinates(hh)[,"lat"]>=20  &
          coordinates(hh)[,"lat"]<=40,]
# hh <- hh[hh$lat>=20 & hh$lat<=40 & hh$lon>=240 & hh$lon<=265,]
  bordes <- bbox(hh)
  hh_sf <- st_as_sf(hh)
  crsglobal <- CRS('+proj=longlat +datum=WGS84')
  
  
  nlevelsMRA <- 3 #Number of MRA levels (M)
  npartitions <- 4 #Number of partitions per level (J)
  nknots <- 10 #Number of knots per partition
  npartitions_r <- 2 
  npartitions_c <- npartitions/npartitions_r
  partitions <- list()
  nc <- nr <- NULL
  nc_n <- nr_n <- NULL
  
  for(i in 0:nlevelsMRA){ 
    # Sequence of partitions per level
    partitions[[i+1]] <- seq(1,npartitions^i)
    nc[i+1] <- npartitions_c^i
    nr[i+1] <- npartitions_r^i
  }
  
  nc_n <- c(1,rep(npartitions_c,nlevelsMRA-1))
  nr_n <- c(1,rep(npartitions_r,nlevelsMRA-1))
  nn <- length(nc)
  
  globraster <- list()
  Qlist_sf <- list()
  Qlist_sp <- list()
  for(i in 0:nlevelsMRA){ 
    #Qlist_sf: Random knots per partition block in sf format
    globraster[[i+1]] <- raster(xmn=bordes[1],ymn=bordes[2],xmx=bordes[3],
                         ymx=bordes[4],val=partitions[[i+1]],
                         crs=crsglobal,ncols=nc[i+1],nrows=nr[i+1])
    
    globraster_stars <- st_as_stars(globraster[[i+1]])
    globraster_sf <- st_as_sf(globraster_stars)
    
    Qlist_sf[[i+1]] <-  suppressMessages(do.call(
      rbind,purrr::map(.x = partitions[[i+1]],
      function(l) st_sf(st_sample(globraster_sf %>% 
      filter(layer==l),size = nknots)))))
    Qlist_sp[[i+1]] <- suppressMessages(as_Spatial(Qlist_sf[[i+1]]))
  }
  
  for(i in 0:nlevelsMRA){
    # Knots location on each MRA level
    indiceshh  <- raster::extract(globraster[[i+1]],hh,cellnumbers=TRUE)[,1]
    hh_sf      <- hh_sf %>% mutate(!!paste0('R',i+1):=indiceshh)
    for(j in 0:nlevelsMRA){
      indiceshh <- raster::extract(globraster[[i+1]],
                  Qlist_sp[[j+1]],cellnumbers=TRUE)[,1]
      Qlist_sf[[j+1]] <- Qlist_sf[[j+1]] %>% 
                  mutate(!!paste0('R',i+1):=indiceshh)
    }
  }
  
  Qlist_sf[[nlevelsMRA+2]] <- hh_sf 
  
  for(i in 1:(nlevelsMRA+1)) {colnames(Qlist_sf[[i]])[1]  <- 
    'geometry';st_geometry(Qlist_sf[[i]]) <- 'geometry'}
  
#  comparacion_grupos <- combinations(n = nlevelsMRA+2, r = 2, 
#  repeats.allowed = F, v = 1:(nlevelsMRA+2))
#  comparacion_grupos <- rbind(cbind(1:(1+nlevelsMRA),
#  1:(1+nlevelsMRA)),comparacion_grupos)

Tm_matrices <- list() 
#MRA-block modulating function. Formula (2.7) in K&G-2020.
TSm_matrices <- list() 
#MRA-block modulating function in the case of S as first argument.

for(m in 0:nlevelsMRA){
  Tm_matrices[[m+1]] <- list()
  TSm_matrices[[m+1]] <- list()
  for(l in 0:m){
    Tm_matrices[[m+1]][[l+1]] <- list()
    if(l==0){
      k <- -1
      Tm_matrices[[m+1]][[l+1]][[k+2]] <- 1
      TSm_matrices[[m+1]][[k+2]] <- 1
    }else{
      for(k in 0:(l-1)){
        Tm_matrices[[m+1]][[l+1]][[k+2]] <- 0+outer(unlist(
          st_drop_geometry(Qlist_sf[[m+1]] %>% 
          dplyr::select(!!paste0('R',k+1)))),
          unlist(st_drop_geometry(Qlist_sf[[l+1]] %>% 
          dplyr::select(!!paste0('R',k+1)))),"==")
        
        Tm_matrices[[m+1]][[l+1]][[k+2]] <- Matrix(kronecker(
          Tm_matrices[[m+1]][[l+1]][[k+2]],
          matrix(rep(1,nCov_v*nCov_v),nrow = nCov_v)))
        
        TSm_matrices[[m+1]][[k+2]] <- 0+outer(unlist(
          st_drop_geometry(Qlist_sf[[nlevelsMRA+2]] %>% 
          dplyr::select(!!paste0('R',k+1)))),
          unlist(st_drop_geometry(Qlist_sf[[m+1]] %>% 
          dplyr::select(!!paste0('R',k+1)))),"==") 
        
        TSm_matrices[[m+1]][[k+2]] <- Matrix(kronecker(
          TSm_matrices[[m+1]][[k+2]],
          matrix(rep(1,nCov_v*nCov_v),nrow = nCov_v)))
      }
    }
  }
}
  return(list(Tm_matrices,TSm_matrices,Qlist_sf,nlevelsMRA))
}
  
