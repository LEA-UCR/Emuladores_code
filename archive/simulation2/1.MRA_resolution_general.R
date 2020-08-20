library(raster) 
library(purrr)
library(tidyr)
library(dplyr)
library(sf)

gen_resolution <- function(dataset_file){
load(datasetfile) # data loading
# parameters of the raster using the spatial structure
bordes <- bbox(hh)
crsglobal <- CRS('+proj=longlat +datum=WGS84')
# MRA definitions: levels and partitions (4 levels and 
# 4 partitions: default)
nlevelsMRA <- 2 #Number of MRA levels
npartitions <- 4
npartitions_r <- 2
npartitions_c <- npartitions/npartitions_r
partitions <- list()
nc <- nr <- NULL
nc_n <- nr_n <- NULL
for(i in 1:nlevelsMRA){
  partitions[[i]] <- seq(1,npartitions^(i-1))
  nc[i] <- npartitions_c^(i-1)
  nr[i] <- npartitions_r^(i-1)
}
nc_n <- c(1,rep(npartitions_c,nlevelsMRA-1))
nr_n <- c(1,rep(npartitions_r,nlevelsMRA-1))
nn <- length(nc)

f <- function(nc, nr, partitions, nc_n, nr_n, bordes, crsglobal){
  globraster <- raster(xmn=bordes[1],ymn=bordes[2],xmx=bordes[3],
                       ymx=bordes[4],val=partitions,
                       crs=crsglobal,ncols=nc,nrows=nr)
  indicesreg     <- raster::extract(globraster,hh, 
                                    cellnumbers=TRUE)[,1]
  cellloc.reg    <- rowColFromCell(globraster,indicesreg)
  indicesregtemp <- as.data.frame(cellloc.reg) %>% 
    mutate(rown = row%%nr_n,coln=col%%nc_n) %>% 
    mutate(rown=ifelse(rown==0,nr_n,rown),
           coln=ifelse(coln==0,nc_n,coln))
  indexmatrix <- as.data.frame(expand.grid(1:nr_n,1:nc_n))
  indexmatrix <- indexmatrix %>% dplyr::select(rown=Var1,coln=Var2)%>%
    mutate(celln=1:(nr_n*nc_n))
  indicesregK <- as.numeric((indicesregtemp %>% 
                 left_join(indexmatrix,by = c('rown','coln')) %>%
                 dplyr::select(celln))$celln)
  return(indicesregK)
}

indicesregK <- pmap(list(nc, nr, partitions,nc_n,nr_n ), f,
            bordes = bordes, crsglobal = crsglobal)

# Indexing by cell (iK...)
indicesregK   <- data.frame(matrix(unlist(indicesregK), ncol=nn,
                                   nrow=length(indicesregK[[1]])))
names(indicesregK) <- as.character(unlist(lapply(1:nn,
                                   function(i)paste0("iK",i))))

hh <- st_as_sf(hh)
# Indexing by partition (iP...)

indicesW.make <- function(indice){
  combs <- indicesregK[,1:indice] %>% crossing() %>%
    mutate(ttemp=1:n())
  colnames(combs) <- c(paste0('iK',1:indice),paste0('iP',indice))
  return(combs)
}

indicesW <- purrr::map(1:nlevelsMRA,indicesW.make)

# Everything together: iP and iK
tablaindicesW <- indicesW %>% reduce(left_join)

#Update of data indexing with iP
indicesregK <- indicesregK %>% left_join(tablaindicesW) %>%
  mutate(indice_m=1:n())
#Update of spatial structure with data indexing
hh <- hh %>% bind_cols(indicesregK)
#Random generation of MRA knots 
generate_samples <- function(data,knots){ 
  suppressMessages(st_sample(st_as_sfc(st_bbox(data)), size = knots))
}

library(rlang)
# Tomado de https://www.natedayta.com/2018/03/04/split-a-tidyverse-incarnation-of-split/
split_ <- function(data, ..., .drop = TRUE) {
  vars <- ensyms(...)
  vars <- purrr::map(vars, function(x) eval_tidy(x, data))
  split(data, vars, drop = .drop)
}

#Random generation of MRA knots and partition identification 
create_knots <- function(partition, knots){
  vv <-  c(quo('iP1'),quo('iP2'),quo('iP3'),quo('iP4'),quo('iP5'),
           quo('iP6'),quo('iP7'),quo('iP8'),quo('iP9'),quo('iP10'))
  ##Lo de arriba definitivamente debe mejorarse
  points <- purrr::map(hh %>%split_(!!vv[[partition]]), 
                         generate_samples,knots)
  points <- imap(points, 
                   ~st_sf(tibble(iP = 
                  rep(.y, length(.x))),geometry = .x))
  points <- do.call(rbind, points)
  points <- points %>% group_by(iP) %>% summarise()
  points %>% mutate(n_points = map_int(geometry, nrow))
  return(points)
}

# Number of knots according to Katzfuss et al, 2017.
# nknots <- floor(dim(dataset)[1]/(4^4))+1
nknots <- 10
show(paste0('El numero de nodos X particion es: ',nknots))
# Random generation of knots per MRA level

order.knots <- function(indice,nknots){
  knots<-create_knots(indice,nknots)
  knots_tb <- as.data.frame(st_coordinates(knots))
  colnames(knots_tb) <- c('lon','lat',paste0('iP',indice))
  return(knots_tb)
}

knots_tb <- purrr::map(1:(nlevelsMRA-1),order.knots,nknots=nknots)


# iP hierarchy: needed to compute the loglikelihood recursively 

knots.parent <- function(indice){
  if(indice==1){
    knots_tb[[indice]] <- knots_tb[[indice]] %>% mutate(indice_m=1:n())
    knots_tb[[indice]] <- st_as_sf(knots_tb[[indice]],coords = c(1,2))
  }else{
    if(indice==nlevelsMRA){
      knots_tb[[indice]] <- hh
    }else{
      tablapadres <- tablaindicesW %>% dplyr::select(starts_with('iP')) %>%
        dplyr::select(num_range('iP',1:indice))%>% distinct() 
      knots_tb[[indice]] <- knots_tb[[indice]] %>% left_join(tablapadres) %>% mutate(indice_m=1:n())
      knots_tb[[indice]] <- st_as_sf(knots_tb[[indice]],coords = c(1,2))
    }
  }
  return(knots_tb[[indice]])
}

knotsMRA <- purrr::map(1:nlevelsMRA,knots.parent)

## remove all values that are duplicated
rm(list=ls()[! ls() %in% c("bordes","indicesW", "knotsMRA",
                           "nn",'hh')])
return(list(bordes,indicesW,knotsMRA,nn,hh))
}


