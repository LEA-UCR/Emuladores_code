library(raster) 
library(purrr)
library(tidyr)
library(dplyr)
library(sf)
load("dataset_general.Rdata") # data loading


# parameters of the raster using the spatial structure
bordes <- bbox(hh)
crsglobal <- CRS('+proj=longlat +datum=WGS84')

# MRA definitions: levels and partitions (4 levels and 
# 4 partitions: default)
partitions <- list(1,
                   c(1:4),
                   seq(1,16,1), 
                   seq(1,64,1))
length(partitions) # how many levels
nc <- list(1, 2, 4, 8)
nr <- list(1, 2, 4, 8)
nc_n <- c(1,2,2,2) 
nr_n <- c(1,2,2,2) 
nn <- length(nc)


# tic()
# for(i in 1:nn){
#   globraster <- raster(xmn=bordes[1],ymn=bordes[2],xmx=bordes[3],
#                        ymx=bordes[4],val=partitions[[i]],
#                        crs=crsglobal,ncols=nc[[i]],nrows=nr[[i]])
#   
#   indicesglob[[i]]    <- raster::extract(globraster,globalpoints, 
#                                          cellnumbers=TRUE)[,1]
#   cellloc.glob[[i]]   <- rowColFromCell(globraster,indicesglob[[i]])
#   
#   
#   indicesreg[[i]]     <- raster::extract(globraster,regionalpoints, 
#                                          cellnumbers=TRUE)[,1]
#   cellloc.reg[[i]]    <- rowColFromCell(globraster,indicesreg[[i]])
#   
#   indicesglobtemp <- as.data.frame(cellloc.glob[[i]]) %>% 
#     mutate(rown = row%%nr_n[i],coln=col%%nc_n[i]) %>% 
#     mutate(rown=ifelse(rown==0,nr_n[i],rown),
#            coln=ifelse(coln==0,nc_n[i],coln))
#   
#   indicesregtemp <- as.data.frame(cellloc.reg[[i]]) %>% 
#     mutate(rown = row%%nr_n[i],coln=col%%nc_n[i]) %>% 
#     mutate(rown=ifelse(rown==0,nr_n[i],rown),
#            coln=ifelse(coln==0,nc_n[i],coln))
#   
#   indexmatrix <- as.data.frame(expand.grid(1:nr_n[i],1:nc_n[i]))
#   indexmatrix <- indexmatrix %>% dplyr::select(rown=Var1,coln=Var2)%>%
#     mutate(celln=1:(nr_n[i]*nc_n[i]))
#   
#   indicesglobK[[i]] <- as.numeric((indicesglobtemp %>% 
#                             left_join(indexmatrix,by = c('rown','coln')) %>%
#                             dplyr::select(celln))$celln)
#   indicesregK[[i]] <- as.numeric((indicesregtemp %>% 
#                             left_join(indexmatrix,by = c('rown','coln')) %>%
#                             dplyr::select(celln))$celln)
# }
# toc()  ##22s


# Indexing of the spatial structure of data over the MRA partitions
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
  return(list(indicesregK))
}

#tic()
INDICES <- pmap(list(nc, nr, partitions,nc_n,nr_n ), f,bordes = bordes, crsglobal = crsglobal)
indicesregK <- list(INDICES[[1]][[1]],INDICES[[2]][[1]],INDICES[[3]][[1]], INDICES[[4]][[1]])
#toc()
rm(INDICES)


# Indexing by cell (iK...)
indicesregK   <- data.frame(matrix(unlist(indicesregK), ncol=nn,
                                   nrow=length(indicesregK[[1]])))
names(indicesregK) <- as.character(unlist(lapply(1:nn,
                                                 function(i)paste0("iK",i))))

hh <- st_as_sf(hh)

# Indexing by partition (iP...)
indicesW <- list()

indicesW[[1]] <- indicesregK %>% 
  expand(iK1) %>% 
  arrange(iK1)%>%
  mutate(iP1 = 1:n())
indicesW[[2]] <- indicesregK %>% 
  expand(iK1,iK2) %>% 
  arrange(iK1,iK2) %>% 
  mutate(iP2 = 1:n())
indicesW[[3]] <- indicesregK %>% 
  expand(iK1,iK2,iK3) %>%
  arrange(iK1,iK2,iK3) %>% 
  mutate(iP3 = 1:n())
indicesW[[4]] <- indicesregK %>% 
  expand(iK1,iK2,iK3,iK4) %>%
  arrange(iK1,iK2,iK3,iK4) %>% 
  mutate(iP4 = 1:n())

# Everything together: iP and iK
tablaindicesW <- indicesW[[4]]%>%
  left_join(indicesW[[3]]) %>% left_join(indicesW[[2]]) %>%
  left_join(indicesW[[1]])

#Update of data indexing with iP
indicesregK <- indicesregK %>% left_join(tablaindicesW) 

#Update of spatial structure with data indexing
hh <- hh %>% bind_cols(indicesregK)

#Random generation of MRA knots 
generate_samples <- function(data,knots) 
  suppressMessages(st_sample(st_as_sfc(st_bbox(data)), size = knots))

#Random generation of MRA knots and partition identification 
create_knots <- function(partition, knots){
  if(partition==1){
    points <- purrr::map(hh %>%split(.$iP1), 
                         generate_samples,knots)
    points <- imap(points, 
                   ~st_sf(tibble(iP = 
                                   rep(.y, length(.x))),geometry = .x))}
  if(partition==2){
    points <- purrr::map(hh %>%split(.$iP2), 
                         generate_samples,knots)
    points <- imap(points, 
                   ~st_sf(tibble(iP = 
                                   rep(.y, length(.x))),geometry = .x))}
  if(partition==3){
    points <- purrr::map(hh %>%split(.$iP3), 
                         generate_samples,knots)
    points <- imap(points, 
                   ~st_sf(tibble(iP = 
                                   rep(.y, length(.x))),geometry = .x))}
  
  points <- do.call(rbind, points)
  points <- points %>% group_by(iP) %>% summarise()
  points %>% mutate(n_points = map_int(geometry, nrow))
  return(points)
}

# Number of knots according to Katzfuss et al, 2017.
nknots <- floor(dim(dataset)[1]/(4^4))+1
show(paste0('El numero de nodos X particion es: ',nknots))

# Random generation of knots per MRA level
knots1<-create_knots(1,nknots)
knots2<-create_knots(2,nknots)
knots3<-create_knots(3,nknots)

knots1_tb <- as.data.frame(st_coordinates(knots1))
knots2_tb <- as.data.frame(st_coordinates(knots2))
knots3_tb <- as.data.frame(st_coordinates(knots3))

colnames(knots1_tb) <- c('lon','lat','iP1')
colnames(knots2_tb) <- c('lon','lat','iP2')
colnames(knots3_tb) <- c('lon','lat','iP3')

# iP hierarchy: needed to compute the loglikelihood recursively 
knots1_tb <- knots1_tb
tablapadres1 <- tablaindicesW %>% dplyr::select(starts_with('iP')) %>%
  dplyr::select(-iP4,-iP3)%>% distinct()
knots2_tb <- knots2_tb %>% left_join(tablapadres1)
tablapadres2 <- tablaindicesW %>% dplyr::select(starts_with('iP')) %>%
  dplyr::select(-iP4)%>% distinct()
knots3_tb <- knots3_tb %>% left_join(tablapadres2)
knots4_tb <- hh #IMPORTANT: last level uses all the available data as knots

# Knots structure: 
knotsMRA <- list()

knotsMRA[[1]] <- st_as_sf(knots1_tb,coords = c(1,2))
knotsMRA[[2]] <- st_as_sf(knots2_tb,coords = c(1,2))
knotsMRA[[3]] <- st_as_sf(knots3_tb,coords = c(1,2))
knotsMRA[[4]] <- knots4_tb

## remove all values that are duplicated
rm(list=ls()[! ls() %in% c("bordes","indicesW", "knotsMRA",
                           "nn",'hh')])

