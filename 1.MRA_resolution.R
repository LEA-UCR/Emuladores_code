library(sf)
library(raster) 
library(purrr)
library(tidyr)
source("0.generate_data.R") # generating betas takes a while

## Data generated:

### y is the response variable with regional resolution
### hh4 are     the xs with global resolution

#dev.off()
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
#dev.off()
### In order to use MRA to analyze this dataset, 
### we need to create the partitions and knots:

# create a raster with the points and a data table:
bordes <- bbox(hh2)
crsglobal <- CRS('+proj=longlat +datum=WGS84')

### Aquí pueden cambiarse las particiones dependiendo del problema.
partitions <- list(1,
                   c(1:4),
                   seq(1,64,1), 
                   seq(1,256,1))
length(partitions) # how many levels
nc <- list(1, 2, 8, 16)
nr <- list(1, 2, 8, 16)
nc_n <- c(1,2,2,2) 
nr_n <- c(1,2,2,2) 



indicesglob  <- list()
indicesreg   <- list()
cellloc.glob <- list()
cellloc.reg  <- list()
indicesglobK  <- list()
indicesregK   <- list()

# ¿Cuántas veces queremos partir el dominio y cuál es el borde?
nn <- length(nc)

for(i in 1:nn){
  globraster <- raster(xmn=bordes[1],ymn=bordes[2],xmx=bordes[3],
                       ymx=bordes[4],val=partitions[[i]],
                       crs=crsglobal,ncols=nc[[i]],nrows=nr[[i]])
  
  indicesglob[[i]]    <- raster::extract(globraster,globalpoints, cellnumbers=TRUE)[,1]
  cellloc.glob[[i]]   <- rowColFromCell(globraster,indicesglob[[i]])
  
  
  indicesreg[[i]]     <- raster::extract(globraster,regionalpoints, cellnumbers=TRUE)[,1]
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
}

# table for regional indices:
indicesregK   <- data.frame(matrix(unlist(indicesregK), ncol=nn,
                                   nrow=length(indicesregK[[1]])))
names(indicesregK) <- as.character(unlist(lapply(1:nn,
                                                 function(i)paste0("iK",i))))

all <- cbind(dataset,indicesregK)
str(all)


regionalpoints <- SpatialPointsDataFrame(
  cbind(all$coo$lonreg,all$coo$latreg), 
  as.data.frame(cbind(longlo= all$coo$longlo,
                      latglo= all$coo$latglo,
                      IDglo = all$coo$IDglo,
                      Yresp = all$resp,
                      Xcov  = all$covariate,
                      iK1   = all$iK1,
                      iK2   = all$iK2,
                      iK3   = all$iK3,
                      iK4   = all$iK4)),
  proj4string = CRS('+proj=longlat +datum=WGS84'),
  bbox = bordes)

regionalpoints = st_as_sf(regionalpoints)

# table for global indices:
datasetglob <- dataset$coo %>% dplyr::select(IDglo,latglo,longlo) %>% 
  distinct(IDglo,latglo,longlo)

indicesglobK   <- data.frame(matrix(unlist(indicesglobK), ncol=nn,
                                    nrow=length(indicesglobK[[1]])))
names(indicesglobK) <- as.character(unlist(lapply(1:nn,
                                                  function(i)paste0("iK",i))))

allglob <- cbind(datasetglob,indicesglobK)

globalpoints <- SpatialPointsDataFrame(cbind(allglob$longlo, 
                                             allglob$latglo), 
                                       as.data.frame(cbind(
                                         longlo= allglob$longlo,
                                         latglo= allglob$latglo,
                                         iK1   = allglob$iK1,
                                         iK2   = allglob$iK2,
                                         iK3   = allglob$iK3,
                                         iK4   = allglob$iK4)),
                                       proj4string = CRS('+proj=longlat +datum=WGS84'),
                                       bbox = bordes)

globalpoints = st_as_sf(globalpoints)

indicesW <- list()
indicesW[[1]] <- globalpoints %>% expand(iK1) %>% arrange(iK1)%>%
  mutate(iP1 = 1:n())
indicesW[[2]] <- globalpoints %>% expand(iK1,iK2) %>% 
  arrange(iK1,iK2) %>% mutate(iP2 = 1:n())
indicesW[[3]] <- globalpoints %>% expand(iK1,iK2,iK3) %>%
  arrange(iK1,iK2,iK3) %>% mutate(iP3 = 1:n())
indicesW[[4]] <- globalpoints %>% expand(iK1,iK2,iK3,iK4) %>%
  arrange(iK1,iK2,iK3,iK4) %>% mutate(iP4 = 1:n())

tablaindicesW <- indicesW[[4]]%>%
  left_join(indicesW[[3]]) %>% left_join(indicesW[[2]]) %>%
  left_join(indicesW[[1]])


globalpoints <- globalpoints %>% left_join(tablaindicesW) 
regionalpoints <- regionalpoints %>% left_join(tablaindicesW) 

## How many nodes per partition?

#Knots for each cell:
generate_samples <- function(data,knots) 
  suppressMessages(st_sample(data, size = knots))
create_knots <- function(partition, knots){
  if(partition==1){
points <- purrr::map(globalpoints %>%split(.$iP1), 
              generate_samples,knots)
points <- imap(points, 
               ~st_sf(tibble(iP = 
              rep(.y, length(.x))),geometry = .x))}
  if(partition==2){
    points <- purrr::map(globalpoints %>%split(.$iP2), 
                  generate_samples,knots)
    points <- imap(points, 
                   ~st_sf(tibble(iP = 
                   rep(.y, length(.x))),geometry = .x))}
  if(partition==3){
    points <- purrr::map(globalpoints %>%split(.$iP3), 
                  generate_samples,knots)
    points <- imap(points, 
                   ~st_sf(tibble(iP = 
                  rep(.y, length(.x))),geometry = .x))}
  if(partition==4){
    points <- purrr::map(regionalpoints %>%split(.$iP4), 
                  generate_samples,knots)
    points <- imap(points, 
                   ~st_sf(tibble(iP = 
                   rep(.y, length(.x))),geometry = .x))}
points <- do.call(rbind, points)
points <- points %>% group_by(iP) %>% summarise()
points %>% mutate(n_points = map_int(geometry, nrow))
return(points)
}

#iP1 -> sample
knots1<-create_knots(1,100)
#iP2 -> sample
knots2<-create_knots(2,50)
#iP3 -> sample
knots3<-create_knots(3,25)
#iP4 all
knots4<-create_knots(4,5)

#points <- knots3
#pp<-ggplot() + 
#  geom_sf(data = points%>% filter(iP==4|iP==64), 
#          aes(colour = iP,
#              fill = iP),
#          size = 1) + 
#  scale_color_brewer(type = "div", palette = 4) + 
#  scale_fill_brewer(type = "div", palette = 4)

#print(pp)

# How to extract the knots from the main table:
# Knots (1,2,3,4) are geometries, and as such can be compared
# to the main table's geometries:
#which((globalpoints$geometry %in% st_cast(knots4[knots4$iP==3,]$geometry, "POINT")))

# Now I need to now what's the final table format, to make it
# work with this nodes.


 

knots1_tb <- as.data.frame(st_coordinates(knots1))
knots2_tb <- as.data.frame(st_coordinates(knots2))
knots3_tb <- as.data.frame(st_coordinates(knots3))
knots4_tb <- as.data.frame(st_coordinates(knots4))

colnames(knots1_tb) <- colnames(knots2_tb) <- colnames(knots3_tb) <- colnames(knots4_tb) <- c('longlo','latglo','ID')

knots1_tb <- knots1_tb %>% left_join(regionalpoints) %>% distinct(longlo,latglo,.keep_all = T)
knots2_tb <- knots2_tb %>% left_join(regionalpoints)%>% distinct(longlo,latglo,.keep_all = T)
knots3_tb <- knots3_tb %>% left_join(regionalpoints)%>% distinct(longlo,latglo,.keep_all = T)
knots4_tb <- knots4_tb %>% left_join(regionalpoints)%>% distinct(longlo,latglo,.keep_all = T)

knotsMRA <- list()

knotsMRA[[1]] <- st_sf(knots1_tb)
knotsMRA[[2]] <- st_sf(knots2_tb)
knotsMRA[[3]] <- st_sf(knots3_tb)
knotsMRA[[4]] <- st_sf(knots4_tb)

##Covariate matrix

