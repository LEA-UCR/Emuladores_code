library(sf)
library(raster) 
library(purrr)
source("0.generate_data.R") # generating betas takes a while

## Data generated:

### y is the response variable with regional resolution
### hh4 are     the xs with global resolution

dev.off()
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
dev.off()
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

indicesglob  <- list()
indicesreg   <- list()
cellloc.glob <- list()
cellloc.reg  <- list()

# ¿Cuántas veces queremos partir el dominio y cuál es el borde?
nn <- length(nc)
i<-5
temp <- raster(xmn=bordes[1],ymn=bordes[2],xmx=bordes[3],
               ymx=bordes[4],val=partitions[[i]],
               crs=crsglobal,ncols=nc[[i]],nrows=nr[[i]])
#plot(temp)
plot(rasterToPolygons(temp),add=T)
#points(globalpoints, cex=0.05,col="blue")

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

regionalpoints <- SpatialPointsDataFrame(
  cbind(all$coo$lonreg,all$coo$latreg), 
  as.data.frame(cbind(longlo= all$coo$longlo,
                      latglo= all$coo$latglo,
                      IDglo = all$coo$IDglo,
                      Yresp = all$resp,
                      Xcov  = all$covariate,
                      iP1   = all$iP1,
                      iP2   = all$iP2,
                      iP3   = all$iP3,
                      iP4   = all$iP4)),
  proj4string = CRS('+proj=longlat +datum=WGS84'),
  bbox = bordes)

regionalpoints = st_as_sf(regionalpoints)
plot(regionalpoints)
plot(regionalpoints %>% filter(iP4==1) )

globalpoints <- SpatialPointsDataFrame(
  cbind(all$coo$longlo,all$coo$latglo), 
  as.data.frame(cbind(lonreg= all$coo$lonreg,
                      latreg= all$coo$latreg,
                      IDglo = all$coo$IDglo,
                      Yresp = all$resp,
                      Xcov  = all$covariate,
                      iP1   = all$iP1,
                      iP2   = all$iP2,
                      iP3   = all$iP3,
                      iP4   = all$iP4)),
  proj4string = CRS('+proj=longlat +datum=WGS84'),
  bbox = bordes)

globalpoints = st_as_sf(globalpoints)
plot(globalpoints)
plot(globalpoints %>% filter(iP4==4) )

## How many nodes per partition?

#Knots for each cell:
generate_samples <- function(data,knots) 
  suppressMessages(st_sample(data, size = knots))
create_knots <- function(partition, knots){
  if(partition==1){
points <- map(globalpoints %>%split(.$iP1), 
              generate_samples,knots)
points <- imap(points, 
               ~st_sf(data_frame(iP = 
              rep(.y, length(.x))),geometry = .x))}
  if(partition==2){
    points <- map(globalpoints %>%split(.$iP2), 
                  generate_samples,knots)
    points <- imap(points, 
                   ~st_sf(data_frame(iP = 
                   rep(.y, length(.x))),geometry = .x))}
  if(partition==3){
    points <- map(globalpoints %>%split(.$iP3), 
                  generate_samples,knots)
    points <- imap(points, 
                   ~st_sf(data_frame(iP = 
                  rep(.y, length(.x))),geometry = .x))}
  if(partition==4){
    points <- map(globalpoints %>%split(.$iP3), 
                  generate_samples,knots)
    points <- imap(points, 
                   ~st_sf(data_frame(iP = 
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

points <- knots3
pp<-ggplot() + 
  geom_sf(data = points%>% filter(iP==4|iP==64), 
          aes(colour = iP,
              fill = iP),
          size = 1) + 
  scale_color_brewer(type = "div", palette = 4) + 
  scale_fill_brewer(type = "div", palette = 4)

print(pp)

# How to extract the knots from the main table:
# Knots (1,2,3,4) are geometries, and as such can be compared
# to the main table's geometries:
which((globalpoints$geometry %in% st_cast(knots4[knots4$iP==3,]$geometry, "POINT")))

# Now I need to now what's the final table format, to make it
# work with this nodes.


