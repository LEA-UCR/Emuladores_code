library(raster) 
library(purrr)
library(tidyr)
library(dplyr)
library(sf)
library(RandomFields)

dataset_file <- 'sim_data/datasetSVCExponential1.Rdata'
load(dataset_file)
bordes <- bbox(hh)


# Partition over the simulation area ----

npartitions <- 100
npartitions_r <- 10
npartitions_c <- npartitions/npartitions_r
partitions <- seq(1,npartitions)

crsglobal <- CRS('+proj=longlat +datum=WGS84')

globraster <- raster(
  xmn = bordes[1],
  ymn = bordes[2],
  xmx = bordes[3],
  ymx = bordes[4],
  val = partitions,
  crs = crsglobal,
  ncols = npartitions_c,
  nrows = npartitions_r
)

indicesreg     <- raster::extract(globraster,hh, 
                                  cellnumbers=TRUE)[,1]
cellloc.reg    <- rowColFromCell(globraster,indicesreg)
hh <- st_as_sf(hh)
hh <- hh %>% mutate(indices=indicesreg)


coeficientes_hh <- function(i){
  subset_hh <- hh %>% 
    filter(indices==i)
  
  modelo_subset <- lm(Y~scale(X),data = subset_hh)
  
  return(coefficients(modelo_subset))
}
indice_resultados <- 1:npartitions
names(indice_resultados) <- 'indices'
coeff <- purrr::map_dfr(indice_resultados,coeficientes_hh)
coeff <- data.frame(t(coeff))
rownames(coeff) <- NULL
colnames(coeff) <- c('Beta0','Beta1')


grid_sf <- st_as_sf(as(globraster,'SpatialPolygonsDataFrame'))
centroides <- st_centroid(grid_sf)


indicescent     <- raster::extract(globraster,as_Spatial(centroides), 
                                  cellnumbers=TRUE)[,1]
centroides <- centroides %>% mutate(indices=indicescent) %>%
  mutate(beta0=coeff$Beta0,beta1=coeff$Beta1)

centroides_fractalB1 <- centroides %>% dplyr::select(beta1)
coord <- st_coordinates(centroides_fractalB1)
fractalD <- RFfractaldim(data = as_Spatial(centroides_fractalB1),bin=seq(0,1,0.01))

#RFhurst(x = coord[,1],y=coord[,2],data = centroides_fractalB1$beta1,method='dfa')
#rr <- RFvariogram(data = as_Spatial(centroides_fractalB1))

