library(tidyverse)
library(sf)
library(sp)
library(raster)
library(RandomFields)



load('../datos/NARCCAP/PRtotal.RData')

base_cruda <- PRtot
rm(PRtot)
Year_p <- 1968:1999
month_p <- 1:12
fechas_p <- expand.grid(Year_p,month_p)
colnames(fechas_p) <- c('Year_p','month_p')

fechas_p <- fechas_p %>% arrange(Year_p,month_p)

compute.D <- function(indice){
  base_filtrada <- base_cruda %>% filter(Year==fechas_p[indice,1],Month==fechas_p[indice,2]) %>%
    mutate(Y=PRregional,X=PRglobal) %>% dplyr::select(-PRregional,-PRglobal)
  
  datos_filt <- base_filtrada %>% dplyr::select(Y,X)
  coordenadas <- base_filtrada %>% dplyr::select(lon,lat)
  
  hh <- SpatialPointsDataFrame(coords = coordenadas,data = datos_filt)
  proj4string(hh) <- '+proj=longlat +datum=WGS84'
  bordes <- bbox(hh)
  
  
  # Partition over the simulation area ----
  
  npartitions <- 50
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
    
    if(dim(subset_hh)[1]==0){
      return(c(NA,NA))
    }else{
      modelo_subset <- lm(Y~scale(X),data = subset_hh)
      return(coefficients(modelo_subset))
    }
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
  centroides_fractalB1 <- na.omit(centroides_fractalB1)
  coord <- st_coordinates(centroides_fractalB1)
  fractalD <- RFfractaldim(data = as_Spatial(centroides_fractalB1),bin=seq(0,60,1),mode = 'nographics')
  #RFhurst(x = coord[,1],y=coord[,2],data = centroides_fractalB1$beta1,method='dfa')
  rr <- RFvariogram(data = as_Spatial(centroides_fractalB1))
  D <- fractalD$vario$D
  H <- 3-D
  return(c(D,H))
}


resultados <- NULL
for(j in 1:dim(fechas_p)[1]){
  show(j)
  resultadost <- compute.D(j)
  resultados <- rbind(resultados,resultadost)
}

colnames(resultados)<- c('D','H')   
resultados <- cbind(fechas_p[-384,],resultados)
rownames(resultados) <- NULL

library(lubridate)
resultados <-
  resultados %>% mutate(fechas = ISOdate(year = Year_p, month = month_p, day = 1))

plotD <- ggplot(data = resultados,mapping = aes(x = fechas,y = D))+
  geom_line(lwd=0.75)+theme_bw()+xlab('Year')+ylab('Fractal Index')+
  ggtitle('Precipitation',)
