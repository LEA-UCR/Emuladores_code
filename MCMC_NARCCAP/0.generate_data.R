library(tidyverse)
library(sf)
library(sp)

variable_narccap <- 'Temp' # Temp/Prec

if(variable_narccap=='Temp'){
  load('../datos/NARCCAP/TStotal.RData')
  base_cruda <- TStot
  rm(TStot)
}else{
  load('../datos/NARCCAP/PRtotal.RData')
  base_cruda <- PRtot
  rm(PRtot)
}


Year_p <- 1968:1999
month_p <- 1:12
fechas_p <- expand.grid(Year_p,month_p)
colnames(fechas_p) <- c('Year_p','month_p')

fechas_p <- fechas_p %>% arrange(Year_p,month_p)

indice <- dim(fechas_p)[1]
if(variable_narccap=='Temp'){
  base_filtrada <- base_cruda %>% filter(Year==fechas_p[indice,1],Month==fechas_p[indice,2]) %>%
    mutate(Y=TSregional,X=TSglobal) %>% dplyr::select(-TSregional,-TSglobal)
}else{
  base_filtrada <- base_cruda %>% filter(Year==fechas_p[indice,1],Month==fechas_p[indice,2]) %>%
    mutate(Y=PRregional,X=PRglobal) %>% dplyr::select(-PRregional,-PRglobal)
}

datos_filt <- base_filtrada %>% dplyr::select(Y,X)
coordenadas <- base_filtrada %>% dplyr::select(lon,lat)

dataset <- datos_filt
hh <- SpatialPointsDataFrame(coords = coordenadas,data = datos_filt)
proj4string(hh) <- '+proj=longlat +datum=WGS84'
bordes <- bbox(hh)

save(dataset, hh, file=paste0("data_narccap/dataset",
                              variable_narccap,".Rdata"))
