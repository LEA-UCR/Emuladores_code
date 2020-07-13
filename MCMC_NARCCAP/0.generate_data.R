library(tidyverse)
library(sf)
library(sp)

variable_narccap <- 'Prec' # Temp/Prec

if(variable_narccap=='Temp'){
  load('./data_narccap/TStotal.RData')
  base_cruda <- TS_tot
  rm(TS_tot)
}else{
  load('./data_narccap/PRtotal.RData')
  base_cruda <- pr_tot
  rm(pr_tot)
}


if(variable_narccap=='Temp'){
  datos_filt <- base_cruda %>% rename(Y=ts) %>% 
    select(Y,TREFHT,OMEGA,PSL,U,V) %>% mutate(Y=log(log(Y)))
  coordenadas <- base_cruda %>% select(lon,lat)
}else{
  datos_filt <- base_cruda %>% rename(Y=pr) %>% 
    select(Y,TREFHT,OMEGA,PSL,U,V) %>% mutate(Y=log(Y))
  coordenadas <- base_cruda %>% select(lon,lat)
}

dataset <- datos_filt

hh <- SpatialPointsDataFrame(coords = coordenadas,data = datos_filt)
proj4string(hh) <- '+proj=longlat +datum=WGS84'
bordes <- bbox(hh)

save(dataset, hh, file=paste0("data_narccap/dataset",
                              variable_narccap,".Rdata"))
