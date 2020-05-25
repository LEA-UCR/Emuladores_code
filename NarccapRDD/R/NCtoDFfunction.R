#' Convert the narccap nc file in a data frame
#'
#' This  function coverts the .nc downloaded  from NARCCAP regional models in a data.frame object.
#' This proccess may require high RAM memory capacity

#' @keywords datasets, download,  data
#' @export
#'
#' @param PATH Is the directory where is located the .nc files of a same variable
#' @import tidyr dplyr stringr lubridate ncdf4 sp
#' @examples
#' NC2DFR("~/RegionalModelPrecipitation")

NC2DFR <- function(PATH){
  dirbaseregional <- PATH
  listfilesregional <- list.files(path = dirbaseregional)[str_detect(list.files(path = dirbaseregional),'nc$')]
  varreg <- NULL
  blatitudet <- NULL
  blongitudet <- NULL

  c1 <- separate(tibble(listfilesregional),1, sep="_", as.character(c(1:4))) %>% select("1")
  c1 <- paste(c1[1,1])

  for(i in 1:length(listfilesregional)){
    show(paste0('Construccion datos mensuales-Regional-',i))
    regional <- ncdf4::nc_open(paste0(dirbaseregional,listfilesregional[i]))  ##Leemos el archivo .nc utilizando la lista previamente creada
    varreg_pre <- ncdf4::ncvar_get(regional,c1)
    lonvar <- ncdf4::ncvar_get(regional,'lon')
    latvar <- ncdf4::ncvar_get(regional,'lat')
    xcvar <- ncdf4::ncvar_get(regional,'xc')
    ycvar <- ncdf4::ncvar_get(regional,'yc')
    timevar <- ncdf4::ncvar_get(regional,'time')  ##Hasta acá lo que se hace es obtener las variables espaciales, temporales y la de interés del archivo de datos
    blatitude <- c(min(latvar)-1.4,floor(max(latvar))+1.4)
    blatitudet <- rbind(blatitudet,blatitude)
    blongitude <- c(min(lonvar)-1.4,floor(max(lonvar))+1.4)
    blongitudet <- rbind(blongitudet,blongitude)

    fechabase <- lubridate::ymd('1968-01-01') ##Se establece la fecha desde la que comienza a contar el tiempo en el modelo
    timevar <- fechabase+lubridate::as.period(ddays(timevar)) ##Transforma la variable de tiempo a formato Año mes día
    dimnames(varreg_pre)[[1]] <- xcvar
    dimnames(varreg_pre)[[2]] <- ycvar
    dimnames(varreg_pre)[[3]] <- as.character(timevar) ##Define nombres a los ejes del conjunto de los datos

    varreg_pre <- reshape2::melt(varreg_pre) ##Reorganiza los datos

    dimnames(latvar)[[1]] <- xcvar
    dimnames(lonvar)[[1]] <- xcvar
    dimnames(latvar)[[2]] <- ycvar
    dimnames(lonvar)[[2]] <- ycvar

    latlondata <- expand.grid(xcvar,ycvar)
    latlondatapre <- t(sapply(X=1:dim(latlondata)[1], FUN = function(x) return(c(latvar[as.character(latlondata[x,1]),as.character(latlondata[x,2])],lonvar[as.character(latlondata[x,1]),as.character(latlondata[x,2])]))))
    latlondata <- cbind(latlondata,latlondatapre)
    colnames(latlondata) <- c('xc','yc','lat','lon')
    colnames(varreg_pre) <- c('xc','yc','Time',c1) ##Se ponen Nombres a las columnas

    varreg_pre <- varreg_pre %>% mutate(Time=ymd_hms(as.character(Time))) %>%
      mutate(Year = year(Time),Month=month(Time)) %>%
      group_by(Year,Month,xc,yc) %>% summarise(m=mean(eval(as.name(c1)))) %>% ungroup()  ##Convierte la variable time en tres variables que indican el año y mes y saca la media de las observaciones para ese mes.

    varreg_pre <- varreg_pre %>% left_join(latlondata,by = c('xc', 'yc')) %>%
      select(-xc,-yc) ##Une al cuadro las variables espaciales

    varreg_pre <- varreg_pre %>% mutate(ID = rep(i, dim(varreg_pre)[1]))

    varreg <- bind_rows(varreg,varreg_pre) ##Uno los datos del archivo con los que fueron trabajados anteriormente en el ciclo.
  }


  a <- group_by(varreg,Year, Month, lat,lon) %>% summarize(n= n())
  repetidos <- subset(a, n!=1)

  a <- unique(repetidos$Year)
  ids <- as.numeric(unique(varreg$ID))[-(length(a)+1)]
  for (i in 1:length(ids)){
   varreg <- subset(varreg, !(Year==a[i]&Month==1&ID==ids[i]))
  }
  colnames(varreg) <- c("Year", "Month", c1, "lat", "lon", "ID")
  assign(x=paste(c1, "Regional"),value=varreg, envir = .GlobalEnv)
}



