#' Convert the narccap nc file in a data frame
#'
#' This  function coverts the .nc downloaded  from NARCCAP global models in a data.frame object.
#' This proccess may require high RAM memory capacity
#'
#' @keywords datasets, download, data
#' @export
#' @param PATH Is the directory where is located the .nc files of a same variable
#' @param VARNAME String, is the name of the variable in the dataset
#' @examples
#' NC2DF("~/PrecipitationGlobalModelData")
#' @import dplyr stringr  ncdf4 lubridate reshape2 sp

NC2DFG <- function(PATH, VARNAME){
  blatitude <- c(19.12639, 74.40000)
  blongitude <- c(198.6576, 326.4000)
  dirbase <- PATH
  listfilesg <-
    list.files(path = dirbase)[str_detect(list.files(path = dirbase), 'nc$')]
  var <- NULL

  c <- sub("\\_.*", "", listfilesg[1])
  c <- VARNAME
  fecha <- '1870-01-01'
  for (i in 1:length(listfilesg)) {
    show(paste0('Construccion datos mensuales-', i))
    vglobal <- nc_open(paste0(dirbase, listfilesg[i]))
    var_pre <- ncvar_get(vglobal, c)
    lonvar <- ncvar_get(vglobal, 'lon')
    latvar <- ncvar_get(vglobal, 'lat')
    timevar <- ncvar_get(vglobal, 'time')

    fechabase <- ymd(fecha)
    timevar <- fechabase + as.period(ddays(timevar))

    dimnames(var_pre)[[1]] <- lonvar
    dimnames(var_pre)[[2]] <- latvar
    dimnames(var_pre)[[3]] <- as.character(timevar)

    var_pre <- melt(var_pre)

    colnames(var_pre) <- c('lon', 'lat', 'Time', c)

    var_pre <-
      var_pre %>% filter(lat >= blatitude[1],
                         lat <= blatitude[2],
                         lon >= blongitude[1],
                         lon <= blongitude[2]) %>%
      mutate(Time = ymd(as.character(Time))) %>% filter(Time >= ymd('1968-01-01')) %>%
      mutate(Year = year(Time), Month = month(Time)) %>% group_by(Year, Month, lon, lat) %>%
      summarise(m = mean(eval(as.name(c)))) %>% ungroup()

    var <- bind_rows(var, var_pre)
  }
colnames(var) <- c("Year", "Month", "lon", "lat", VARNAME )
assign(paste(c, "Global"),var,envir = .GlobalEnv)
rm(var)

}
