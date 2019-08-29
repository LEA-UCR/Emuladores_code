#' Join the regional and global datasets
#'
#' This function joins the processed data from regional and global model in a Data Frame that matches the points using grids
#'
#' @keywords joining, global model, regional model, data
#' @export
#' @param var is the dataframe produced by the function NC2DFG()
#' @param varreg is the dataframe produced by the function NC2DFG()
#' @examples
#'
#' ###Not run
#' joinRGdata(globaldata, regionaldata)

joinRGdata <- function(var,varreg){
  c <- colnames(varreg)[3]
  c1 <- colnames(var)[5]
  gridglobal <- var %>% select(lon,lat) %>% distinct(lat,lon) ##Genera la cuadricula de los datos del modelo global.
  gridglobalsp <- SpatialPixels(SpatialPoints(gridglobal),tolerance = 1e-4)

  pointsregional <- varreg %>% select(lon,lat) %>% distinct(lat,lon) ##Genera la cuadricula de los datos del modelo regional.
  pointsregionalsp <- SpatialPoints(pointsregional)

  indicesgrid <- over(pointsregionalsp,gridglobalsp) ##Se genera un indice para cada cuadrícula para agrupar los datos globales y regionales dentro de una misma cuadricula

  pointsregional <- pointsregional %>% mutate(indicegrid = indicesgrid)
  gridglobal <- gridglobal %>% mutate(indicegrid=seq(1,dim(gridglobal)[1]))

  var <- var %>% left_join(gridglobal,by = c("lon", "lat")) ##Le agrega el indice a ambos datos para posteriormente unirlos
  varreg <- varreg %>% left_join(pointsregional,by = c("lon", "lat"))
  colnames(var)[5] <- paste(c,'global')
  colnames(varreg)[3] <- paste(c1,'regional')

  VTot <- varreg %>% left_join(var,by = c("Year", "Month", "indicegrid")) %>%
    select(Year,Month,paste(c1,'regional'),lat=lat.x,lon=lon.x,paste(c,'global'), indicegrid) ##Une los archivos por medio del año, mes e indice de cuadricula

  assign("JoinnedDF",VTot, envir=.GlobalEnv)
  rm(VTot)

}
