#' Visualize Narccap data in a Map
#'
#' narccapMAP function generates a map of the chosen variable in a specific month and year

#' @keywords Map, Visualization, Climate data
#' @export
#'
#' @param data Data frame objet generated with NC2DFR or NC2DFG functions
#' @param var Name of the variable chosen to be visualized
#' @param year Numeric value of the year chosen to be visualized
#' @param month Number of the month you chosen to be visualized
#' @import dplyr stringr  ncdf4 lubridate reshape2 sp ggmap ggplot2
#' @examples
#'
#' ###Not run
#' narccapMAP(`RegionalData`, "sic", 1985, 1)
#'

narccapMAP <- function(data, var, year, month) {
  lat <- c(min(data$lat), max(data$lat))
  long <- c(min(data$lon), max(data$lon)) - 360
  bbox <- make_bbox(long, lat, f = 0.05)
  b <-
    get_map(bbox,
            maptype = "toner-lite",
            source = "stamen",
            color = "bw")
  ggmap(b)
  YY <- year
  MM <- month
  mes <- c(
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec"
  )
  G6 <-
    ggmap(b) + geom_point(data = data[data$Year == YY &
                                         data$Month == MM, ],
                          aes(lon - 360, lat, color = eval(as.name(var))),
                          alpha = 0.9) +
    scale_color_viridis_c(var) +
    labs(
      x = "Longitud",
      y = "Latitud",
      title = paste(var, "at", mes[MM], "of", YY)
    ) +
    theme_classic()
  return(G6)
}
