#' Download .nc datasets from Narccap
#'
#' This function downloads the specified datasets from NARCCAP
#' @keywords datasets, download, data
#' @export
#' @import rlang dplyr  rvest  tidyr curl
#' @param table Is the table produced with the function ShowData()
#' @param path The directory where you are going to download the data, we recomend to create a new directory for each variable.
#' @param AbName The abbreviate name of the variable you are going to download
#' @param FromYear The initial year data sets
#' @param ToYear The last year  of data sets
#'
#'
#' @examples
#' ###Not run
#' table <- ShowData()
#' DownloadDataF(table=table, path="~/RegionalModelPrecipitation", AbName="pr", FromYear=1968, ToYear=1990)

DownloadDataF <-
  function(table, path=" ",
           AbName,
           FromYear = 1968,
           ToYear = 2071) {
    cuadro1 <-
      table %>% filter(Ab.Name == AbName,
                       Year >= FromYear,
                       YearEnd <= ToYear,
                       Time == Time)
    for (i in 1:length(cuadro1$link)) {
      curl_fetch_disk(
        paste0(
          "https://tds.ucar.edu/thredds/fileServer/datazone/narccap/DATA/CRCM/",
          cuadro1$M3[i],
          "/",
          cuadro1$Table[i],
          "/",
          cuadro1$link[i]
        ),
        paste0("~/", cuadro1$link[i])
      )
    }
  }
