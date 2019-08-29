#' Show Narccap datasets available to download
#'
#' This function creates a table with the information of the diferent datasets
#' available to download from R in NARCCAP web site.
#' The columns shows
#' + Nombre de la variable
#' + Abreviatura de la variable
#' + Año de Inicio de los datos contenidos
#' + Año final de los datos contenidos
#' + Tamaño del archivo
#' + Current or future

#' @keywords datasets, download, data
#' @export
#' @import rlang dplyr rvest tidyr curl xml2
#' @examples
#'
#' ###Not run
#' ShowData()

ShowData <- function(){
  urls <- "https://www.earthsystemgrid.org/dataset/narccap.crcm.output.ccsm.html"
  wp <- read_html(urls)

  links1 <- html_nodes(wp, ".paddingTop div a")
  links1 <- html_text(links1)
  links1 <- tibble(links1) %>% separate(links1, into=as.character(c(1:6)))%>%
    unite("ext3", c("5","6"), sep="", remove=F) %>% unite("ext2", c("3","4"), sep="-",remove=F) %>%
    unite("ext1", c("1","2"), sep=".", remove=F) %>%
    unite("ext0", c("ext2","ext3"), sep=".",remove=F) %>%
    unite("ext", c("ext1","ext0"), sep=".",remove=F) %>%
    select(c("ext","1","2", "ext2","4" ,"ext3")) %>%
    mutate_at(.vars=c("ext", "ext3"),tolower) %>% mutate_at(.vars=c("2"), toupper) ##Transformaciones para unificar formato
  cuadro <- NULL
  for ( i in 1: length(links1$ext)){
    url <- paste0("https://www.earthsystemgrid.org/dataset/",links1$ext[i],"/file.html")
    webpage <- read_html(url)

    links_html <- html_nodes(webpage,'td:nth-child(2) a') ##LEE DE HTML
    links <- html_text(links_html) ##PASO DE HTML A TEXTO

    mb_html <- html_nodes(webpage,"td:nth-child(3)")
    mb <- html_text(mb_html)
    Size <- mb[-1]

    var_html <- html_nodes(webpage,'.checkbox label')  ##LEE DEL HTML LAS VARIABLES
    var <- html_text(var_html)  ##PASA A TEXTO
    var <- gsub("\n","",var) ##TRANSFORMACIONES
    var <- str_trim(var)
    var <- tibble(var) %>% separate(var, into=c("Name", "Full Name"), sep=6)
    var$Name <- gsub(" ", "", var$Name)
    var$`Full Name` <- str_squish(var$`Full Name`)

    Table <- rep(links1$ext3[i], length(links))
    Time <- rep(links1$`4`[i], length(links))
    M1<- rep(links1$`1`[i], length(links))
    M2 <- rep(links1$`2`[i], length(links))
    M3 <- rep(links1$ext2[i], length(links))

    df <- tibble(links) %>% separate(links, as.character(c(1:6)), remove = F)
    options(warn=-1)
    if(is.na(df$`6`[1])==T){
      df <- df %>% select(-c("2","3","5","6")) %>% separate("4", c("Year","M"),4) %>% select(-"M") %>%
        mutate(`4`="NA")
    } else{
      df <- df %>% select(-c("2","3","6")) %>% separate("5", c("Year","M"),4) %>% select(-"M")
    }
    df$Year <- as.numeric(df$Year)
    df$YearEnd <- apply(as.matrix(df$Year), 1, function(x) ifelse(x==1968, 1970, x+5))

    df <- mutate(df,FullVarName=factor(df$`1`))
    levels(df$FullVarName) <- factor(var$`Full Name`)
    df <- cbind(df,Size, Table, Time, M1,M2,M3) ##Crea el cuadro de esta tabla
    cuadro <- bind_rows(cuadro, df) ##Pega los datos de la tabla anterior con la recien procesada
  }

  colnames(cuadro) <- c("link", "Ab.Name", "Year", "P#","YearEnd","Full Name","Size","Table", "Time","M1", "M2", "M3" ) ##Da nombre a las columnas
  return(cuadro)
}


