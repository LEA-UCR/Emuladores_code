library(stringr)
library(ncdf4)
library(lubridate)
library(reshape2)
library(dplyr)
library(tidyr)

# Acceso a la Tabla 4 regional (fijos): https://www.narccap.ucar.edu/data/table4/ orog es la altitud.
setwd("~/Dropbox/Emuladores")

altitude <- ncdf4::nc_open("orog_CRCM.nc")
print(altitude)
values <- ncdf4::ncvar_get(altitude, "orog")
lat <- ncdf4::ncvar_get(altitude, "lat")
lon <- ncdf4::ncvar_get(altitude, "lon")
