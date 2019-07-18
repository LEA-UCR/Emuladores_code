## Archivos necesarios:

library(fields)
library(spatstat)
library(maps)
library(mapdata)

load("domainfinal.Rdata")
load("coordinates.Rdata")
load("altitude.Rdata")


plot(final$lonval-360,final$latval, cex=0.07,
     ylim=c(10,80),xlim=c(-170,-20),
     xlab="Easting", ylab="Northing",
     main="")

map("worldHires",xlim=c(-170,-20), ylim=c(10,80), 
col=1, add=TRUE)
points(t(aa)[,2]-360,t(aa)[,1],col="blue", cex=0.05)

legend("bottomleft", title="Coordinates for Model",
       c("CCSM - Global","CRCM - Regional"),
       fill=c("black","blue"))
rect(-135,25,-100,60,border=5,cex=1.2)

quilt.plot(altitude[,2]-360,altitude[,1],altitude[,3],
           ylim=c(10,80),xlim=c(-170,-20),lwd=0.1,
           xlab="Easting", ylab="Northing",
           main="")
map("worldHires",xlim=c(-170,-20), ylim=c(10,80), 
    col=1, add=TRUE)
rect(-135,25,-100,60,border="red",lwd=3.5)
