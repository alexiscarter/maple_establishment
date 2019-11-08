## Plot sampling points with elevation
## Made with help from https://github.com/KevCaz/mapsWithR

## Load libraries and data points
library(sp)
library(raster)
library(GISTools)  
library(maps) 

tmp <- read.csv('map/GPS_megantic.csv', sep = ';')
meg <- SpatialPointsDataFrame(
  tmp[c('lon', 'lat')],
  data = tmp[c('id', 'elev')],
  proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
)

## Get coordinates
coord <- coordinates(meg)

## Get altitude
altCAN <- raster::getData(name="SRTM", lon = coord[1,1], lat = coord[1,2], path="") # more precise

## Canadian boundaries
bouCAN2 <- raster::getData(country='CAN', level=2, path="")

## Québec only
bouCAN2_Q <- bouCAN2[bouCAN2@data$NAME_1 == "Québec",]

## Get info for sampling points
al_ov <- over(meg, bouCAN2_Q)

## Get district name of Mount saint joseph
id <- bouCAN2_Q@data$NAME_2 %in% al_ov$NAME_2

## Elevation of the district
ra_elv <- rasterize(bouCAN2_Q[id,], crop(altCAN, bouCAN2_Q[id,]@bbox), mask=TRUE)

## Longitude and latitude extent
box <- matrix(c(-71.110, 45.440, -71.110, 45.469), nrow=2, ncol=2)

## Make the plot plot
plot(t(box), type = 'n', asp = 1, xlab = 'Longitude (°)', ylab = 'Latitude (°)') 
contour(ra_elv, add=T, nlevels = 10, labcex = 1)
plot(meg, add=T, pch = c(15, 19, 17)[tmp$forest], cex= 1, col = 'black')
text(meg, labels = tmp$block, cex= 0.7, pos = 3)
north.arrow(xb=-71.1045, yb=45.445, len=0.0005, lab="N") 
map.scale(x=-71.107, y=45.442, ratio=FALSE, relwidth=0.1) 
# saved in pdf protrait 5*7in
