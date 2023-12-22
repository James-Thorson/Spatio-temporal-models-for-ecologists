
library(raster)
library(sp)

LL <- locator(100)
# saveRDS(LL, file=paste0(root_dir,"LL.rds"))
region_extent <- data.frame(long=LL$x, lat=LL$y)
## Need to duplicate a point so that it is connected
region_extent <- rbind(region_extent, region_extent[1,])

#### Turn it into a spatial polygon object
## https://www.maths.lancs.ac.uk/~rowlings/Teaching/Sheffield2013/cheatsheet.html
poly <- Polygon(region_extent)
polys <- Polygons(list(poly), ID='all')
sps <- SpatialPolygons(list(polys))
#sps <- SpatialPolygonsDataFrame(sps, data.frame(Id=factor('all'), F_AREA=1, row.names='all'))
proj4string(sps)<- CRS("+proj=longlat +datum=WGS84")

raster::shapefile( sps, filename="red_snapper_extent", overwrite=TRUE )
plot(sps,add=TRUE, col=rgb(0,0,1,0.2))
