
setwd( R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Chap_10)' )

library(rasterdiv)
library(terra)

NDVI = rast(copNDVI)
terra::writeRaster( NDVI, filename="NDVI.tif" )

