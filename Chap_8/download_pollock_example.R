
devtools::install_github("afsc-gap-products/coldpool")
devtools::install_github("James-Thorson/FishData")

# Download pollock data
library(FishData)
EBS = download_catch_rates( survey="Eastern_Bering_Sea", species_set = "Gadus chalcogrammus",
                      localdir = "C:/Users/James.Thorson/Desktop/Work files/Collaborations/2022 -- Spatial textbook" )
NBS = download_catch_rates( survey="Northern_Bering_Sea", species_set = "Gadus chalcogrammus",
                      localdir = "C:/Users/James.Thorson/Desktop/Work files/Collaborations/2022 -- Spatial textbook" )

# Combine and save
pollock = rbind( cbind("survey"="east",EBS), cbind("survey"="north",NBS) )
saveRDS( pollock, file="C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_8/pollock.rds" )

# Download cold pool extent
coldpool = coldpool:::cold_pool_index[,1:2]
saveRDS( coldpool, file="C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_8/coldpool.rds" )

#
EBS = st_read( "C:/Users/James.Thorson/Desktop/Git/FishStatsUtils/inst/region_shapefiles/EBSshelf/EBSshelf.shp" )
NBS = st_read( "C:/Users/James.Thorson/Desktop/Git/FishStatsUtils/inst/region_shapefiles/NBS/NBS.shp" )
survey_domain = st_union( EBS, NBS )
survey_domain = st_transform( survey_domain, crs="+proj=longlat +datum=WGS84" )
# remove holes from polygon for survey area
survey_domain = st_multipolygon(lapply(survey_domain[[7]], function(x) x[1]))
#st_write( BS, "C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_8/BS.shp" )
saveRDS( survey_domain, "C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_8/survey_domain.rds" )

