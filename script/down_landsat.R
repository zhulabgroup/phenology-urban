library(devtools)
library(raster)
library(rgdal)
# library(gdalUtils)
library(sf)
library(mapview)
library(sp)
#install packages
#install.packages("landsat")
#install.packages("LST")
devtools::install_github("frederikzi/getSpatialData") 

library(getSpatialData)
#explanation for getspatialData:
#https://www.rdocumentation.org/packages/getSpatialData/versions/0.1.0
library(landsat) #process landsat data
library(LST) #claculate LST from landsat bands 10 and 11
#check if server is available
services()
#login to USGS
login_USGS(username = "JialiUMICH", 
           password="19980126_Lily",
           n_retry = 10,
           verbose = T)
login_earthdata(username = "jializ98", 
                password = "19980126_Lily",
                n_retry = 10,
                verbose=T)

#set archive directory
set_archive("~/phenology-urban/data/raw/NY/LST_L8")

roi_sf <- sf::read_sf('data/raw/NY/boundary/nybb_dissolved.shp')
ny_boundary <- st_transform(roi_sf, crs = st_crs(4326))
ny_geometry <- st_geometry(ny_boundary)
ny_spatial <- as(ny_geometry, "Spatial")
ny_sp <- as(ny_spatial, "SpatialPolygons")

# Set AOI
set_aoi(aoi = ny_sp)

#download 
get_products("Landsat")
getLandsat_records(time_range = c("2018-01-01", "2018-12-31"),
                   products = "landsat_ot_c2_l2",
                   check_products = T)
getLandsat_data()
