library(raster)
tif_files <- list.files(path="~/phenology-urban/data/raw/NY/dem", pattern="\\.tif$", full.names=TRUE)

library(sf)
polygons <- st_read("~/phenology-urban/data/raw/NY/NY_building_footprint/geo_export_56deb8a8-0463-4f00-8f45-99e751e3b8e1.shp")
polygons$dem <- polygons$groundelev + polygons$heightroof
polygons <- na.omit(polygons)




for (file in tif_files) {
  raster_data <- raster(file)
  raster_extent <- extent(raster_data)
  raster_proj <- projection(raster_data)
  raster_resolution <- res(raster_data)
  
  raster_output <- raster(raster_extent, res = raster_resolution, crs = raster_proj)
  rasterized <- rasterize(polygons, raster_output, field = "dem", fun = "count")
  rasterized2 <- rasterize(polygons, raster_data, field = "dem", fun = "count", update=TRUE) 
  
  outpath <- sub("be_", "bld2_be_", file)
  writeRaster(rasterized2, outpath, format = "GTiff")
  
  
}


