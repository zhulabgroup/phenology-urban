library(sf)
library(raster)

# Read shapefile
boundary <- st_read("~/lab-data/phenology/phenology_urban/NY/boundary/nybb_dissolved.shp")

for (i in 2012:2022) {
  # Read raster file
  input_dir <- stringr::str_c("~/lab-data/phenology/phenology_urban/NY/nlcd_tcc/raw/nlcd_tcc_conus_",i,"_v2021-4.tif")
  raster_file <- raster(input_dir)
  # same crs
  boundary <- st_transform(boundary, crs(raster_file))
  # clip
  clipped_raster <- crop(raster_file, extent(boundary))
  
  
  #write
  output_dir <- stringr::str_c("~/lab-data/phenology/phenology_urban/NY/nlcd_tcc/clipped/nlcd_tcc_NY_",i,".tif")
  writeRaster(clipped_raster, output_dir, format = "GTiff", overwrite = TRUE)
}


# visualization
plot(lst_rs)
plot(boundary, col = "red")  # Add the boundary shapefile to the plot