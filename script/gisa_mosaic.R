library(terra)
library(stringr)

raster_list <- list()
for (lon in seq(70, 125, 5)) {
  for (lat in seq(25, 50, 5)) {
    file_str <- str_c("/Volumes/seas-zhukai/datasets/imagery/GlobalImperviousSurface2.0/images/GISA_v02_1972_2019_", lat, "N_", lon, "W.tif")
    # Append the raster to the list
    if (file.exists(file_str)) {
      print(file_str)
      single_gisa <- rast(file_str)
      raster_list[[length(raster_list) + 1]] <- single_gisa
    } else {
    }
    
  }
}

# Mosaic all the resampled rasters
us_gisa_file <- do.call(mosaic, raster_list)
us_gisa_file
# View the resulting mosaic raster
plot(us_gisa_file)
