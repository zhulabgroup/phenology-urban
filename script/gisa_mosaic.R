library(terra)
library(stringr)

raster_list <- list()
for (lon in seq(70, 125, 5)) {
  for (lat in seq(25, 50, 5)) {
    file_str <- str_c("/Volumes/seas-zhukai/datasets/imagery/GlobalImperviousSurface2.0/images/GISA_v02_1972_2019_", lat, "N_", lon, "W.tif")
    # Append the raster to the list
    if (file.exists(file_str)) {
      #print(file_str)
      single_gisa <- rast(file_str)
      raster_list[[length(raster_list) + 1]] <- single_gisa
    } else {
    }
    
  }
}


# Mosaic all the resampled rasters
us_gisa_file1 <- do.call(mosaic, raster_list[1:16])
writeRaster(us_gisa_file1, "~/phenology-urban/data/proc/conus/us_gisa1.tif", overwrite = TRUE)
us_gisa_file2 <- do.call(mosaic, raster_list[17:32])
writeRaster(us_gisa_file2, "~/phenology-urban/data/proc/conus/us_gisa2.tif", overwrite = TRUE)
us_gisa_file3 <- do.call(merge, raster_list[33:48])
writeRaster(us_gisa_file3, "~/phenology-urban/data/proc/conus/us_gisa3.tif", overwrite = TRUE)
us_gisa_file4 <- do.call(merge, raster_list[49:64])
writeRaster(us_gisa_file4, "~/phenology-urban/data/proc/conus/us_gisa4.tif", overwrite = TRUE)

# Merge the four rasters
single_gisa <- merge(us_gisa_file1, us_gisa_file2, us_gisa_file3, us_gisa_file4)
writeRaster(single_gisa, "~/phenology-urban/data/proc/conus/us_gisa.tif", overwrite = TRUE)
