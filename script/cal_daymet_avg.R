library(terra)
library(sf)

folder.raw <- "data/raw/CONUS/Daymet/annual1980-2023/"
output_folder <- "~/phenology-urban/data/proc/daymet/conus/"

calculate_and_save_tavg_tif <- function(year) {
  tmin_file <- paste0(folder.raw, "tmin_annavg_", year, "_ncss.nc")
  tmax_file <- paste0(folder.raw, "tmax_annavg_", year, "_ncss.nc")
  tmin.nc <- terra::rast(tmin_file)
  tmax.nc <- terra::rast(tmax_file)
  usboun_reproj <- sf::st_transform(us_boundary, crs = sf::st_crs(terra::crs(tmax.nc, proj = TRUE)))
  
  tmax <- terra::project(terra::mask(terra::crop(tmax.nc, terra::ext(usboun_reproj)), usboun_reproj), terra::crs(us_boundary))
  tmin <- terra::project(terra::mask(terra::crop(tmin.nc, terra::ext(usboun_reproj)), usboun_reproj), terra::crs(us_boundary))
  
  tavg_ann <- (tmax[[1]] + tmin[[1]]) / 2
  output_file <- paste0(output_folder, "tavg_ann_", year, ".tif")
  terra::writeRaster(tavg_ann, output_file, overwrite = TRUE)
}

years <- 1980:2023
lapply(years, calculate_and_save_tavg_tif)