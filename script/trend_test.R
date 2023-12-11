library(trend)
library(raster)

#read raster
raster_list <- lapply(1:10, function(i) raster(paste0("~/phenology-urban/data/raw/NY/LST_L8/year_mean/", i+2012, ".tif")))

#prepare the blank tif
mk_result <- raster(raster_list[[1]])  # 
mk_result[] <- NA  # 

#mk_test for every pixel
for (i in 1:ncell(mk_result)) {
  values_at_pixel <- sapply(raster_list, function(x) getValues(x)[i])
  if (all(!is.na(values_at_pixel))) {
    mk_test_result <- mk.test(values_at_pixel)
    mk_result[i] <- ifelse(mk_test_result$p.value < 0.05, mk_test_result$z, 0)
  }
}

writeRaster(mk_result, "path/to/mk_result.tif", format = "GTiff", overwrite = TRUE)
