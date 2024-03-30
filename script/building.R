library(raster)
tif_files <- list.files(path="~/phenology-urban/data/raw/NY/dem", pattern="\\.tif$", full.names=TRUE)

library(sf)
polygons <- st_read("~/phenology-urban/data/raw/NY/NY_building_footprint/geo_export_56deb8a8-0463-4f00-8f45-99e751e3b8e1.shp")
polygons$dem <- polygons$groundelev + polygons$heightroof
building_filtered <- polygons[!is.na(polygons$dem) & polygons$dem != 0, ]



for (file in tif_files) {
  raster_data <- raster(file)
  raster_extent <- extent(raster_data)
  raster_proj <- projection(raster_data)
  raster_resolution <- res(raster_data)
  
  raster_output <- raster(raster_extent, res = raster_resolution, crs = raster_proj)
  outpath <- sub("be_", "bld_be_", file)
  
  rasterized <- rasterize(polygons, raster_output, field = "dem", fun='last', update=TRUE, filename=outpath)
  
}


