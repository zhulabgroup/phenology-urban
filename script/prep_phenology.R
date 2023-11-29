# Clip phenology index(sos,eos from Li, 2019) based on specific boundary and visualization


library(sf)
library(raster)

# Read shapefile
boundary <- st_read("~/lab-data/phenology/phenology_urban/NY/boundary/nybb_dissolved.shp")
raster_brick <- raster::brick(paste("/home/jializ/phenology-urban/data/raw/Landsat/US_uCluster_142_SOS_1985_2015.tif",
                                    sep = ""))

for (i in 30:31) {
  # Read raster file
  raster_file <- raster_brick[[i]]
  # same crs
  boundary <- st_transform(boundary, crs(raster_file))
  # clip
  clipped_raster <- crop(raster_file, extent(boundary))
  #write
  output_dir <- stringr::str_c("~/lab-data/phenology/phenology_urban/NY/phenology_Li/SOS_",i+1984,".tif")
  writeRaster(clipped_raster, output_dir, format = "GTiff", overwrite = TRUE)
}

# visualization
# Set breaks for color scheme
breaks <- seq(0, 366, by = 1)
# Set a colorscheme by assigning colors to breaks
colors <- colorRampPalette(c("white", "green", "red"))(length(breaks) - 1)

plot(clipped_raster,breaks = breaks, col = colors)
plot(boundary, add = TRUE, col = "red")  # Add the boundary shapefile to the plot