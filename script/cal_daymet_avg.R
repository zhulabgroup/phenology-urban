library(terra)
library(sf)
# download the Daymet data
library(daymetr)  # for getting Daymet data
library(ncdf4)  # for processing .nc data
library(janitor)  # for detecting cols with all records empty
library(gstat)  # for interpolation
library(sp)
library(sf)
library(raster)


city = "NY"
folder.raw <- paste0("~/phenology-urban/data/raw/", city, "/Daymet/daily2016-2024")
download_Daymet <- TRUE

################# based on the boundary of the city #################
us_boundary <- sf::read_sf("~/phenology-urban/data/raw/NY/boundary/nybb_dissolved.shp")
us_boundary <- us_boundary %>% st_transform(4326)

################# larege area #################

# if (download_Daymet) {
#   bbox <- as.vector(terra::ext(us_boundary))
#   n_rows <- 5
#   n_cols <- 10
#   x_seq <- seq(bbox["xmin"], bbox["xmax"], length.out = n_cols + 1)
#   y_seq <- seq(bbox["ymin"], bbox["ymax"], length.out = n_rows + 1)
#   for (i in seq_len(length(x_seq) - 1)) {
#     for (j in seq_len(length(y_seq) - 1)) {
#       sub_bbox <- c(
#         ymax = y_seq[j + 1],
#         ymin = y_seq[j],
#         xmax = x_seq[i + 1],
#         xmin = x_seq[i]
#       )
#       folder.sub <- paste0(folder.raw, "/region_", i, "_", j)
#       if (!dir.exists(folder.sub)) {
#         dir.create(folder.sub, recursive = TRUE)
#       }
#       download_daymet_ncss(location = sub_bbox[c("ymax", "xmin", "ymin", "xmax")],
#                            start = 2016, end = 2024,
#                            frequency = "daily",
#                            param = c("srad"),
#                            path = folder.sub)
#     }
#   }
# }

############# city area #############

if (download_Daymet){
  # get LST extent
  bbox <- as.vector(terra::ext(us_boundary))
  # download Daymet data
  download_daymet_ncss(location = bbox[c('ymax', 'xmin', 'ymin', 'xmax')],
                       start = 2016, end = 2023,
                       frequency = "daily",
                       param = c("tmin", "tmax", "prcp"),
                       path=folder.raw)
}

############# neon area #############
library(ggplot2)
library(sf)
library(dplyr)

# # 计算 bbox
# bbox <- neon_tree_location %>%
#   filter(site == "BART") %>%
#   vect() %>%
#   terra::buffer(width = 1000) %>%
#   terra::ext() %>%
#   as.vector()  # 转换为数值向量 c(xmin, xmax, ymin, ymax)
# 
# # 将 bbox 转换为 sf 矩形
# bbox_sf <- st_as_sf(st_sfc(st_polygon(list(matrix(
#   c(bbox[1], bbox[3],  # xmin, ymin
#     bbox[2], bbox[3],  # xmax, ymin
#     bbox[2], bbox[4],  # xmax, ymax
#     bbox[1], bbox[4],  # xmin, ymax
#     bbox[1], bbox[3]
#   ), ncol = 2, byrow = TRUE))),
#   crs = st_crs(neon_tree_location)))
# 
# points = neon_tree_location %>%
#   filter(site == "BART")
# 
# # 绘图
# ggplot() +
#   # geom_raster(data = clim_df_transformed, aes(x = X, y = Y, fill = value)) +
#   geom_sf(data = bbox_sf, fill = NA, color = "red", size = 1) +  # 绘制 bbox
#   geom_sf(data = points, color = "blue", size = 2) +  # 绘制点数据
#   theme_minimal() +
#   labs(title = paste("Bounding Box and Points for Site"))




for (i in names(neon_site_list)){
  bbox <-  neon_tree_location %>%
    filter(site == i) %>%
    vect() %>%
    terra::buffer(width = 1000) %>%
    terra::ext() %>%
    as.vector()
  points = neon_tree_location %>%
    filter(site == i) 
  folder.raw <- paste0("~/phenology-urban/data/raw/CONUS/Daymet/NEON/", i) 
  if (!dir.exists(folder.raw)) {
    dir.create(folder.raw, recursive = TRUE)
  }
  download_daymet_ncss(location = bbox[c('ymax', 'xmin', 'ymin', 'xmax')],
                       start = 2016, end = 2023,
                       frequency = "daily",
                       param = c("tmin", "tmax", "prcp", "srad"),
                       path=folder.raw)
}

# 
# 
# folder.raw <- "data/raw/CONUS/Daymet/annual1980-2023/"
# output_folder <- "~/phenology-urban/data/proc/daymet/conus/"
# 
# 
# calculate_and_save_tavg_tif <- function(year) {
#   tmin_file <- paste0(folder.raw, "tmin_annavg_", year, "_ncss.nc")
#   tmax_file <- paste0(folder.raw, "tmax_annavg_", year, "_ncss.nc")
#   tmin.nc <- terra::rast(tmin_file)
#   tmax.nc <- terra::rast(tmax_file)
#   usboun_reproj <- sf::st_transform(us_boundary, crs = sf::st_crs(terra::crs(tmax.nc, proj = TRUE)))
#   
#   tmax <- terra::project(terra::mask(terra::crop(tmax.nc, terra::ext(usboun_reproj)), usboun_reproj), terra::crs(us_boundary))
#   tmin <- terra::project(terra::mask(terra::crop(tmin.nc, terra::ext(usboun_reproj)), usboun_reproj), terra::crs(us_boundary))
#   
#   tavg_ann <- (tmax[[1]] + tmin[[1]]) / 2
#   output_file <- paste0(output_folder, "tavg_ann_", year, ".tif")
#   terra::writeRaster(tavg_ann, output_file, overwrite = TRUE)
# }
# 
# years <- 1980:2023
# lapply(years, calculate_and_save_tavg_tif)




############# clean srad #############

data_dir <- paste0("~/phenology-urban/data/raw/", city, "/Daymet/daily2016-2024/")
file_list <- list.files(data_dir, pattern = "srad_daily_\\d{4}_ncss\\.nc$", full.names = TRUE)

result <- NULL

# points_in_buffer <- read_csv(paste0("~/phenology-urban/data/proc/urban/", city, "/tree_WU_500_buffer_PS.csv"), show_col_types = FALSE)
# points <- points_in_buffer %>%
#   group_by(id) %>%
#   summarise(geometry = first(geometry))

points <- neon_tree_location # urban_tree_location

coords <- do.call(rbind, lapply(points$geometry, function(x) {
  as.numeric(unlist(strsplit(gsub("c\\(|\\)", "", x), ", ")))
}))
points_vect <- vect(coords, crs = "EPSG:4326")

for (file in file_list) {
  srad_nc <- terra::rast(file)
  extracted_values <- terra::extract(srad_nc, points_vect) %>%
    dplyr::select(-ID)
  colnames(extracted_values) <- as.character(srad_nc@ptr[["time"]])
  temp_result <- cbind(points["id"], extracted_values) %>%
    pivot_longer(
      cols = -id,
      names_to = "time",
      values_to = "value"
    ) %>%
    mutate(
      date = as.Date(as.POSIXct(as.numeric(time), origin = "1970-01-01", tz = "UTC"))
    )
  result <- bind_rows(result, temp_result)
}

saveRDS(result, paste0("~/phenology-urban/data/raw/", city, "/Daymet/daily2016-2024/srad_daily_2016-2023.rds"))


############# clean srad for prediction #############
city = "NY"
data_dir <- paste0("~/phenology-urban/data/raw/", city, "/Daymet/daily2016-2024/")
file_list <- list.files(data_dir, pattern = "srad_daily_(2022|2023)_ncss\\.nc$", full.names = TRUE)

result <- NULL

# points_in_buffer <- read_csv(paste0("~/phenology-urban/data/proc/urban/", city, "/tree_WU_500_buffer_PS.csv"), show_col_types = FALSE)
# points <- points_in_buffer %>%
#   group_by(id) %>%
#   summarise(geometry = first(geometry))

points <- read_csv(paste0("~/urban-cooling/data/raw/WU/", city, "/location.csv"), show_col_types = FALSE) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)
points_vect <- vect(points)
points <- st_drop_geometry(points)
for (file in file_list) {
  srad_nc <- terra::rast(file)
  extracted_values <- terra::extract(srad_nc, points_vect) %>%
    dplyr::select(-ID)
  colnames(extracted_values) <- as.character(time(srad_nc))
  temp_result <- cbind(points["Site"], extracted_values) %>%
    pivot_longer(
      cols = -Site,
      names_to = "time",
      values_to = "value"
    )
  result <- bind_rows(result, temp_result)
}

saveRDS(result, paste0("~/phenology-urban/data/raw/", city, "/Daymet/daily2016-2024/srad_wu_2023.rds"))
