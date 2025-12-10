library(zhulabtools)
load_packages(c("tidyverse", "ggplot2", "sf", "lubridate", "terra", "patchwork", "scales", "viridis", "leaflet"))
select <- dplyr::select
# output directory
out_dir <- "~/phenology-urban/temp/"
tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

# Load the basic grid:
dir_path <- "~/lab-data/datasets/climate/UHAT/grid"
tifs <- list.files(
  dir_path,
  pattern = "^MYD11A1\\..*LST_Day_1km\\.tif$",
  full.names = TRUE
)
r_template <- terra::rast(tifs[1])

# helper that builds one raster for each date
build_UHAT_tif <- function(date, out_dir, type) {
  # # 0) check location
  # city_out_dir <- paste0(out_dir, "/", city, "/")
  # if (!dir.exists(city_out_dir)) {
  #   dir.create(city_out_dir, recursive = TRUE)
  # }
  # fname <- file.path(
  #   city_out_dir,
  #   paste0(sprintf("%s_%s.tif", date, type))
  # )
  # if (file.exists(fname)) {
  #   cat(sprintf("\n File already exists for city: %s, date: %s, type: %s. Skipping.\n", 
  #               city, date, type))
  #   return(FALSE)  # Skip the current iteration of the loop
  # }
  
  # 1) slice data
  df <- UHAT_climate %>%
    filter(time_str == date) %>%
    transmute(lon, lat, temp = if (type == "day") t2m_pred_day else t2m_pred_night)
  
  if (nrow(df) == 0) {
    message(sprintf("[skip] %s %s: no rows", date))
    return(invisible(NULL))
  }
  
  # 2) points -> template cells
  pts_ll  <- vect(df, geom = c("lon","lat"), crs = "EPSG:4326")
  pts_prj <- project(pts_ll, crs(r_template))  # to template CRS
  
  xy_prj  <- crds(pts_prj, df = TRUE)
  cell_id <- terra::cellFromXY(r_template, as.matrix(xy_prj))
  
  pts_with_cell <- df |>
    mutate(.cell = cell_id) |>
    filter(!is.na(.cell))
  
  if (nrow(pts_with_cell) == 0) {
    message(sprintf("[skip] %s %s: all points outside template", yr, seas))
    return(invisible(NULL))
  }
  
  # 3) aggregate duplicates per cell
  cell_values <- pts_with_cell |>
    group_by(.cell) |>
    summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop")
  
  # 4) fill new raster on template grid
  r_out <- rast(r_template)
  values(r_out) <- NA_real_
  r_out[cell_values$.cell] <- cell_values$temp
  
  # 5) trim NA borders and reproject to WGS84 lon/lat
  r_trim <- r_out |>
    trim()
  
  rx <- res(r_trim)[1]
  ry <- res(r_trim)[2]
  target_xmin <- min(xy_prj[,1]) - rx/2
  target_ymin <- min(xy_prj[,2]) - ry/2
  cur_ext <- ext(r_trim)
  cur_xmin <- cur_ext[1]
  cur_ymin <- cur_ext[3]
  dx_shift <- target_xmin - cur_xmin
  dy_shift <- target_ymin - cur_ymin
  r_aligned <- terra::shift(r_trim, dx = dx_shift, dy = dy_shift)
  crs(r_aligned) <- crs(r_template)  # ensure correct CRS
  
  # 6) write GeoTIFF
  
  # writeRaster(
  #   r_aligned, fname,
  #   overwrite = TRUE, NAflag = -9999,
  #   gdal = c("COMPRESS=LZW"), datatype = "FLT4S"
  # )
  # message(sprintf("[ok] extracted %s", date))
  # invisible(fname)
  # return(TRUE)
  
  return(r_aligned)
}

################## Save the seasonal file #################

for (city in cities[20:29]) {
  message(paste0("Processing city: ", city))
  uhat_file <- paste0("~/lab-data/datasets/climate/UHAT/raw/", city, ".csv")
  if (!file.exists(uhat_file)) {
    message(paste0("No UHAT found for: ", city))
    next  # Skip the current iteration of the loop
  }
  raw_uhat <- read_csv(uhat_file)
  
  UHAT_climate <- raw_uhat |>
    filter(time >= "2016-09-01") |>
    mutate(time_str = as.character(time)) |>
    mutate(time_str = gsub("-", "_", time_str)) |>
    # from K to C
    mutate(across(c(t2m_pred_day, t2m_pred_night), ~ . - 273.15)) |>
    select(lon, lat, time_str, t2m_pred_day, t2m_pred_night)
  
  day_dates_to_do <- UHAT_climate |>
    select(lon, lat, time_str, t2m_pred_day) |>
    filter(t2m_pred_day != "NA") |>
    distinct(time_str) |>
    pull()
  
  night_dates_to_do <- UHAT_climate |>
    select(lon, lat, time_str, t2m_pred_night) |>
    filter(t2m_pred_night != "NA") |>
    distinct(time_str) |>
    pull()
  
  # run for all dates
  # day
  nc_fname_day <- file.path(out_dir, paste0(city, "_UHAT_t2m_day.nc"))
  day_list <- list()
  for (date in day_dates_to_do) {
    message(paste0("Processing date", date))
    day_list[[date]] <-  build_UHAT_tif(date, out_dir, type = "day")
  }
  ext_list <- lapply(day_list, ext)
  union_ext <- ext(min(sapply(ext_list, function(e) e$xmin)), 
                   max(sapply(ext_list, function(e) e$xmax)),
                   min(sapply(ext_list, function(e) e$ymin)),
                   max(sapply(ext_list, function(e) e$ymax)))
  r_all_day <- lapply(day_list, function(r) {
    extend(r, union_ext)
  }) |> rast()
  time(r_all_day) <- as.Date(names(day_list), format = "%Y_%m_%d")
  terra::writeCDF(r_all_day, filename = nc_fname_day, overwrite = TRUE)
  
  # night
  nc_fname_night <- file.path(out_dir,paste0(city, "_UHAT_t2m_night.nc"))
  night_list <- list()
  for (date in night_dates_to_do) {
    message(paste0("Processing date", date))
    night_list[[date]] <- build_UHAT_tif(date, out_dir, type = "night")
  }
  ext_list <- lapply(night_list, ext)
  union_ext <- ext(min(sapply(ext_list, function(e) e$xmin)), 
                   max(sapply(ext_list, function(e) e$xmax)),
                   min(sapply(ext_list, function(e) e$ymin)),
                   max(sapply(ext_list, function(e) e$ymax)))
  r_all_night <- lapply(night_list, function(r) {
    extend(r, union_ext)
    }) |> rast()
  time(r_all_night) <- as.Date(names(night_list), format = "%Y_%m_%d")
  terra::writeCDF(r_all_night, filename = nc_fname_night, overwrite = TRUE)
}


################## Check results ##################


plot_one <- function(city, out_dir, tree_dir,
                             date = "", type = "day") {
  f <- file.path(out_dir, paste0(city,"/", date,"_", type,".tif"))
  stopifnot(file.exists(f))
  r_aligned <- terra::rast(f)
  dir_path <- "~/lab-data/datasets/climate/UHAT/grid"
  tifs <- list.files(
    dir_path,
    pattern = "^MYD11A1\\..*LST_Day_1km\\.tif$",
    full.names = TRUE
  )
  r_template <- terra::rast(tifs[1])
  crs(r_aligned) <- crs(r_template)  # ensure correct CRS
  
  tree_rds_path <- file.path(tree_dir, paste0(city, ".rds"))
  stopifnot(file.exists(tree_rds_path))
  tree_sf <- read_rds(tree_rds_path) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
    st_transform(crs(r_aligned))
  # tree_v <- vect(tree_sf_proj)
  # e <- ext(tree_v)
  # r <- crop(r, e, snap = "out")
  vals <- raster::values(r_aligned)
  pal <- colorNumeric(viridis(256), domain = vals, na.color = "transparent")

  UHAT_climate_point <-  read_csv(paste0("~/lab-data/datasets/climate/UHAT/raw/", city, ".csv"), show_col_types = FALSE) |>
    distinct(lon, lat) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
    st_transform(crs = crs(r_aligned))
  
  df <- as.data.frame(r_aligned, xy = TRUE, na.rm = TRUE)
  val_col <- names(df)[3]
  ggplot() +
    geom_sf(data = tree_sf, color = "black", fill = "black",
            size = 0.1, alpha = 0.8, shape = 16) +
    geom_sf(data = UHAT_climate_point, color = "red", fill = "red",
            size = 2, alpha = 0.7, shape = 16) +
    geom_tile(data = df, aes(x = x, y = y, fill = .data[[val_col]]), alpha = 0.8) +
    # coord_equal() +
    scale_fill_viridis(oob = scales::squish, name = val_col) +
    labs(
      title = str_replace(basename(f), ".tif", ""),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  # 
  # leaflet() |>
  #   addProviderTiles(providers$OpenStreetMap,
  #                    options = providerTileOptions(opacity = 0.9)) |>
  #   addRasterImage(r_aligned |> project("EPSG:4326"), colors = pal, opacity = 0.7, project = FALSE) |>
  #   addCircleMarkers(
  #     data = tree_sf,
  #     lng = ~lon, lat = ~lat,
  #     radius = 2,
  #     fillColor = "black", color = "black",
  #     fillOpacity = 0.4, opacity = 0.4,
  #     stroke = FALSE
  #   ) |>
  #   addCircleMarkers(
  #     data = UHAT_climate_point,
  #     lng = ~lon, lat = ~lat,
  #     radius = 3,
  #     fillColor = "red", color = "red",
  #     fillOpacity = 0.6, opacity = 0.6,
  #     stroke = FALSE
  #   ) |>
  #   addScaleBar(position = "bottomleft", options = scaleBarOptions(metric = TRUE, imperial = FALSE))
  
}
# plot_one(city = "Albuquerque", out_dir, tree_dir, date = "2013_04_02", type = "day")
# 
# city_plots <- lapply(cities, function(city) {
#   plot_one(city = city, out_dir = out_dir, tree_dir = tree_dir, yr = 2014, seas = "winter")
# })
# final_plot <- wrap_plots(city_plots, ncol = 11) &
#   theme(plot.margin = margin(6, 6, 6, 6))
# final_plot

