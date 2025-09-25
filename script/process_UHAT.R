library(zhulabtools)
load_packages(c("tidyverse", "ggplot2", "sf", "lubridate", "terra", "patchwork", "scales", "viridis", "leaflet"))
select <- dplyr::select
# output directory
out_dir <- "~/lab-data/datasets/climate/UHAT/clean"
tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

################## Save the seasonal file #################

for (city in cities) {
  message(paste0("Processing city: ", city))
  uhat_file <- paste0("~/lab-data/datasets/climate/UHAT/raw/", city, ".csv")
  if (!file.exists(uhat_file)) {
    message(paste0("No UHAT found for: ", city))
    next  # Skip the current iteration of the loop
  }
  raw_uhat <- read_csv(uhat_file)
  UHAT_climate <- raw_uhat |>
    mutate(season = case_when(
      month(time) %in% 3:5 ~ "spring",
      month(time) %in% 6:8 ~ "summer",
      month(time) %in% 9:11 ~ "fall",
      month(time) %in% c(12, 1, 2) ~ "winter"
    )) |>
    mutate(year = year(time),
           year = if_else(month(time) == 12, year + 1, year)) |>
    group_by(lon, lat, year, season) |>
    summarise(
      tmax = mean(t2m_pred_day, na.rm = TRUE),
      tmin = mean(t2m_pred_night, na.rm = TRUE),
      tavg = (tmax + tmin) / 2,
      .groups = "drop"
    ) |>
    # from K to C
    mutate(across(c(tmax, tmin, tavg), ~ . - 273.15)) |>
    select(lon, lat, year, season, tmax, tmin, tavg)

  
  # Match to centroids
  
  # Save the basic grid:
  dir_path <- "~/lab-data/datasets/climate/UHAT/grid"
  tifs <- list.files(
    dir_path,
    pattern = "^MYD11A1\\..*LST_Day_1km\\.tif$",
    full.names = TRUE
  )
  r_template <- terra::rast(tifs[1])
  
  seasons_to_do <- c("summer", "winter")
  
  # years present for those seasons
  years_to_do <- UHAT_climate %>%
    filter(season %in% seasons_to_do) %>%
    pull(year) %>%
    unique() %>%
    sort()

  # helper that builds one raster for (yr, seas)
  build_UHAT_tif <- function(yr, seas, out_dir) {
    # 1) slice data
    df <- UHAT_climate %>%
      filter(year == yr, season == seas) %>%
      transmute(lon, lat, temp = if ("tavg" %in% names(.)) tavg else temp)
    
    if (nrow(df) == 0) {
      message(sprintf("[skip] %s %s: no rows", yr, seas))
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
    fname <- file.path(
      out_dir,
      paste0(city, sprintf("_%s_%s.tif", seas, yr))
    )
    writeRaster(
      r_aligned, fname,
      overwrite = TRUE, NAflag = -9999,
      gdal = c("COMPRESS=LZW"), datatype = "FLT4S"
    )
    message(sprintf("[ok] wrote %s", fname))
    invisible(fname)
  }
  # run for all (year, season) pairs
  for (yr in years_to_do) {
    message(paste0("Processing year", yr))
    for (seas in seasons_to_do) {
      build_UHAT_tif(yr, seas, out_dir)
    }
  }
}


################## Check results ##################


plot_one <- function(city, out_dir, tree_dir,
                             yr = 2020, 
                             seas = "summer") {
  f <- file.path(out_dir, paste0(city,"_",seas,"_",yr,".tif"))
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
  
  UHAT_climate_point <-  read_csv(paste0("~/lab-data/datasets/climate/UHAT/raw/", city, ".csv")) |>
    mutate(season = case_when(
      month(time) %in% 3:5 ~ "spring",
      month(time) %in% 6:8 ~ "summer",
      month(time) %in% 9:11 ~ "fall",
      month(time) %in% c(12, 1, 2) ~ "winter"
    )) |>
    mutate(year = year(time),
           year = if_else(month(time) == 12, year + 1, year)) |>
    group_by(lon, lat, year, season) |>
    summarise(
      tmax = mean(t2m_pred_day, na.rm = TRUE),
      tmin = mean(t2m_pred_night, na.rm = TRUE),
      tavg = (tmax + tmin) / 2,
      .groups = "drop"
    ) |>
    # from K to C
    mutate(across(c(tmax, tmin, tavg), ~ . - 273.15)) |>
    select(lon, lat, year, season, tmax, tmin, tavg) |>
    filter(year == yr, season == seas) |>
    distinct(lon, lat) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
    st_transform(crs = crs(r_aligned))
  
  df <- as.data.frame(r_aligned, xy = TRUE, na.rm = TRUE)
  val_col <- names(df)[3]
  ggplot() +
    geom_sf(data = tree_sf, color = "black", fill = "black",
            size = 0.1, alpha = 0.8, shape = 16) +
    geom_sf(data = UHAT_climate_point, color = "red", fill = "red",
            size = 0.5, alpha = 0.7, shape = 16) +
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
plot_one(city = "Aurora", out_dir, tree_dir)

city_plots <- lapply(cities, plot_one)
final_plot <- wrap_plots(city_plots[15:17], ncol = 8) &
  theme(plot.margin = margin(6, 6, 6, 6))
final_plot

