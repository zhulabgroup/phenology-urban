##### Save the tif of seasonal average temperature #####
library(zhulabtools)
load_packages(c("tidyverse", "ggplot2", "sf", "lubridate", "terra"))
select <- dplyr::select

tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

for (city in cities[40:76]) {
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
  
  test0 <- UHAT_climate %>%
    distinct(year, season, lon, lat) |>
    group_by(year, season) |>
    summarise(n = n(), .groups = "drop") |>
    filter(year >= 2016)
  test <- raw_uhat %>%
    distinct(time, lon, lat) |>
    group_by(time) |>
    summarise(n = n(), .groups = "drop") |>
    filter(time >= as.Date("2016-01-01"))
  
  
  # output directory
  out_dir <- "~/lab-data/datasets/climate/UHAT/clean"
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
      trim() |>
      project("EPSG:4326")  # continuous variable → bilinear used by default
    
    # 6) write GeoTIFF
    fname <- file.path(
      out_dir,
      paste0(city, sprintf("_%s_%s.tif", seas, yr))
    )
    writeRaster(
      r_trim, fname,
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



library(terra)        # 读tif
library(dplyr)
library(ggplot2)
library(patchwork)    # wrap_plots
library(stringr)
library(scales)       # oob = squish
# 可选：更好看的连续色带
# install.packages("viridis")
library(viridis)

tif_files <- list.files(out_dir, pattern = "_winter_2022\\.tif$", full.names = TRUE)
stopifnot(length(tif_files) > 0)

rng <- do.call(rbind, lapply(tif_files, function(f) {
  m <- terra::minmax(terra::rast(f))
  data.frame(min = m[1, 1], max = m[2, 1])
}))

plot_one <- function(f) {
  r <- terra::rast(f)
  r1 <- r[[1]]
  df <- as.data.frame(r1, xy = TRUE, na.rm = TRUE)
  val_col <- names(df)[3]
  
  ggplot(df, aes(x = x, y = y, fill = .data[[val_col]])) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis(
      oob = scales::squish,
      name = val_col
    ) +
    labs(
      title = str_replace(basename(f), "_2020_summer.tif", ""),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}


city_plots <- lapply(tif_files, plot_one)
final_plot <- wrap_plots(city_plots, ncol = 8) &
  theme(plot.margin = margin(6, 6, 6, 6))

test <- Moab |>
  group_by(lat, lon) |>
  summarize(count = n())




