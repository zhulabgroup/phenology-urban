library(zhulabtools)
load_packages(c("tidyverse", "ggplot2", "sf", "lubridate", "terra", "patchwork", "scales", "viridis", "leaflet"))
select <- dplyr::select
# output directory
out_dir <- "~/lab-data/datasets/climate/UHAT/clean"
tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))


################## Save the seasonal file #################


# helper that builds one raster for (yr, seas)
build_UHAT_tif <- function(yr, seas, out_dir, temp_type, UHAT_climate,
                           r_template, city, overwrite = FALSE) {
  # 1) slice data
  df <- UHAT_climate %>%
    filter(year == yr, season == seas) %>%
    transmute(lon, lat, temp = .data[[temp_type]])
  
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
    paste0(city, sprintf("_%s_%s_%s.tif", seas, yr, temp_type))
  )
  if (file.exists(fname) && !overwrite) {
    return(invisible(fname))
  }
  writeRaster(
    r_aligned, fname,
    overwrite = overwrite, NAflag = -9999,
    gdal = c("COMPRESS=LZW"), datatype = "FLT4S"
  )
  message(sprintf("[ok] wrote %s", fname))
  invisible(fname)
}

summarise_UHAT_seasons <- function(raw_uhat, season_type = "full_season") {
  season_type <- match.arg(season_type, c("full_season", "early_season", "later_season"))
  winter_months <- switch(season_type,
                          full_season = c(12, 1, 2), early_season = c(12, 1), later_season = c(1, 2)
  )
  summer_months <- switch(season_type,
                          full_season = c(6, 7, 8), early_season = c(6, 7), later_season = c(7, 8)
  )
  raw_uhat |>
    mutate(season = case_when(
      month(time) %in% winter_months ~ "winter",
      month(time) %in% summer_months ~ "summer",
      season_type == "full_season" & month(time) %in% 3:5 ~ "spring",
      season_type == "full_season" & month(time) %in% 9:11 ~ "fall",
      TRUE ~ NA_character_
    )) |>
    filter(!is.na(season)) |>
    mutate(year = year(time),
           year = if_else(month(time) == 12, year + 1, year)) |>
    group_by(lon, lat, year, season) |>
    summarise(
      tmax = mean(t2m_pred_day, na.rm = TRUE),
      tmin = mean(t2m_pred_night, na.rm = TRUE),
      tavg = (tmax + tmin) / 2,
      .groups = "drop"
    ) |>
    mutate(across(c(tmax, tmin, tavg), ~ . - 273.15)) |>
    select(lon, lat, year, season, tmax, tmin, tavg)
}

process_UHAT_city <- function(city, season_type = "full_season", overwrite = FALSE) {
  season_type <- match.arg(season_type, c("full_season", "early_season", "later_season"))
  out_dir <- paste0("~/lab-data/datasets/climate/UHAT/", switch(season_type,
                                                    full_season = "clean", early_season = "clean_early_season",
                                                    later_season = "clean_later_season"
  ))
  message(paste0("UHAT: Processing city: ", city))
  uhat_file <- paste0("~/lab-data/datasets/climate/UHAT/raw/", city, ".csv")
  if (!file.exists(uhat_file)) {
    stop("No UHAT found for: ", city)
  }
  raw_uhat <- read_csv(uhat_file)
  UHAT_climate <- summarise_UHAT_seasons(raw_uhat, season_type)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  # Match to centroids
  
  # Save the basic grid:
  dir_path <- "~/lab-data/datasets/climate/UHAT/grid"
  tifs <- list.files(
    dir_path,
    pattern = "^MYD11A1\\..*LST_Day_1km\\.tif$",
    full.names = TRUE
  )
  r_template <- terra::rast(tifs[1])
  
  seasons_to_do <- if (season_type == "full_season") {
    c("spring", "summer", "fall", "winter")
  } else {
    c("summer", "winter")
  }
  temp_type <- c("tavg")
  # years present for those seasons
  years_to_do <- UHAT_climate %>%
    filter(season %in% seasons_to_do) %>%
    pull(year) %>%
    unique() %>%
    sort()
  
  # run for all (year, season) pairs
  for (yr in years_to_do) {
    message(paste0("Processing year", yr))
    for (seas in seasons_to_do) {
      for (i in temp_type){
        build_UHAT_tif(yr, seas, out_dir, i, UHAT_climate, r_template, city, overwrite)
      }
    }
  }
}

if (sys.nframe() == 0L) {
  for (city in cities) {
    tryCatch(process_UHAT_city(city), error = function(e) message("Skipping ", city, ": ", e$message))
  }
}