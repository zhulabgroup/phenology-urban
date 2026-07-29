pacman::p_load(terra)
# output directory
out_dir <- here("~/lab-data/datasets/climate/UHAT/clean")
tree_dir <- here("release/input/tree_list/")
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

################## Save the seasonal file #################


# helper that builds one raster for (yr, seas)
build_UHAT_tif <- function(yr, seas, out_dir, temp_type) {
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
  writeRaster(
    r_aligned, fname,
    overwrite = TRUE, NAflag = -9999,
    gdal = c("COMPRESS=LZW"), datatype = "FLT4S"
  )
  message(sprintf("[ok] wrote %s", fname))
  invisible(fname)
}

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
  
  seasons_to_do <- c("spring", "summer", "fall", "winter")
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
        build_UHAT_tif(yr, seas, out_dir, i)
      }
    }
  }
}
