##### Save the tif of seasonal average temperature #####
library(zhulabtools)
load_packages(c("tidyverse", "ggplot2", "sf", "spBayes", "INLA", "lubridate", "terra"))
select <- dplyr::select

t2m_prediction_NewYork_Newark_NY <- read_csv("~/lab-data/datasets/climate/UHAT/raw/t2m_prediction_NewYork_Newark_NY.csv")
UHAT_climate <- t2m_prediction_NewYork_Newark_NY |>
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

# test0 <- UHAT_climate %>%
#   distinct(year, season, lon, lat) |>
#   group_by(year, season) |>
#   summarise(n = n(), .groups = "drop")
# 
# 
# test <- t2m_prediction_NewYork_Newark_NY %>%
#   distinct(time, lon, lat) |>
#   group_by(time) |>
#   summarise(n = n(), .groups = "drop")


# output directory
out_dir <- "~/lab-data/datasets/climate/UHAT/clean"
# helper that builds one raster for (yr, seas)
build_UHAT_tif <- function(yr, seas) {
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
    sprintf("NY_%s_%s.tif", seas, yr)
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
  for (seas in seasons_to_do) {
    build_UHAT_tif(yr, seas)
  }
}



