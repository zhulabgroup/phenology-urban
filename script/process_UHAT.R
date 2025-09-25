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
library(viridis)
library(leaflet)


tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

plot_one <- function(city) {
  f <- paste0(out_dir, "/", city, "_winter_2020.tif")
  r <- terra::rast(f)
  r1 <- r[[1]]
  df <- as.data.frame(r1, xy = TRUE, na.rm = TRUE)
  val_col <- names(df)[3]
  tree_rds <- read_rds(paste0(tree_dir, city, ".rds"))
  tree_sf <- tree_rds |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
    
  ggplot() +
    geom_tile(data = df, aes(x = x, y = y, fill = .data[[val_col]])) +
    coord_equal() +
    scale_fill_viridis(oob = scales::squish, name = val_col) +
    geom_sf(data = tree_sf, color = "black", fill = "black",
            size = 0.6, alpha = 0.8, shape = 16) +
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
}


city_plots <- lapply(cities, plot_one)
final_plot <- wrap_plots(city_plots[15:17], ncol = 8) &
  theme(plot.margin = margin(6, 6, 6, 6))
final_plot


test <- Moab |>
  group_by(lat, lon) |>
  summarize(count = n())

plot_one_leaflet <- function(city, out_dir, tree_dir,
                             tif_suffix = "_winter_2020.tif",
                             raster_opacity = 0.8,
                             point_radius = 3,
                             point_color = "black",
                             point_opacity = 0.8) {
  
  f <- file.path(out_dir, paste0(city, tif_suffix))
  stopifnot(file.exists(f))
  r <- terra::rast(f)
  r1 <- r[[1]]
  r1_ll <- terra::project(r1, "EPSG:4326")   # 投影到经纬度
  rr <- raster::raster(r1_ll)                # 转为 raster 对象以便 addRasterImage
  
  # 2) 颜色映射（使用 viridis）
  vals <- raster::values(rr)
  pal <- colorNumeric(viridis(256), domain = vals, na.color = "transparent")
  
  # 3) 读树点并确保为 WGS84
  tree_rds_path <- file.path(tree_dir, paste0(city, ".rds"))
  stopifnot(file.exists(tree_rds_path))
  tree_rds <- readr::read_rds(tree_rds_path)
  
  # 如果列名不是 lon/lat，请改成你实际的列名
  tree_sf <- st_as_sf(tree_rds, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  # 4) 计算地图范围（用栅格范围）
  ex <- terra::ext(r1_ll)
  # xmin, xmax, ymin, ymax
  xmin <- ex[1]; xmax <- ex[2]; ymin <- ex[3]; ymax <- ex[4]
  
  # 5) 组装 Leaflet 地图
  leaflet() |>
    addProviderTiles(providers$OpenStreetMap) |>
    addRasterImage(rr, colors = pal, opacity = raster_opacity, project = FALSE) |>
    addCircleMarkers(
      data = tree_sf,
      lng = ~lon, lat = ~lat,
      radius = point_radius,
      fillColor = point_color, color = point_color,
      fillOpacity = point_opacity, opacity = point_opacity,
      stroke = FALSE
    ) |>
    fitBounds(lng1 = xmin, lat1 = ymin, lng2 = xmax, lat2 = ymax) |>
    addControl(html = htmltools::tags$div(
      style = "background: rgba(255,255,255,0.85); padding: 6px 8px; border-radius: 4px; font-size: 12px;",
      htmltools::tags$b(str_replace(basename(f), ".tif", ""))
    ),
    position = "topleft")
}
plot_one_leaflet(city = "Moab",out_dir,tree_dir)



