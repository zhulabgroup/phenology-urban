pacman::p_load(terra, ggspatial)

source(here("script/04_stats_model/fit_spatial/prepare_model_input.R"))

fix_map_ratio <- function(xlim, ylim, ratio = 1.5) {
  xmid <- mean(xlim)
  ymid <- mean(ylim)
  
  xr <- diff(xlim)
  yr <- diff(ylim)
  
  current_ratio <- xr / yr
  
  if (current_ratio < ratio) {
    xr_new <- yr * ratio
    xlim <- c(xmid - xr_new / 2, xmid + xr_new / 2)
  } else {
    yr_new <- xr / ratio
    ylim <- c(ymid - yr_new / 2, ymid + yr_new / 2)
  }
  
  list(xlim = xlim, ylim = ylim)
}

search_boundary <- function(city) {
  boundary <- switch(
    city,
    "NewYork" = st_read(here("data/raw/Example_city/NY/boundary/nybb_dissolved.shp")),
    "DT" = st_read(here("data/raw/Example_city/DT/City_of_Detroit_Boundary/City_of_Detroit_Boundary.shp")),
    "Seattle" = st_read(here("data/raw/Example_city/ST/seattle_boundary/CityUGA.shp")) %>%
      filter(CITY_NM == "Seattle"),
    "Houston" = st_read(here("data/raw/Example_city/HT/HGAC_City_Boundaries/HGAC_City_Boundaries.shp")) %>%
      filter(NAME == "Houston") |>
      st_make_valid(),  # change to some distance
    "TP" = st_read(here("data/raw/Example_city/TP/Tampa_Planning_Districts/Tampa_Planning_Districts.shp")) %>%
      st_make_valid() |> st_union(),
    "Denver" = st_read(here("data/raw/Example_city/DV/DRCOG_Municipalities/DRCOG_Municipalities.shp")) %>%
      filter(city == "Denver"),
    "Austin" = st_read(here("data/raw/Example_city/AT/TxDOT_City_Boundaries/Cities.shp")) %>%
      filter(CITY_NM == "Austin"),
    "Rochester" = st_read(here("data/raw/Example_city/Rochester/City_of_Rochester_Border/City_of_Rochester_Border.shp")),
    "Stockton" = st_read(here("data/raw/Example_city/Stockton/City_Limits/geo_export_95d77263-4f44-4046-b5b9-23c00947fd1f.shp")) |>
      filter(name == "STOCKTON"),
    "LosAngeles" = st_read(here("data/raw/Example_city/LosAngeles/City_Boundary.shp")),
    "SanFrancisco" = st_read(here("data/raw/Example_city/SF/Planning_Neighborhood_Groups_Map/geo_export_a4e97979-9d5a-4474-8af3-073be4dffae4.shp")) |>
      st_make_valid() |> st_union(),
    "Knoxville" = st_read(here("data/raw/Example_city/Knoxville/TN_City_Boundaries/TN_City_Boundaries.shp")) %>%
      filter(NAME == "Knoxville"),
    "WashingtonDC" = st_read(here("data/raw/Example_city/WashingtonDC/DC_Boundary/DC_Boundary.shp")),
    "Columbus" = st_read(here("data/raw/Example_city/Columbus/Corporate_Boundary/Corporate_Boundary.shp")),
    "Arlington" = st_read(here("data/raw/Example_city/Arlington/Municipal_Boundary/Municipal_Boundary.shp"))|> 
      filter(MUNICIPALN == "ARLINGTON"),
    "Sacramento" = st_read(here("data/raw/Example_city/Sacramento/City_and_County_Boundary/City_and_County_Boundary.shp")) |>
      filter(CITY_NAME == "SACRAMENTO") |>
      st_make_valid(),
    "Atlanta" = st_read(here("data/raw/Example_city/Atlanta/Quadrants/Quadrants.shp")) |>
      st_make_valid() |> st_union()  |>
      nngeo::st_remove_holes(),
    "NewOrleans" = st_read(here("data/raw/Example_city/NewOrleans/Orleans Parish Boundary_20260416/geo_export_2201ea0f-3a1c-4066-aec3-a8e995d07a7a.shp")),
    "Portland" = st_read(here("data/raw/Example_city/Portland/City_Boundaries/City_Boundaries.shp")) |>
      filter(CITYNAME == "Portland") |>
      st_make_valid(),
    "Bloomington" = st_read(here("data/raw/Example_city/Bloomington/Bloomington Municipal Boundary_20260421/geo_export_951c8648-2ab5-41e2-820c-1a4d268e2537.shp")),
    "Norfolk" = st_read(here("data/raw/Example_city/Norfolk/City_Limit/City_Limit_-_City_of_Norfolk.shp")),
    stop("Invalid city abbreviation.")
  )
  boundary <- st_transform(boundary, crs = 4326)  # ensure WGS84
  return(boundary)
}

plot_species_raw_phe <- function(species, color_palette, phe_type, color_limits = NULL, highlight_city = NULL, year){
  all_city_input <- read_rds(here(paste0("release/input/phenology/no_outlier/nooutlier_", phe_type, "_input.rds"))) |>
    filter(year == !!year)
  phe_title <- case_when(
    phe_type == "SOS" ~ "Start-of-season",
    phe_type == "EOS" ~ "End-of-season"
  )
  
  species_city <- all_city_input |>
    filter(canonicalname == species)|>
    dplyr::group_by(city) |>
    dplyr::summarise(
      median_doy = median(doy, na.rm = TRUE),
      sd_doy = sd(doy, na.rm = TRUE),
      IQR_d = quantile(doy, 0.75, na.rm = TRUE) - quantile(doy, 0.25, na.rm = TRUE),
      lower_IQR_d = max(quantile(doy, 0.25, na.rm = TRUE), min(doy, na.rm = TRUE)),
      upper_IQR_d = min(quantile(doy, 0.75, na.rm = TRUE), max(doy, na.rm = TRUE)),
      .groups = "drop") |>
    arrange(median_doy) |>
    mutate(city = factor(city, levels = city))
  
  integ_city <- all_city_input |>
    filter(canonicalname == species)|>
    dplyr::group_by(city) |>
    dplyr::summarise(
      median_doy = median(doy, na.rm = TRUE),
      .groups = "drop") |>
    summarise(
      median_doy_cities = median(median_doy, na.rm = TRUE),
      sd_doy_cities = sd(median_doy, na.rm = TRUE),
      IQR_d = quantile(median_doy, 0.75, na.rm = TRUE) - quantile(median_doy, 0.25, na.rm = TRUE),
      lower_IQR_d = max(quantile(median_doy, 0.25, na.rm = TRUE), min(median_doy, na.rm = TRUE)),
      upper_IQR_d = min(quantile(median_doy, 0.75, na.rm = TRUE), max(median_doy, na.rm = TRUE)),
      .groups = "drop")
  
  sp_boxplot <- ggplot() +
    geom_rect(
      data = integ_city,
      aes(xmin = 0.5, xmax = length((all_city_input |>
                                       filter(canonicalname == species) |>
                                       distinct(city) |>
                                       pull())),
          ymin = lower_IQR_d, ymax = upper_IQR_d),
      fill = "grey90", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_abline(slope = 0, intercept = integ_city$median_doy_cities, linetype = "dotted", color = "grey80", size = 1) +
    geom_errorbar(data = species_city,
                  aes(x = city, ymin = lower_IQR_d, ymax = upper_IQR_d), color = "grey60",
                  width = 0.2, linewidth = 0.5, alpha = 0.7, inherit.aes = FALSE) +
    geom_errorbar(data = species_city |> dplyr::filter(city %in% highlight_city),
      aes(x = city, ymin = lower_IQR_d, ymax = upper_IQR_d), color = "grey50",
      width = 0.2, linewidth = 1.2, alpha = 1, inherit.aes = FALSE ) +
    geom_point(data = species_city,
               aes(x = city, y = median_doy, color = median_doy), size = 2,
               alpha = 0.5, inherit.aes = FALSE) +
    scale_color_viridis_c(name = paste0(phe_title, "\n(DOY)"),
                           limits = color_limits,
                          breaks = color_limits,
                          labels = round(color_limits, 0),
                          guide = guide_colorbar(
                            direction = "vertical",
                            label.position = "right"
                          )) +
    coord_flip(ylim = case_when(
      phe_type == "SOS" ~ c(0, 150),
      phe_type == "EOS" ~ c(200, 350)
    )) +
    labs(
      x = "City"
    ) +
    theme(
      plot.title = ggtext::element_markdown()
    ) +
    labs(y = case_when(
      phe_type == "SOS" ~ "Start-of-season (DOY)",
      phe_type == "EOS" ~ "End-of-season (DOY)"
    )) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 15),
      legend.position = "right",
      legend.key.height = unit(0.3, "cm"),   
      legend.key.width  = unit(0.5, "cm"),   
      legend.title = element_text(size = 14, margin = margin(r = 10, b=10)), 
      legend.text  = element_text(size = 12),
      legend.spacing.y = unit(0.1, "cm"),
      panel.background = element_rect(fill = "white", color = "black"),
      # Grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Facet/strip placement
      strip.placement = "outside"
    )
  
  
  inventory_centroids <- readRDS(here("release/input/tree_list/summary/inventory_centroids.rds")) |>
    left_join(species_city, by = c("site" = "city")) |>
    filter(!is.na(median_doy))
  boundary_us <- st_read(here("release/input/boundary/US/us_boundarydissolve.shp"), quiet = T) |> 
    st_transform(4326)
  
  doy_map <- ggplot() +
    geom_sf(data = boundary_us, fill = "grey90", color = "white") +
    geom_sf(data = inventory_centroids |> dplyr::filter(site %in% highlight_city), 
            fill = NA, color="grey50", size = 5,stroke=1.5, alpha = 1, shape = 21) +
    geom_sf(data = inventory_centroids, aes(color = median_doy), size = 5, alpha = 0.6) +
    coord_sf(expand = FALSE) +
    labs(x = "Longitude", y = "Latitude", title = paste0("<i>", species, "</i>")) +
    scale_color_viridis_c(name = paste0("City median\n", phe_title),
                          limits = color_limits,
                          breaks = color_limits,
                          labels = round(color_limits, 0),
                          guide = guide_colorbar(
                            direction = "horizontal",
                            label.position = "bottom"
                          )) +
    scale_radius(name = paste0("City SD\n", phe_title),
                 range = c(1, 20)) +
    theme_minimal() +
    theme(
      plot.title = ggtext::element_markdown(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 15),
      axis.ticks   = element_line(linewidth = 0.2),
      title = element_text(size = 16),
      legend.position = "none",
      legend.key.height = unit(0.3, "cm"),   
      legend.key.width  = unit(0.5, "cm"),   
      legend.title = element_text(size = 14, margin = margin(r = 10, b=10)), 
      legend.text  = element_text(size = 12),
      legend.spacing.y = unit(0.1, "cm"),
    )
  
  return(list(sp_boxplot = sp_boxplot, doy_map = doy_map))
}

plot_phe_one_city <- function(city = "NewYork", species = "Acer rubrum", phe_type = "SOS", year = 2022, color_limits = NULL){
  data_input <- prepare_data(city = city, phe_type = phe_type) |>
    filter(canonicalname == !!species) |>
    filter(year == !!year)
  phe_metric <- "doy"
  x <- data_input[[phe_metric]]
  q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- diff(q)
  lo  <- q[1]
  hi  <- q[2]
  city_input <- data_input  |>
    filter(is.na(.data[[phe_metric]]) |
             dplyr::between(.data[[phe_metric]], lo, hi))
  
  city_boundary <- search_boundary(city)
  
  tree_location <- city_input |>
    filter(!(lon == "NaN") & !(lat == "NaN")) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F)
  
  tree_location <- tree_location|>
    st_transform(crs(city_boundary)) |>
    st_intersection(city_boundary)
  bbox_tree <- st_bbox(tree_location)
  bbox_city <- st_bbox(city_boundary)
  
  xlim <- c(bbox_city["xmin"], bbox_city["xmax"])
  ylim <- c(bbox_city["ymin"], bbox_city["ymax"])
  lims <- fix_map_ratio(xlim, ylim, ratio = 1.3)
  
  if (phe_type == "SOS") {
    legend_phe_name <- "Start-of-season\n(DOY)"
    legend_temp_name <- "Mean winter temperature\n(°C)"
  } else if (phe_type == "EOS") {
    legend_phe_name <- "End-of-season\n(DOY)"
    legend_temp_name <- "Mean summer temperature\n(°C)"
  } else {
    stop("Invalid phe_type. Use 'SOS' or 'EOS'.")
  }
  map_phe <- ggplot() +
    geom_sf(data = city_boundary, fill = "grey90", color = "black") +
    geom_sf(data = tree_location, aes(color = doy),  size = 2, alpha = 0.6) +
    scale_color_viridis_c(name = legend_phe_name,
                          limits = color_limits,
                          breaks = color_limits,
                          labels = round(color_limits, 0),
                          guide = guide_colorbar(
                            direction = "horizontal",
                            label.position = "bottom"
                          )) +
    scale_y_continuous(
      n.breaks = 3
    ) +
    scale_x_continuous(
      n.breaks = 3
    ) +
    labs(title = paste0("<i>", species, "</i>"))+
    theme_minimal() +
    coord_sf(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE,clip = "off"   # important
    ) +
    annotation_scale(
      location = "tr",
      width_hint = 0.3,
      text_cex = 1.2,
      unit_category = "metric",
      pad_x = unit(0.1, "cm"),
      pad_y = unit(-0.4, "cm")
    ) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_blank(),
      title = element_text(size = 5, color="white"),
      legend.position = "none",
      # Panel background and border
      panel.background = element_rect(fill = "white", color = "black"),
      # Grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Facet/strip placement
      strip.placement = "outside"
    )
  return(map_phe)
}

