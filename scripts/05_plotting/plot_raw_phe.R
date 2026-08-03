pacman::p_load(terra, ggspatial)

source(here("scripts/04_stats_model/fit_spatial/prepare_model_input.R"))
source(here("scripts/05_plotting/plot_setting_city_map.R"))

plot_species_raw_phe <- function(species, color_palette, phe_type, color_limits = NULL, highlight_city = NULL, year){
  all_city_input <- read_rds(here(paste0("data/input/phenology/no_outlier/nooutlier_", phe_type, "_input.rds"))) |>
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
  
  
  inventory_centroids <- readRDS(here("data/input/tree_list/summary/inventory_centroids.rds")) |>
    left_join(species_city, by = c("site" = "city")) |>
    filter(!is.na(median_doy))
  boundary_us <- st_read(here("data/input/boundary/US/us_boundarydissolve.shp"), quiet = T) |> 
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

