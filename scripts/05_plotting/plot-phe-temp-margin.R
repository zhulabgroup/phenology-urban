source(here("scripts/05_plotting/plot-setting-city-map.R"))
source(here("scripts/04_stats_model/fit_spatial/prepare_model_input.R"))
pacman::p_load(terra, ggspatial)

plot_phe_temp_margin <- function(city = "NewYork", species = "Acer rubrum", phe_type = "SOS", year = 2022){
  model_input_sf <-  prepare_data(city = city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)|>
    filter(canonicalname == species) |>
    filter(year == !!year) |>
    st_drop_geometry()
  phe_metric <- "doy"
  x <- model_input_sf[[phe_metric]]
  q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- diff(q)
  lo  <- q[1] - 1.5 * iqr
  hi  <- q[2] + 1.5 * iqr
  city_input <- model_input_sf  |>
    filter(is.na(.data[[phe_metric]]) |
             dplyr::between(.data[[phe_metric]], lo, hi))
  
  city_boundary <- search_boundary(city)
  
  tree_location <- city_input |>
    filter(!(lon == "NaN") & !(lat == "NaN")) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F)
  
  if (phe_type == "SOS") {
    season <- "winter"
  } else if (phe_type == "EOS") {
    season <- "summer"
  } else {
    stop("Invalid phe_type. Use 'SOS' or 'EOS'.")
  }
  file <- paste0(here("release/input/climate/UHAT/clean/"),city, "_", season, "_", year, ".tif")
  temp_tif <- terra::rast(file)
  tifs <- list.files(
    here("release/input/climate/UHAT/grid"),
    pattern = "^MYD11A1\\..*LST_Day_1km\\.tif$",
    full.names = TRUE
  )
  r_template <- terra::rast(tifs[1])
  terra::crs(temp_tif) <- terra::crs(r_template)  # ensure correct CRS
  tree_location <- tree_location|>
    st_transform(terra::crs(temp_tif))
  
  clipped_over <- temp_tif |>
    # crop(temp_tif, vect(tree_location)) |>
    terra::project(terra::crs(city_boundary), res = 0.009)  # transform to WGS84 for plotting
  city_v <- vect(city_boundary)
  clipped <- crop(clipped_over, city_v) |> mask(city_v)
  temp_df <- as.data.frame(clipped, xy = TRUE) |>
    rename_with(~ "value", .cols = 3)
  
  tree_location <- tree_location|>
    st_transform(terra::crs(city_boundary))  
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
    geom_sf(data = tree_location, aes(color = doy), size = 0.5) +
    scale_color_viridis_c(name = legend_phe_name,
                          breaks = range(tree_location$doy, na.rm = TRUE),
                          labels = round(range(tree_location$doy, na.rm = TRUE), 0),
                          guide = guide_colorbar(
                            direction = "horizontal",
                            label.position = "bottom"
                          )) +
    scale_y_continuous(
      n.breaks = 5 
    ) +
    scale_x_continuous(
      n.breaks = 3,
      labels = function(x) {
        lab <- ifelse(
          x < 0,
          paste0(sprintf("%.2f", abs(x)), "°W"),
          paste0(sprintf("%.2f", x), "°E")
        )
        lab[seq_along(lab) > 2] <- ""
        lab
      }
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    coord_sf(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE,
             clip = "off"   # important
    ) +
    annotation_scale(
      location = "tr",
      width_hint = 0.4,
      text_cex = 1.3,
      unit_category = "metric",
      pad_y = unit(-0.7, "cm")
    ) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.position = "bottom",
      legend.key.height = unit(0.3, "cm"),   
      legend.key.width  = unit(0.5, "cm"),   
      legend.title = element_text(size = 14, margin = margin(r = 10, b=10)), 
      legend.text  = element_text(size = 12),
      legend.spacing.y = unit(0.1, "cm"),
    )
  
  map_temp <- ggplot() +
    geom_sf(data = city_boundary, fill = "grey90", color = "black") +
    geom_tile(data = temp_df, aes(x = x, y = y, fill = value)) +
    geom_sf(data = tree_location, color = "grey78", size = 0.01, alpha = 0.9) +
    scale_fill_viridis_c(
      name = legend_temp_name,
      option = "inferno",
      breaks = range(temp_df$value, na.rm = TRUE),
      labels = round(range(temp_df$value, na.rm = TRUE), 1),
      guide = guide_colorbar(
        direction = "horizontal",
        label.position = "bottom"
      )) +
    scale_y_continuous(n.breaks = 5) +
    scale_x_continuous(
      n.breaks = 3,
      labels = function(x) {
        lab <- ifelse(
          x < 0,
          paste0(sprintf("%.2f", abs(x)), "°W"),
          paste0(sprintf("%.2f", x), "°E")
        )
        lab[seq_along(lab) > 2] <- ""
        lab
      }
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    coord_sf(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 18),
      legend.position = "bottom",
      legend.key.height = unit(0.3, "cm"),   
      legend.key.width  = unit(0.5, "cm"),   
      legend.title = element_text(size = 14, margin = margin(r = 10, b=10)), 
      legend.text  = element_text(size = 12),
      legend.spacing.y = unit(0.1, "cm"),
    )
  
  ## point to grid integration
  
  doy_rast <- {
    x <- terra::extract(
      clipped_over,
      vect(tree_location |> filter(!is.na(doy))),
      cells = TRUE
    ) |>
      bind_cols(
        st_drop_geometry(
          tree_location |> filter(!is.na(doy))
        )
      ) |>
      group_by(cell) |>
      summarise(doy = median(doy, na.rm = TRUE), .groups = "drop") |>
      filter(!is.na(cell))
    
    r <- rast(clipped_over)
    values(r) <- NA_real_
    values(r)[x$cell] <- x$doy
    r
  }
  phe_gridded_df <- as.data.frame(doy_rast, xy = TRUE, na.rm = TRUE) |>
    rename_with(~ "value", .cols = 3)
  
  map_phe_gridded <- ggplot() +
    geom_sf(data = city_boundary, fill = "grey90", color = "black") +
    geom_tile(data = phe_gridded_df, aes(x = x, y = y, fill = value)) +
    geom_sf(data = tree_location, color = "grey78", size = 0.01, alpha = 0.9) +
    scale_fill_viridis_c(name = legend_phe_name,
                         breaks = range(phe_gridded_df$value, na.rm = TRUE),
                         labels = round(range(phe_gridded_df$value, na.rm = TRUE), 0),
                         guide = guide_colorbar(
                           direction = "horizontal",
                           label.position = "bottom"
                         )) +
    scale_y_continuous(
      n.breaks = 5 
    ) +
    scale_x_continuous(
      n.breaks = 3,
      labels = function(x) {
        lab <- ifelse(
          x < 0,
          paste0(sprintf("%.2f", abs(x)), "°W"),
          paste0(sprintf("%.2f", x), "°E")
        )
        lab[seq_along(lab) > 2] <- ""
        lab
      }
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    coord_sf(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE,
             clip = "off"   # important
    ) +
    annotation_scale(
      location = "tr",
      width_hint = 0.4,
      text_cex = 1.3,
      unit_category = "metric",
      pad_y = unit(-0.7, "cm")
    ) +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),
        legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),   
        legend.key.width  = unit(0.5, "cm"),   
        legend.title = element_text(size = 14, margin = margin(r = 10, b=10)), 
        legend.text  = element_text(size = 12),
        legend.spacing.y = unit(0.1, "cm"),
      )
  
  ## change to model fit and point density
  
  phe_fitted <- read_rds(paste0(here("release/intermediate/spatial/inla/spatial_spde/"), city, "/canonicalname/", gsub(" ", "_", species), "_", phe_type, ".rds"))
  
  b0     <- phe_fitted$res$summary.fixed["Intercept", "mean"]
  b_temp <- phe_fitted$res$summary.fixed["temp",      "mean"]
  lwr_temp <- phe_fitted$res$summary.fixed["temp", "0.025quant"]
  upr_temp <- phe_fitted$res$summary.fixed["temp", "0.975quant"]
  b_prcp <- phe_fitted$res$summary.fixed["prcp",      "mean"]
  b_srad <- phe_fitted$res$summary.fixed["srad",      "mean"]
  u_year <- phe_fitted$res$summary.random$year[phe_fitted$res$summary.random$year$ID==year,"mean"]
  sd_b0     <- phe_fitted$res$summary.fixed["Intercept", "sd"]
  sd_b_temp <- phe_fitted$res$summary.fixed["temp",      "sd"]
  sd_b_prcp <- phe_fitted$res$summary.fixed["prcp",      "sd"]
  sd_b_srad <- phe_fitted$res$summary.fixed["srad",      "sd"]
  u_year_sd <- phe_fitted$res$summary.random$year[phe_fitted$res$summary.random$year$ID == year, "sd"]
  
  temp_col <- grep("^temp", names(city_input), value = TRUE)[1]
  prcp_col <- grep("^prcp", names(city_input), value = TRUE)[1]
  srad_col <- grep("^srad", names(city_input), value = TRUE)[1]
  
  mean_prcp <- mean(city_input[[prcp_col]], na.rm = TRUE)
  mean_srad <- mean(city_input[[srad_col]], na.rm = TRUE)
  
  temp_seq <- seq(
    min(city_input[[temp_col]], na.rm = TRUE),
    max(city_input[[temp_col]], na.rm = TRUE),
    length.out = 200
  )
  map_phe_gridded <- map_phe_gridded + theme(axis.title.y = element_blank(),
                                             axis.text.y = element_blank())
  map_phe <- map_phe + theme(axis.title.y = element_blank(),
                             axis.text.y = element_blank())
  
  if (b_temp > 0){text_color <- "#0D0887FF"} else {text_color <- "#FCA636FF"}
  
  interpret_text <- if (b_temp > 0) {
    if (phe_type == "SOS") {
      "Warmer microclimate, later spring."
    } else if (phe_type == "EOS") {
      "Warmer microclimate, later fall."
    } else {
      stop("Invalid phe_type. Use 'SOS' or 'EOS'.")
    }
  } else if (b_temp < 0) {
    if (phe_type == "SOS") {
      "Warmer microclimate, earlier spring."
    } else if (phe_type == "EOS") {
      "Warmer microclimate, earlier fall."
    } else {
      stop("Invalid phe_type. Use 'SOS' or 'EOS'.")
    }
  } else {
    paste0("In ", city, ", there is no association between ", season, " temperature and ", phe_type, ".")
  }
  city_text = city_label[city]
  
  map_panels_gridded <- map_temp + map_phe_gridded + 
    plot_annotation(
      title = bquote(.(city_text)* "," ~ .(year) ~~~ .(interpret_text)),
      subtitle = bquote(hat(beta)[T] == .(round(b_temp, 2)) ~
                          "[" * .(round(lwr_temp, 2)) * "," ~ .(round(upr_temp, 2)) * "]" ),
      theme = theme(
        plot.title = element_text(size = 22),
        plot.subtitle = element_text(color = text_color, size = 22)
      ))+
    plot_layout(guides = "collect")&
    theme(legend.position = "bottom")
  
  map_panels <- map_temp + map_phe + 
    plot_annotation(
      title = bquote(.(city_text)* "," ~ .(year) ~~~ .(interpret_text)),
      subtitle = bquote(hat(beta)[T] == .(round(b_temp, 2)) ~
                          "[" * .(round(lwr_temp, 2)) * "," ~ .(round(upr_temp, 2)) * "]" ),
      theme = theme(
        plot.title = element_text(size = 22),
        plot.subtitle = element_text(color = text_color, size = 22)
      ))+
    plot_layout(guides = "collect")&
    theme(legend.position = "bottom")
  
  return(list(
    map_panels = map_panels,
    map_panels_gridded = map_panels_gridded
  ))
}
