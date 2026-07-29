plot_phe_temp_margin <- function(city = "NewYork", species = "Acer rubrum", phe_type = "SOS", year = 2022){
  source(here("script/spatial_stats_model/prepare_model_input.R"))
  model_input_sf <-  prepare_data(city = city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)|>
    filter(canonicalname == species) |>
    filter(year == !!year)
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
  file <- paste0("~/lab-data/datasets/climate/UHAT/clean/",city, "_", season, "_", year, ".tif")
  temp_tif <- terra::rast(file)
  tifs <- list.files(
    "~/lab-data/datasets/climate/UHAT/grid",
    pattern = "^MYD11A1\\..*LST_Day_1km\\.tif$",
    full.names = TRUE
  )
  r_template <- terra::rast(tifs[1])
  crs(temp_tif) <- crs(r_template)  # ensure correct CRS
  tree_location <- tree_location|>
    st_transform(crs(temp_tif))
  
  clipped_over <- temp_tif |>
    # crop(temp_tif, vect(tree_location)) |>
    project(crs(city_boundary), res = 0.009)  # transform to WGS84 for plotting
  city_v <- vect(city_boundary)
  clipped <- crop(clipped_over, city_v) |> mask(city_v)
  temp_df <- as.data.frame(clipped, xy = TRUE) |>
    rename_with(~ "value", .cols = 3)
  
  tree_location <- tree_location|>
    st_transform(crs(city_boundary))  
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
  
  phe_fitted <- read_rds(paste0("~/phenology-urban/data/proc/model_output/final_models_clean_outlier/inla/spatial_spde/", city, "/canonicalname/", gsub(" ", "_", species), "_", phe_type, ".rds"))
  
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
  
  ####################### Prediction map #######################
  
  pred_tree_location <- read_rds(paste0("~/phenology-urban/data/raw/CONUS/Daymet/example_pred_city/", city, "/pred_monthly_2022.rds")) |>
    mutate(
      doy_pred = b0 +
        b_temp * value +
        b_prcp * prcp +
        b_srad * srad +
        u_year
    ) |>
    st_as_sf(coords = c("x", "y"), crs = 4326)
  doy_rast_pred <- {
    x <- terra::extract(
      clipped,
      vect(pred_tree_location |> filter(!is.na(doy_pred))),
      cells = TRUE
    ) |>
      bind_cols(
        st_drop_geometry(
          pred_tree_location |> filter(!is.na(doy_pred))
        )
      ) |>
      group_by(cell) |>
      summarise(doy_pred = median(doy_pred, na.rm = TRUE), .groups = "drop") |>
      filter(!is.na(cell))
    
    r <- rast(clipped)
    values(r) <- NA_real_
    values(r)[x$cell] <- x$doy_pred
    r
  }
  pred_phe_gridded_df <- as.data.frame(doy_rast_pred, xy = TRUE, na.rm = TRUE) |>
    rename_with(~ "value", .cols = 3)
  
  map_phe_gridded_pred <- ggplot() +
    geom_sf(data = city_boundary, fill = "grey90", color = "black") +
    geom_tile(data = pred_phe_gridded_df, aes(x = x, y = y, fill = value)) +
    geom_sf(data = tree_location, color = "grey78", size = 0.01, alpha = 0.9) +
    scale_fill_viridis_c(name = paste0("Modelled ", legend_phe_name),
                         breaks = range(pred_phe_gridded_df$value, na.rm = TRUE),
                         labels = round(range(pred_phe_gridded_df$value, na.rm = TRUE), 0),
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
    coord_sf(
      xlim = lims$xlim, ylim = lims$ylim,
      expand = FALSE,
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
  
  
  # 5. Predicted y for year
  y_mean_year <- b0 + b_temp * temp_seq +
    b_prcp * mean_prcp +
    b_srad * mean_srad +
    u_year   # <-- add random effect
  eta_var <- sd_b0^2 +
    (temp_seq^2)       * sd_b_temp^2 +
    (mean_prcp^2)      * sd_b_prcp^2 +
    (mean_srad^2)      * sd_b_srad^2 +
    u_year_sd^2
  
  eta_sd <- sqrt(eta_var)
  
  # 95% CI
  z_975 <- qnorm(0.975)
  y_lower_year <- y_mean_year - z_975 * eta_sd
  y_upper_year <- y_mean_year + z_975 * eta_sd
  
  plot_data <- data.frame(
    temp   = temp_seq,
    y_mean = y_mean_year,
    y_low  = y_lower_year,
    y_high = y_upper_year
  )
  
  if (b_temp > 0){
    text_color <- "#0D0887FF"
  } else {
    text_color <- "#FCA636FF"}
  
  sp_city_scatter_plot <- ggplot() +
    geom_hex(data = city_input, aes(x = .data[[temp_col]], y = doy), bins = floor(nrow(city_input)/90 + 6), alpha = 0.5, color = "NA") +
    # geom_point(data = city_input, aes(x = .data[[temp_col]], y = doy),alpha = 0.5) +
    # geom_smooth(data = city_input, aes(x = .data[[temp_col]], y = doy), method = "lm", se = FALSE, color = "blue", linewidth = 1) +
    scale_fill_gradient(
      low = "grey70", 
      high = "grey10" 
    ) +
    # facet_wrap(~ year) +
    # geom_ribbon(data = plot_data, aes(x = temp, y = y_mean, ymin = y_low, ymax = y_high), fill = "#FCA636FF", alpha = 0.2) +
    geom_line(data = plot_data, aes(x = temp, y = y_mean), color = text_color, linewidth = 1) +
    theme_bw(base_size = 14)+
    labs(x = paste0("Mean ", season, " temperature (°C)"), 
         y = case_when(
           phe_type == "SOS" ~ "Start-of-season (DOY)",
           phe_type == "EOS" ~ "End-of-season (DOY)"
         )) +
    theme_classic() +
    # annotate(
    #   "text",
    #   x = Inf, y = Inf,     
    #   label = paste0("hat(beta)[T] == ", round(b_temp, 2)),
    #   hjust = 1.1, vjust = 1.1, 
    #   size  = 7,
    #   color = text_color,
    #   parse = TRUE
    # )+
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      title = element_text(size = 14),
      legend.position = "none"
    )
  map_phe_gridded_pred <- map_phe_gridded_pred + theme(axis.title.y = element_blank(),
                                                       axis.text.y = element_blank())
  map_phe_gridded <- map_phe_gridded + theme(axis.title.y = element_blank(),
                                             axis.text.y = element_blank())
  map_phe <- map_phe + theme(axis.title.y = element_blank(),
                             axis.text.y = element_blank())
  
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
  map_panels_gridded_model <- map_temp + map_phe_gridded_pred + 
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
    sp_city_scatter_plot = sp_city_scatter_plot,
    map_panels = map_panels,
    map_panels_gridded = map_panels_gridded,
    map_panels_gridded_model = map_panels_gridded_model
  ))
}
