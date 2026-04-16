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
      filter(NAME == "Houston"),  # change to some distance
    "TP" = st_read(here("data/raw/Example_city/TP/Municipal_Boundary/Municipal_Boundary.shp")) %>%
      filter(NAME == "Tampa"),
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
    "Atlanta" = st_read(here("data/raw/Example_city/Atlanta/Official_Neighborhoods_-_Open_Data/Official_Neighborhoods_-_Open_Data.shp")) |>
      st_make_valid() |> st_union(),
    stop("Invalid city abbreviation.")
  )
  boundary <- st_transform(boundary, crs = 4326)  # ensure WGS84
  return(boundary)
}

plot_phe_temp_margin <- function(city = "NewYork", species = "Acer rubrum", phe_type = "SOS", year = 2023){
  source(here("script/spatial_stats_model/prepare_model_input.R"))
  model_input_sf <-  prepare_data(city = city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)|>
    filter(canonicalname == species) |>
    filter(year == !!year)
  
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
  tree_location <- tree_location|>
    st_transform(crs(city_boundary))
  city_v <- vect(city_boundary)
  clipped <- crop(clipped_over, city_v) |> mask(city_v)
  temp_df <- as.data.frame(clipped, xy = TRUE) |>
    rename_with(~ "value", .cols = 3)
  
  bbox_tree <- st_bbox(tree_location)
  bbox_city <- st_bbox(city_boundary)
  
  
  # the larger bbox
  xlim <- c(min(bbox_tree["xmin"], bbox_city["xmin"]), max(bbox_tree["xmax"], bbox_city["xmax"]))
  ylim <- c(min(bbox_tree["ymin"], bbox_city["ymin"]), max(bbox_tree["ymax"], bbox_city["ymax"]))
  if(city == "Sacramento"){
    xlim <- c(bbox_city["xmin"], bbox_city["xmax"])
    ylim <- c(bbox_city["ymin"], bbox_city["ymax"])
  }
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
      n.breaks = 3
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    coord_sf(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
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
    scale_x_continuous(n.breaks = 3) +
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
      n.breaks = 3
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
  
  ## change to model fit and point density
  
  phe_fitted <- read_rds(paste0("~/phenology-urban/data/proc/model_output/final_models_clean_outlier/inla/spatial_spde/", city, "/canonicalname/", gsub(" ", "_", species), "_", phe_type, ".rds"))
  
  b0     <- phe_fitted$res$summary.fixed["Intercept", "mean"]
  b_temp <- phe_fitted$res$summary.fixed["temp",      "mean"]
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
  
  pred_tree_location <- temp_df |>
    left_join(read_rds("~/phenology-urban/data/raw/CONUS/Daymet/example_pred_city/", city, "prcp_monthly_2022.rds"), by = c("x", "y")) |>
    left_join(read_rds("~/phenology-urban/data/raw/CONUS/Daymet/example_pred_city/", city, "srad_monthly_2022.rds"), by = c("x", "y")) |>
    mutate(
      doy_pred = b0 +
        b_temp * value +
        b_prcp * prcp +
        b_srad * srad +
        u_year
    )
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
      n.breaks = 3
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
    annotate(
      "text",
      x = Inf, y = Inf,     
      label = paste0("hat(beta)[Temp] == ", round(b_temp, 2)),
      hjust = 1.1, vjust = 1.1, 
      size  = 7,
      color = text_color,
      parse = TRUE
    )+
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
  
  if (city == "TP"){
    city <- "Tampa"
  } else if (city == "DT"){
    city <- "Detroit"
  }
  
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
  
  map_panels_gridded_model <- map_temp + map_phe_gridded_pred + 
    plot_annotation(
      title = bquote(.(city)* "," ~ .(year)),
      subtitle = bquote(hat(beta)[Temp] == .(round(b_temp, 2)) ~~~ .(interpret_text)),
      theme = theme(
        plot.title = element_text(size = 22),
        plot.subtitle = element_text(color = text_color, size = 22)
      ))+
    plot_layout(guides = "collect")&
    theme(legend.position = "bottom")
  
  map_panels_gridded <- map_temp + map_phe_gridded + 
    plot_annotation(
      title = bquote(.(city)* "," ~ .(year)),
      subtitle = bquote(hat(beta)[Temp] == .(round(b_temp, 2)) ~~~ .(interpret_text)),
      theme = theme(
        plot.title = element_text(size = 22),
        plot.subtitle = element_text(color = text_color, size = 22)
      ))+
    plot_layout(guides = "collect")&
    theme(legend.position = "bottom")
  
  map_panels <- map_temp + map_phe + 
    plot_annotation(
      title = bquote(.(city)* "," ~ .(year)),
      subtitle = bquote(hat(beta)[Temp] == .(round(b_temp, 2)) ~~~ .(interpret_text)),
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
