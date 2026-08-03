plot_city_coef_map <- function(df_coef) {
  
  boundary_us <- st_read(here("data/input/boundary/US/us_boundarydissolve.shp"), quiet = T) |>
    st_transform(4326)
  
  # df_coef: data frame with columns city_label, phe_type, Estimate, lon, lat
  model_results_loc <- df_coef |>
    group_by(city_label, phe_type) |>
    summarise(
      avg = mean(Estimate, na.rm = TRUE),
      sd  = sd(Estimate, na.rm = TRUE),
      lon = first(lon),
      lat = first(lat),
      .groups = "drop"
    )
  
  sos_range_vals <- range(model_results_loc |> filter(phe_type == "SOS") |> pull(avg), na.rm = TRUE)
  eos_range_vals <- range(model_results_loc |> filter(phe_type == "EOS") |> pull(avg), na.rm = TRUE)
  tau_vals <- unique(model_results_loc$sd) |> na.omit()
  tau_breaks <- pretty(range(tau_vals), n = 5)
  
  chelsa_ras_us <- terra::rast(here("data/input/climate/CHELSA/chelsa_us.tif"))
  bio1 <- chelsa_ras_us[[1]]
  bio12 <- chelsa_ras_us[[2]]
  
  brks <- c(14)
  p_sos = ggplot() +
    geom_sf(data = boundary_us , fill = "white", color = "transparent") +
    tidyterra::geom_spatraster(data = bio1, alpha = 0.4) +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(9, "RdBu")),
      name = "bio1 (°C)",
      na.value = "transparent"
    ) +
    # tidyterra::geom_spatraster_contour(data = terra::aggregate(bio1,fact = 100,fun = mean,na.rm = TRUE), 
    #                                    colour = "#EA7A71", linewidth = 3,alpha = 0.8, breaks = brks) +
    # tidyterra::geom_spatraster_contour_text(data = terra::aggregate(bio1,fact = 100,fun = mean,na.rm = TRUE), 
    #                                         breaks = brks, size =4, stroke = 0, colour = "#C0200E") +
    geom_sf(data = boundary_us , fill = "transparent", color = "grey85") +
    geom_point(data = model_results_loc |> 
                 filter(phe_type == "SOS"),
               aes(x = lon, y = lat,
                   size   = (sd)),
               shape  = 21,
               stroke = 1, # acts as the thickness of the boundary line
               colour = "grey20" # acts as the color of the boundary line
    ) +
    geom_point(
      data = model_results_loc |> 
        filter(phe_type == "SOS"),
      aes(x = lon, y = lat,  color = avg, size = sd), alpha =0.7
    ) +
    scale_color_gradientn(
      name = "Mean estimated\neffect",
      colors = c("#FCA636FF", "white", "#0D0887FF"),
      values = scales::rescale(c(round(sos_range_vals[1], 1), 0, round(sos_range_vals[2], 1))),
      breaks = c(round(sos_range_vals[1], 1), 0, round(sos_range_vals[2], 1)),
      limits = c(round(sos_range_vals[1], 1), round(sos_range_vals[2], 1)),
    )+
    facet_wrap(~ phe_type, ncol = 1) +
    scale_size_continuous(
      name = "Among-species variance",
      breaks = (tau_breaks),
      labels = (tau_breaks)^2,
      range = c(3, 12),
      limits = range(tau_vals)
    ) +
    guides(
      color = guide_colorbar(order = 1),
      size  = guide_legend(order = 2)
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      strip.text = element_text(size = 22, face = "bold"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) +
    labs(x = "Longitude", y = "Latitude") +
    coord_sf(expand = FALSE)
  
  brks <- c(300, 600, 1000, 1400, 2500)
  p_eos = ggplot() +
    geom_sf(data = boundary_us , fill = "white", color = "transparent") +
    tidyterra::geom_spatraster(data = bio12, alpha = 0.4) +
    scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(9, "YlGnBu"),
      name = "bio12 (mm)",
      na.value = "transparent"
    ) +
    # tidyterra::geom_spatraster_contour(data = terra::aggregate(bio1,fact = 100,fun = mean,na.rm = TRUE), 
    #                                    colour = "#EA7A71", linewidth = 3,alpha = 0.8, breaks = brks) +
    # tidyterra::geom_spatraster_contour_text(data = terra::aggregate(bio1,fact = 100,fun = mean,na.rm = TRUE), 
    #                                         breaks = brks, size =4, stroke = 0, colour = "#C0200E") +
    geom_sf(data = boundary_us , fill = "transparent", color = "grey85") +
    geom_point(data = model_results_loc |> 
                 filter(phe_type == "EOS"),
               aes(x = lon, y = lat, 
                   size   = (sd)),
               shape  = 21,
               stroke = 1, # acts as the thickness of the boundary line
               colour = "grey20" # acts as the color of the boundary line
    ) +
    geom_point(
      data = model_results_loc |> 
        filter(phe_type == "EOS"),
      aes(x = lon, y = lat, color = avg, size = sd), alpha =0.7
    ) +
    scale_color_gradientn(
      name = "Mean estimated\neffect",
      colors = c("#FCA636FF", "white", "#0D0887FF"),
      values = scales::rescale(c(round(eos_range_vals[1], 1), 0, round(eos_range_vals[2], 1))),
      breaks = c(round(eos_range_vals[1], 1), 0, round(eos_range_vals[2], 1)),
      limits = c(round(eos_range_vals[1], 1), round(eos_range_vals[2], 1)),
    )+
    facet_wrap(~ phe_type, ncol = 1) +
    scale_size_continuous(
      name = "Among-species variance",
      breaks = (tau_breaks),
      labels = (tau_breaks)^2,
      range = c(3, 12),
      limits = range(tau_vals)
    ) +
    guides(
      color = guide_colorbar(order = 1),
      size  = guide_legend(order = 2)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 16),
      axis.title.x = element_text(size = 20),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 22, face = "bold"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) +
    labs(x = "Longitude") +
    coord_sf(expand = FALSE)
  
  p_sos_fill <- p_sos +
    guides(
      color = guide_colorbar(
        direction = "horizontal",
        order = 1,
        label.theme = element_text(angle = 45, hjust = 1, vjust = 1)),
      size  = "none",
      fill  = guide_colorbar(
        direction = "horizontal",
        order = 1,
        label.theme = element_text(angle = 45, hjust = 1, vjust = 1))
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "left",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10)
    )
  
  p_eos_fill <- p_eos +
    guides(
      color = guide_colorbar(
        direction = "horizontal",
        order = 1,
        label.theme = element_text(angle = 45, hjust = 1, vjust = 1)),
      size  = "none",
      fill  = guide_colorbar(
        direction = "horizontal",
        order = 1,
        label.theme = element_text(angle = 45, hjust = 1, vjust = 1))
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "left",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10)
    )
  
  left <- p_sos_fill + p_eos_fill + 
    plot_layout(ncol = 2)
  
  return(left)
}