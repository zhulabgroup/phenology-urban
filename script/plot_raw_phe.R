plot_species_raw_phe <- function(species, all_city_input, color_palette, phe_type) {
  species_city <- all_city_input |>
    filter(canonicalname == species)|>
    dplyr::group_by(city) |>
    dplyr::summarise(
      median_doy = median(doy, na.rm = TRUE),
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
      IQR_d = quantile(median_doy, 0.75, na.rm = TRUE) - quantile(median_doy, 0.25, na.rm = TRUE),
      lower_IQR_d = max(quantile(median_doy, 0.25, na.rm = TRUE), min(median_doy, na.rm = TRUE)),
      upper_IQR_d = min(quantile(median_doy, 0.75, na.rm = TRUE), max(median_doy, na.rm = TRUE)),
      .groups = "drop")
  
  sp_boxplot_SOS <- ggplot() +
    geom_rect(
      data = integ_city,
      aes(xmin = 0.5, xmax = length((all_city_input |>
                                       filter(canonicalname == species) |>
                                       distinct(city) |>
                                       pull())),
          ymin = lower_IQR_d, ymax = upper_IQR_d),
      fill = "grey80", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_abline(slope = 0, intercept = integ_city$median_doy_cities, linetype = "dotted", color = "grey70", size = 1) +
    geom_errorbar(data = species_city,
                  aes(x = city, ymin = lower_IQR_d, ymax = upper_IQR_d, color = city),
                  width = 0.2, linewidth = 1, alpha = 0.5, inherit.aes = FALSE) +
    geom_point(data = species_city,
               aes(x = city, y = median_doy, color = city),
               size = 2, alpha = 0.5, inherit.aes = FALSE) +
    scale_color_manual(values = color_palette) +
    coord_flip(ylim = case_when(
          phe_type == "SOS" ~ c(0, 150),
          phe_type == "EOS" ~ c(200, 350)
        )) +
    labs(
      x = "City",
      title = paste0("<i>", species, "</i>")
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
      title = element_text(size = 16),
      legend.position = "none"
    )
  
  return(sp_boxplot_SOS)
}

# Archive

# sp_density_SOS <- ggplot() +
#   geom_density(data = all_city_input|>
#                  filter(canonicalname == species),
#                aes(x = doy, fill = city),
#                color = alpha("black", 0.3),
#                linewidth= 0.3, alpha = 0.5) +
#   scale_fill_manual(values = color_palette) +
#   scale_y_continuous(
#     name = "Density"
#   ) +
#   labs(x = case_when(
#     phe_type == "SOS" ~ "Start-of-season (DOY)",
#     phe_type == "EOS" ~ "End-of-season (DOY)"
#   )) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(size = 10),
#     axis.text.y = element_text(size = 10),
#     title = element_text(size = 14),
#     legend.position = "none"
#   ) +
#   ggtitle(paste0( species))+
#   scale_x_continuous(limits = case_when(
#     phe_type == "SOS" ~ c(-50, 250),
#     phe_type == "EOS" ~ c(150, 450)
#   ))+
#   theme(axis.title.x = element_blank(),
#         axis.text.x  = element_blank(),
#         axis.ticks.x = element_blank())

# sp_boxplot_SOS <- ggplot() +
#   geom_errorbar(data = all_city_input |>
#                   filter(canonicalname == species)|>
#                   dplyr::group_by(city) |>
#                   dplyr::summarise(
#                     median_doy = median(doy, na.rm = TRUE),
#                     IQR_d = quantile(doy, 0.75, na.rm = TRUE) - quantile(doy, 0.25, na.rm = TRUE),
#                     lower_IQR_d = max(quantile(doy, 0.05, na.rm = TRUE), min(doy, na.rm = TRUE)),
#                     upper_IQR_d = min(quantile(doy, 0.95, na.rm = TRUE), max(doy, na.rm = TRUE)),
#                     median_temp = median(doy, na.rm = TRUE),
#                     .groups = "drop"),
#                 aes(x = city, ymin = lower_IQR_d, ymax = upper_IQR_d, color = city),
#                 width = 0.2, linewidth = 1, alpha = 0.5) +
#   geom_point(data = all_city_input |>
#                filter(canonicalname == species)|>
#                dplyr::group_by(city) |>
#                dplyr::summarise(
#                  median_doy = median(doy, na.rm = TRUE),
#                  .groups = "drop"),
#              aes(x = city, y = median_doy, color = city),
#              size = 2, alpha = 0.5) +
#   scale_color_manual(values = color_palette) +
#   coord_flip(ylim = case_when(
#     phe_type == "SOS" ~ c(-50, 250),
#     phe_type == "EOS" ~ c(150, 450)
#   )) +
#   labs(x = "City") +
#   labs(y = case_when(
#     phe_type == "SOS" ~ "Start-of-season (DOY)",
#     phe_type == "EOS" ~ "End-of-season (DOY)"
#   )) +
#   theme_classic() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_text(size = 10),
#     title = element_text(size = 14),
#     legend.position = "none"
#   )


# sp_IQR_SOS <- ggplot() +
#   geom_point(data = all_city_input |>
#                   filter(canonicalname == species)|>
#                   dplyr::group_by(city) |>
#                   dplyr::summarise(
#                     median_doy = median(doy, na.rm = TRUE),
#                     IQR_d = quantile(doy, 0.75, na.rm = TRUE) - quantile(doy, 0.25, na.rm = TRUE),
#                     lower_IQR_d = max(quantile(doy, 0.25, na.rm = TRUE), min(doy, na.rm = TRUE)),
#                     upper_IQR_d = min(quantile(doy, 0.75, na.rm = TRUE), max(doy, na.rm = TRUE)),
#                     median_temp = median(doy, na.rm = TRUE),
#                     .groups = "drop"), 
#                aes(x = 1, y = IQR_d, color = city), 
#                size = 2, alpha = 0.5) +
#   scale_color_manual(values = color_palette) + 
#   coord_flip(ylim = case_when(
#     phe_type == "SOS" ~ c(-50, 250),
#     phe_type == "EOS" ~ c(150, 450)
#   )) +
#   labs(y = case_when(
#     phe_type == "SOS" ~ "IQR for Start-of-season (DOY)",
#     phe_type == "EOS" ~ "IQR for End-of-season (DOY)"
#   )) +
#   theme_classic() +
#   theme(
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_text(size = 10),
#     title = element_text(size = 14),
#     legend.position = "none"
#   )


# sum_density_SOS <- ggplot() +
#   geom_density(data = all_city_input |>
#                   filter(canonicalname == species)|>
#                   dplyr::group_by(city) |>
#                   dplyr::summarise(
#                     median_doy = median(doy, na.rm = TRUE),
#                     IQR_d = quantile(doy, 0.75, na.rm = TRUE) - quantile(doy, 0.25, na.rm = TRUE),
#                     lower_IQR_d = max(quantile(doy, 0.25, na.rm = TRUE), min(doy, na.rm = TRUE)),
#                     upper_IQR_d = min(quantile(doy, 0.75, na.rm = TRUE), max(doy, na.rm = TRUE)),
#                     .groups = "drop"),
#                color = alpha("black"),
#                aes(x = median_doy),
#                linewidth= 0.5, fill = "grey60", alpha = 0.5) +
#   scale_y_continuous(
#     name = "Density",
#     breaks = seq(0, 0.04, by = 0.02)
#   ) +
#   labs(x = case_when(
#     phe_type == "SOS" ~ "Median start-of-season by city (DOY)",
#     phe_type == "EOS" ~ "Median end-of-season by city (DOY)"
#   )) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(size = 10),
#     axis.text.y = element_text(size = 10),
#     title = element_text(size = 14),
#     legend.position = "none"
#   ) +
#   scale_x_continuous(limits = case_when(
#     phe_type == "SOS" ~ c(-50, 250),
#     phe_type == "EOS" ~ c(150, 450)
#   )) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x  = element_blank(),
#         axis.ticks.x = element_blank())
# 
# sum_boxplot_SOS <- ggplot() +
#   geom_errorbar(data = all_city_input |>
#                   filter(canonicalname == species)|>
#                   dplyr::group_by(city) |>
#                   dplyr::summarise(
#                     median_doy = median(doy, na.rm = TRUE),
#                     .groups = "drop") |>
#                   summarise(
#                     median_doy_cities = median(median_doy, na.rm = TRUE),
#                     IQR_d = quantile(median_doy, 0.75, na.rm = TRUE) - quantile(median_doy, 0.25, na.rm = TRUE),
#                     lower_IQR_d = max(quantile(median_doy, 0.05, na.rm = TRUE), min(median_doy, na.rm = TRUE)),
#                     upper_IQR_d = min(quantile(median_doy, 0.95, na.rm = TRUE), max(median_doy, na.rm = TRUE)),
#                     .groups = "drop"),
#                 aes(x = "none", ymin = lower_IQR_d, ymax = upper_IQR_d),
#                 color = "grey60", width = 0.2, linewidth = 1, alpha = 0.5) +
#   geom_point(data = all_city_input |>
#                filter(canonicalname == species)|>
#                dplyr::group_by(city) |>
#                dplyr::summarise(
#                  median_doy = median(doy, na.rm = TRUE),
#                  .groups = "drop") |>
#                summarise(
#                  median_doy_cities = median(median_doy, na.rm = TRUE),
#                  IQR_d = quantile(median_doy, 0.75, na.rm = TRUE) - quantile(median_doy, 0.25, na.rm = TRUE),
#                  lower_IQR_d = max(quantile(median_doy, 0.25, na.rm = TRUE), min(median_doy, na.rm = TRUE)),
#                  upper_IQR_d = min(quantile(median_doy, 0.75, na.rm = TRUE), max(median_doy, na.rm = TRUE)),
#                  .groups = "drop"),
#              aes(x = "none", y = median_doy_cities),
#              color = "grey60", size = 2, alpha = 0.5) +
#   # geom_rect(
#   #   data = all_city_input |>
#   #     filter(canonicalname == species) %>%
#   #     summarise(
#   #       lower_IQR = quantile(doy, 0.25, na.rm = TRUE),
#   #       upper_IQR = quantile(doy, 0.75, na.rm = TRUE),
#   #       .groups = "drop"
#   #     ),
#   #   aes(xmin = 0.5, xmax = length((all_city_input |>
#   #                                    filter(canonicalname == species) |>
#   #                                    distinct(city) |>
#   #                                    pull())),
#   #       ymin = lower_IQR, ymax = upper_IQR),
#   #   fill = "grey60", alpha = 0.5, inherit.aes = FALSE
#   # ) +
#   coord_flip(ylim = case_when(
#     phe_type == "SOS" ~ c(-50, 250),
#     phe_type == "EOS" ~ c(150, 450)
#   )) +
#   labs(y = case_when(
#     phe_type == "SOS" ~ "Median start-of-season by city (DOY)",
#     phe_type == "EOS" ~ "Median end-of-season by city (DOY)"
#   )) +
#   theme_classic() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.x = element_text(size = 10),
#     title = element_text(size = 14),
#     legend.position = "none"
#   )
# 

# sum_IQR_SOS <- ggplot() +
#   geom_point(data = all_city_input |>
#                   filter(canonicalname == species)|>
#                   dplyr::group_by(city) |>
#                   dplyr::summarise(
#                     median_doy = median(doy, na.rm = TRUE),
#                     .groups = "drop") |>
#                   summarise(
#                     median_doy_cities = median(median_doy, na.rm = TRUE),
#                     IQR_d = quantile(median_doy, 0.75, na.rm = TRUE) - quantile(median_doy, 0.25, na.rm = TRUE),
#                     lower_IQR_d = max(quantile(median_doy, 0.25, na.rm = TRUE), min(median_doy, na.rm = TRUE)),
#                     upper_IQR_d = min(quantile(median_doy, 0.75, na.rm = TRUE), max(median_doy, na.rm = TRUE)),
#                     .groups = "drop"), 
#                aes(x = "none", y = IQR_d), 
#                color = "grey60", size = 2, alpha = 0.5) +
#   coord_flip(ylim = case_when(
#     phe_type == "SOS" ~ c(-50, 250),
#     phe_type == "EOS" ~ c(150, 450)
#   )) +
#   labs(y = case_when(
#     phe_type == "SOS" ~ "IQR for median start-of-season by city (DOY)",
#     phe_type == "EOS" ~ "IQR for median end-of-season by city (DOY)"
#   )) +
#   theme_classic() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.x = element_text(size = 10),
#     title = element_text(size = 14),
#     legend.position = "none"
#   )

# top <- wrap_plots(
#   sp_histgram_SOS, sp_boxplot_SOS,
#   ncol = 1, heights = c(3, 3)
# )
# 
# bottom <- wrap_plots(
#   sum_histgram_SOS, sum_boxplot_SOS,
#   ncol = 1, heights = c(1, 1)
# )
# 
# final_plot <- wrap_plots(
#   top, bottom,
#   ncol = 1, heights = c(6, 2)
# )
# 

