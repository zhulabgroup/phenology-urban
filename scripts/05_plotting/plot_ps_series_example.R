if (!require("BatchPlanet")){
  install.packages('BatchPlanet', repos = c('https://yiluansong.r-universe.dev', 'https://cloud.r-project.org'))
}
pacman::p_load(BatchPlanet)

plot_ps_series_example <- function() {
  doy <- read_rds(here("data/input/phenology/doy/doy_NewYork_Acerrubrum.rds"))|>
    filter(id == 7639) |>
    filter(year == 2023)
  clean_ts <- read_rds(here("data/input/phenology/clean/clean_NewYork_Acerrubrum.rds")) |>
    filter(id == 7639) |>
    filter(date >= as.Date("2022-10-01") & date <= as.Date("2024-3-31"))

  df_ts <- clean_ts %>% dplyr::mutate(value_smooth = whittaker_smoothing_filling(evi, lambda = 50))
  df_doy <- doy %>%
    dplyr::mutate(date_doy = as.Date(stringr::str_c(year, "-01-01")) + doy - 1,
                  date_start = as.Date(stringr::str_c(year, "-01-01")) + start - 1,
                  date_end = as.Date(stringr::str_c(year, "-01-01")) + end - 1)
  y_min <- min(df_ts$evi, na.rm = TRUE)
  y_max <- max(df_ts$evi, na.rm = TRUE)
  
  ggplot() +
    geom_point(data = df_ts, aes(x = date, y = evi), color = "grey", size = 2, alpha = 0.8) +
    labs(x = "Time", y = "Enhanced Vegetation Index (EVI)") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    ) +
    geom_line(data = df_ts, aes(x = date, y = value_smooth), color = "grey40", linewidth = 1.2)+
    geom_segment(data = df_doy |> filter(direction == "up", thres == 0.5),
                 aes(x = date_doy, xend = date_doy, y = y_min, yend = y_max),
                 linewidth = 1.5, col = "dark green", inherit.aes = F) +
    geom_rect(data = df_doy |> filter(direction == "up", thres == 0.5),
              aes(xmin = date_start, xmax = date_end, ymin = y_min, ymax = y_max),
              alpha = 0.1,fill = "dark green", inherit.aes = F) +
    geom_segment(data = df_doy |> filter(direction == "down", thres == 0.5),
                 aes(x = date_doy, xend = date_doy, y = y_min, yend = y_max),
                 linewidth = 1.5, col = "orange", inherit.aes = F) +
    geom_rect(data = df_doy |> filter(direction == "down", thres == 0.5),
              aes(xmin = date_start, xmax = date_end, ymin = y_min, ymax = y_max),
              alpha = 0.1,fill = "orange", inherit.aes = F) +
    geom_text( data = df_doy |> filter(direction == "up", thres == 0.5),
               aes(x = date_doy + 1, y = y_max, label = "SOS = 110"),color = "darkgreen",
               vjust = -0.5,size = 5, inherit.aes = FALSE
    ) +
    geom_text(data = df_doy |> filter(direction == "down", thres == 0.5),
              aes(x = date_doy + 1, y = y_max, label = "EOS = 298"),
              color = "orange",vjust = -0.5,size = 5,inherit.aes = FALSE
    )
}