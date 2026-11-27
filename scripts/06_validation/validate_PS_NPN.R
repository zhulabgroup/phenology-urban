plot_validation_ps_npn <- function(
    npn = readr::read_csv(
      here("data/input/NPN/NPN1624/individual_phenometrics_data.csv"),
      col_types = readr::cols(Latitude = readr::col_double(), Longitude = readr::col_double())
    ),
    urban_areas = sf::st_read(
      here("data/input/boundary/USCities/tl_2024_us_uac20.shp"),
      quiet = TRUE
    ),
    npn_points = readRDS(
      here("data/input/tree_list/summary/matching_npn.rds")
    ),
    ps_doy_data = dplyr::bind_rows(lapply(
      list.files(here("data/input/NPN/PS_valid"),
                 pattern = "[.]rds$", full.names = TRUE),
      readRDS
    )),
    threshold = 21,
    color_palette = NULL) {
  
  urban_areas <- sf::st_transform(urban_areas, 4326)
  npn <- npn |>
    filter(!Species_Functional_Type %in%
             c("Algae", "Evergreen forb", "Forb", "Graminoid", "Semi-evergreen forb"),
           NumDays_Since_Prior_No > 0, NumDays_Until_Next_No > 0,
           NumDays_Since_Prior_No <= threshold, NumDays_Until_Next_No <= threshold) |>
    sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) |>
    sf::st_filter(urban_areas, .predicate = sf::st_within) |>
    sf::st_drop_geometry() |>
    mutate(
      Species = paste0(Genus, " ", Species),
      Phenophase_Description = ifelse(
        Phenophase_Description == "Breaking leaf buds (lilac/honeysuckle)",
        "Breaking leaf buds", Phenophase_Description
      )
    ) |>
    filter(Phenophase_Description != "Young leaves (forbs)",
           Observed_Status_Conflict_Flag == "-9999")
  
  npn_filter <- npn |>
    group_by(Individual_ID, Phenophase_Description) |>
    arrange(First_Yes_Julian_Date) |>
    mutate(
      gap = First_Yes_Julian_Date - lag(Last_Yes_Julian_Date,
                                        default = first(First_Yes_Julian_Date)),
      new_series = ifelse(is.na(gap) | gap > threshold, 1, 0),
      series_id = cumsum(new_series)
    ) |>
    group_by(Individual_ID, Phenophase_Description, series_id) |>
    summarise(
      First_Yes_Julian_Date = min(First_Yes_Julian_Date),
      Last_Yes_Julian_Date = max(Last_Yes_Julian_Date),
      First_Yes_Year = min(First_Yes_Year),
      Last_Yes_Year = max(Last_Yes_Year),
      NumYs_in_Series = sum(NumYs_in_Series, na.rm = TRUE),
      .groups = "drop"
    )
  
  individual_years <- npn_filter |>
    filter(NumYs_in_Series > 1) |>
    group_by(Individual_ID, Phenophase_Description) |>
    summarise(
      years_covered = list(unique(unlist(purrr::map2(First_Yes_Year, Last_Yes_Year, seq)))),
      .groups = "drop"
    ) |>
    mutate(Observe_Year = purrr::map_int(years_covered, length)) |>
    filter(Observe_Year > 1) |>
    select(Individual_ID, Phenophase_Description, Observe_Year) |>
    tidyr::pivot_wider(names_from = Phenophase_Description, values_from = Observe_Year) |>
    mutate(cover_phenophase = rowSums(!is.na(across(-Individual_ID)))) |>
    filter(!(cover_phenophase == 1 & !is.na(Leaves)))
  
  compute_npnphe <- function(start, end) {
    if (!all(c(start, end) %in% names(individual_years))) {
      stop("Missing required NPN phenophases: ",
           paste(setdiff(c(start, end), names(individual_years)), collapse = ", "))
    }
    npn_filter |>
      filter(
        Individual_ID %in% (individual_years |>
                              filter(.data[[start]] >= 1, .data[[end]] >= 1) |>
                              pull(Individual_ID)),
        Phenophase_Description %in% c(start, end)
      ) |>
      mutate(year = ifelse(Phenophase_Description == start, First_Yes_Year, Last_Yes_Year)) |>
      group_by(Individual_ID, year, Phenophase_Description) |>
      summarise(
        doy = if (first(Phenophase_Description) == start) {
          min(First_Yes_Julian_Date, na.rm = TRUE)
        } else {
          max(Last_Yes_Julian_Date, na.rm = TRUE)
        },
        .groups = "drop"
      ) |>
      summarise(npn_doy = mean(doy, na.rm = TRUE), n = n(),
                .by = c(Individual_ID, year)) |>
      filter(n == 2) |>
      select(-n) |>
      mutate(npn_doy = lubridate::yday(lubridate::as_date(npn_doy - 2440588,
                                                          origin = "1970-01-01"))) |>
      left_join(npn |> distinct(Individual_ID, Species, Species_Functional_Type,
                                Longitude, Latitude), by = "Individual_ID")
  }
  
  match_ps <- function(npn_phe, direction, min_doy, max_doy) {
    npn_phe |>
      left_join(npn_points, by = c("Longitude" = "lon", "Latitude" = "lat")) |>
      left_join(ps_doy_data |> filter(.data$direction == .env$direction),
                by = c("new_id" = "id", "year")) |>
      filter(npn_doy >= min_doy, npn_doy <= max_doy,
             doy >= min_doy, doy <= max_doy) |>
      tidyr::drop_na()
  }
  
  npn_sos_ps <- match_ps(compute_npnphe("Leaves", "Increasing leaf size"),
                         "up", 0, 150)
  npn_eos_ps <- match_ps(compute_npnphe("Colored leaves", "Falling leaves"),
                         "down", 150, 365)
  
  summarize_phe <- function(data, by_species = FALSE) {
    if (by_species) {
      data <- data |>
        summarise(npn_doy = mean(npn_doy), doy = mean(doy),
                  .by = c(Species, city, year, Longitude, Latitude))
      grouping <- c("Species", "city", "year")
    } else {
      grouping <- c("city", "year", "Longitude", "Latitude")
    }
    data |>
      summarise(
        n = n(), npn_doy_mean = mean(npn_doy, na.rm = TRUE),
        npn_doy_sd = sd(npn_doy, na.rm = TRUE),
        ps_doy_mean = mean(doy, na.rm = TRUE),
        ps_doy_sd = sd(doy, na.rm = TRUE),
        .by = all_of(grouping)
      ) |>
      mutate(
        npn_low_bound = npn_doy_mean - npn_doy_sd,
        npn_up_bound = npn_doy_mean + npn_doy_sd,
        ps_low_bound = ps_doy_mean - ps_doy_sd,
        ps_up_bound = ps_doy_mean + ps_doy_sd
      )
  }
  
  loc_sos <- summarize_phe(npn_sos_ps)
  loc_eos <- summarize_phe(npn_eos_ps)
  sp_sos <- summarize_phe(npn_sos_ps, by_species = TRUE)
  sp_eos <- summarize_phe(npn_eos_ps, by_species = TRUE)
  all_cities <- sort(unique(c(loc_sos$city, loc_eos$city)))
  all_cities <- all_cities[!is.na(all_cities)]
  if (!length(all_cities)) stop("No matching NPN and PlanetScope observations found")
  if (is.null(color_palette)) {
    color_palette <- stats::setNames(
      withr::with_seed(20260320, randomcoloR::distinctColorPalette(length(all_cities))),
      all_cities
    )
  }
  if (!all(all_cities %in% names(color_palette))) {
    stop("color_palette must be named and contain every city")
  }
  
  plot_panel <- function(data, phe) {
    data <- data |> filter(is.finite(npn_doy_mean), is.finite(ps_doy_mean))
    if (!nrow(data)) stop("No matching observations to plot for ", phe)
    x <- data$npn_doy_mean
    y <- data$ps_doy_mean
    correlation <- if (length(x) >= 3 && sd(x) > 0 && sd(y) > 0) {
      stats::cor.test(x, y, method = "pearson")
    } else NULL
    label <- sprintf(
      " r = %s\n p = %s\nRMSE = %.2f",
      if (is.null(correlation)) "NA" else sprintf("%.3f", unname(correlation$estimate)),
      if (is.null(correlation)) "NA" else format.pval(correlation$p.value, digits = 3),
      sqrt(mean((y - x)^2))
    )
    x_label <- if (phe == "SOS") {
      "Midpoint of first leaf and\nlast increasing leaf size (NPN)"
    } else {
      "Midpoint of first colored leaves and\nlast falling leaves (NPN)"
    }
    ggplot(data, aes(x = npn_doy_mean, y = ps_doy_mean)) +
      geom_segment(aes(x = npn_low_bound, xend = npn_up_bound,
                       yend = ps_doy_mean, color = city),
                   linewidth = 0.5, alpha = 0.3, na.rm = TRUE) +
      geom_segment(aes(y = ps_low_bound, yend = ps_up_bound,
                       xend = npn_doy_mean, color = city),
                   linewidth = 0.5, alpha = 0.3, na.rm = TRUE) +
      geom_point(aes(color = city, size = n), alpha = 0.4) +
      geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
      annotate("text", x = -Inf, y = Inf, label = label,
               hjust = -0.1, vjust = 1.2, size = 4) +
      scale_color_manual(values = color_palette, limits = all_cities, drop = TRUE) +
      scale_size_continuous(name = "Sample size", range = c(1, 5),
                            breaks = function(x) round(seq(x[1], x[2], length.out = 3))) +
      guides(color = "none", size = guide_legend(override.aes = list(color = "grey40", alpha = 0.8))) +
      labs(x = x_label, y = paste0("Day of 50% ",
                                   if (phe == "SOS") "green-up" else "green-down",
                                   " (PlanetScope)")) +
      theme_minimal(base_size = 14) +
      theme(
        axis.title = element_text(size = 16),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(0.80, 0.2),
        legend.background = element_rect(fill = NA, color = NA),
        legend.key = element_blank(), strip.background = element_blank()
      )
  }
  
  legend_data <- data.frame(city = factor(all_cities, levels = all_cities),
                            x = 1, y = seq_along(all_cities))
  legend_plot <- ggplot(legend_data, aes(x, y, color = city)) +
    geom_linerange(aes(xmin = x - 1, xmax = x + 1), linewidth = 0.5, alpha = 0.3) +
    geom_point(size = 1, alpha = 0.3) +
    scale_color_manual(values = color_palette[all_cities], breaks = all_cities,
                       limits = all_cities, name = "City") +
    guides(color = guide_legend(ncol = 2)) +
    theme_void() +
    theme(legend.position = "right", legend.text = element_text(size = 10),
          legend.title = element_text(size = 12))
  shared_legend <- cowplot::get_legend(legend_plot)
  main_plot <- plot_panel(loc_sos, "SOS") + plot_panel(loc_eos, "EOS") +
    plot_panel(sp_sos, "SOS") + plot_panel(sp_eos, "EOS") +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))
  cowplot::plot_grid(main_plot, NULL, shared_legend, NULL,
                     ncol = 4, rel_widths = c(1, 0.05, 0.18, 0.05))
}
