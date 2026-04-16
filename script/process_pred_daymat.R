library(terra)
# download the Daymet data
library(daymetr)  # for getting Daymet data
library(sf)
library(raster)
library(ggplot2)
library(tidyverse)
library(lubridate)

############# city area #############
source(here("script/plot_phe_temp_margin.R"))
cities <- c("NewYork", "TP", "Stockton", "LosAngeles", "Seattle", "Sacramento", "WashingtonDC", "SanFrancisco")
phe_types = c("SOS", "SOS", "SOS", "SOS", "EOS", "EOS", "EOS", "EOS")
output_dir <- "~/phenology-urban/data/raw/CONUS/Daymet/example_pred_city/"

for (i in c(1:8)) {
  city <- cities[i]
  phe_type <- phe_types[i]
  
  if (phe_type == "SOS") {
    pre_season <- "winter"
  } else if (phe_type == "EOS") {
    pre_season <- "summer"
  } else {
    stop("Invalid phe_type. Use 'SOS' or 'EOS'.")
  }
  
  city_boundary <- search_boundary(city)
  file <- paste0("~/lab-data/datasets/climate/UHAT/clean/",city, "_", season, "_", year, ".tif")
  temp_tif <- terra::rast(file)
  tifs <- list.files(
    "~/lab-data/datasets/climate/UHAT/grid",
    pattern = "^MYD11A1\\..*LST_Day_1km\\.tif$",
    full.names = TRUE
  )
  r_template <- terra::rast(tifs[1])
  crs(temp_tif) <- crs(r_template)  # ensure correct CRS
  
  clipped_over <- temp_tif |>
    # crop(temp_tif, vect(tree_location)) |>
    project(crs(city_boundary), res = 0.009)  # transform to WGS84 for plotting
  city_v <- vect(city_boundary)
  clipped <- crop(clipped_over, city_v) |> mask(city_v)
  temp_df <- as.data.frame(clipped, xy = TRUE) |>
    rename_with(~ "value", .cols = 3) |>
    mutate(id = row_number()) |>
    select(id, x, y, value)
  
  message(paste0("Downloading Daymet data for ", city))
  data_dir <- paste0(output_dir, city,"/")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  results_prcp <- temp_df %>%
    pmap_dfr(~ {
      if (..1 %% 100 == 0) {
        message(paste0("Processing ID: ", ..1))
      }
      
      test <- tryCatch({
        download_daymet(
          ..1,  # id_site
          ..3,  # lon
          ..2,  # lat
          start = 2021,
          end = 2022,
          simplify = TRUE,
          internal = TRUE,
          silent = TRUE
        ) %>%
          filter(measurement == "prcp..mm.day.") %>%
          mutate(
            date = as.Date(yday - 1, origin = paste0(year, "-01-01")),
            id = ..1
          ) %>%
          mutate(
            season = case_when(
              month(date) %in% 3:5 ~ "spring",
              month(date) %in% 6:8 ~ "summer",
              month(date) %in% 9:11 ~ "fall",
              month(date) %in% c(12, 1, 2) ~ "winter"
            )
          ) %>%
          mutate(year = if_else(month(date) == 12, year + 1, year)) %>%
          group_by(id, year, season) %>%
          summarise(
            Sum_prcp = sum(value, na.rm = TRUE),
            prcp_count = sum(!is.na(value)),
            .groups = "drop"
          )
      }, error = function(e) {
        message("Failed ID: ", ..1, " | ", conditionMessage(e))
        tibble(
          id = ..1,
          year = NA_integer_,
          season = NA_character_,
          Sum_prcp = NA_real_,
          prcp_count = NA_integer_
        )
      })
      return(test)
    }) |>
    filter(year == 2022 & season == !!pre_season) |>
    select(id, year, season, prcp = Sum_prcp) |>
    left_join(temp_df, by = "id")
  
  results <- temp_df %>%
    pmap_dfr(~ {
      if (..1 %% 100 == 0) {
        message(paste0("Processing ID: ", ..1))
      }
      
      test <- tryCatch({
        download_daymet(
          ..1,  # id_site
          ..3,  # lon
          ..2,  # lat
          start = 2021,
          end = 2022,
          simplify = TRUE,
          internal = TRUE,
          silent = TRUE
        ) %>%
          filter(measurement == "srad..W.m.2.") %>%
          mutate(
            date = as.Date(yday - 1, origin = paste0(year, "-01-01")),
            id = ..1
          ) %>%
          mutate(
            season = case_when(
              month(date) %in% 3:5 ~ "spring",
              month(date) %in% 6:8 ~ "summer",
              month(date) %in% 9:11 ~ "fall",
              month(date) %in% c(12, 1, 2) ~ "winter"
            )
          ) %>%
          mutate(year = if_else(month(date) == 12, year + 1, year)) %>%
          group_by(id, year, season) %>%
          summarise(
            Avg_srad = mean(value, na.rm = TRUE),
            srad_count = sum(!is.na(value)),
            .groups = "drop"
          )
      }, error = function(e) {
        message("Failed ID: ", ..1, " | ", conditionMessage(e))
        tibble(
          id = ..1,
          year = NA_integer_,
          season = NA_character_,
          Avg_srad = NA_real_,
          prcp_count = NA_integer_
        )
      })
      return(test)
    }) |>
    filter(year == 2022 & season == !!pre_season) |>
    select(id, srad = Avg_srad) |>
    right_join(results_prcp, by = "id")
    
  output_file <- paste0(data_dir, "pred_monthly_2022.rds")
  write_rds(results, output_file, compress = "gz")
}
