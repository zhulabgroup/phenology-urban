# city <- "NewYork"
# source <- "Grid" # "WU_buffer" or "Grid" or "WU_surface"
# taxa_col <- "canonicalname" # "canonicalname" or "genus"
# phe_type <- "SOS" # "SOS", "EOS"

prepare_data <- function(city, source = "Grid", taxa_col = "canonicalname", phe_type) {
  if (source == "Grid") { # UHAT and Daymet
       input_dir <- list.files(,
        path = paste0(here("release/input/tree_climate/"), source, "/"), # For example NY
        pattern = paste0("^", city, ".*\\.rds$"),
        full.names = TRUE)
       tree_climate <- read_rds(input_dir) |>
         select(id, taxa, lon, lat, year, temp_summer = UHAT_temp_summer, temp_winter = UHAT_temp_winter, 
                prcp_summer = Daymet_prcp_summer, prcp_winter = Daymet_prcp_winter,
                srad_summer = Daymet_srad_summer, srad_winter = Daymet_srad_winter) |>
         st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F) |>
         st_transform(crs = 5070)
    } else if (source == "Grid_GSL") { # UHAT and Daymet
      input_dir <- list.files(,
                              path = paste0(here("release/input/tree_climate/"), source, "/"), # For example NY
                              pattern = paste0("^", city, ".*\\.rds$"),
                              full.names = TRUE)
      tree_climate <- read_rds(input_dir) |>
        mutate(prcp_growing = 0.5* (Daymet_prcp_spring + Daymet_prcp_summer), 
               srad_growing = 0.5* (Daymet_srad_spring + Daymet_srad_summer)) |>
        select(id, taxa, lon, lat, year, 
               temp_growing = UHAT_temp_growing,
               prcp_growing, srad_growing) |>
        st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F) |>
        st_transform(crs = 5070)
    } else if (source == "Grid_NPN") { # UHAT and Daymet
      input_dir <- list.files(,
                              path = paste0(here("release/input/tree_climate/"), source, "/"), # For example NY
                              pattern = paste0("^", city, ".*\\.rds$"),
                              full.names = TRUE)
      tree_climate <- read_rds(input_dir) |>
        select(id, lon, lat, year, temp_summer = UHAT_temp_summer, temp_winter = UHAT_temp_winter, 
               prcp_summer = Daymet_prcp_summer, prcp_winter = Daymet_prcp_winter,
               srad_summer = Daymet_srad_summer, srad_winter = Daymet_srad_winter) |>
        st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F) |>
        st_transform(crs = 5070)
    } else {
      stop("Invalid source. Must be 'WU' or 'Grid'.")
    }
  
  #Join to phenology data

    file_list <- list.files(path = here("release/input/phenology/doy"), pattern = paste0("^doy_",city,"_.*\\.rds$"), full.names = TRUE)
    data_list <- file_list %>%
      purrr::map(~ read_rds(.x))
    all_doy <- bind_rows(data_list)
    tree_climate_phe <- tree_climate |>
      left_join(all_doy, by = c("id", "year"))

    if (taxa_col == "canonicalname") {
      model_input_all <- prepare_species_phe(phe_type, tree_climate_phe)
    } else {
      stop("Invalid taxa_col. Must be 'canonicalname' or 'genus'.")
    }
  return(model_input_all)
}


prepare_species_phe <- function(phe_type, tree_climate_phe) {
  if (phe_type == "SOS") {
    phe_metric <- "doy"
    direction <- "up"
    thres <- 0.5
    model_input <- tree_climate_phe %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, canonicalname = taxa, lon, lat, year, any_of(c("start", "end", "hex_id")), doy, ends_with("winter")) %>%
      drop_na()
  } else if (phe_type == "EOS") {
    phe_metric <- "doy"
    direction <- "down"
    thres <- 0.5
    model_input <- tree_climate_phe %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, canonicalname = taxa, lon, lat, year, any_of(c("start", "end", "hex_id")), doy, ends_with("summer")) %>%
      drop_na()
  } else if (str_detect(phe_type, "^SOS_\\d+$")) {
    phe_metric <- "doy"
    direction <- "up"
    thres <- phe_type |>
      str_extract("(?<=SOS_)\\d+") |>
      as.numeric() / 10
    model_input <- tree_climate_phe %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, canonicalname = taxa, lon, lat, year, start, end, doy, hex_id, ends_with("winter")) %>%
      drop_na()
  } else if (str_detect(phe_type, "^EOS_\\d+$")) {
    phe_metric <- "doy"
    direction <- "down"
    thres <- phe_type |>
      str_extract("(?<=EOS_)\\d+") |>
      as.numeric() / 10
    model_input <- tree_climate_phe %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, canonicalname = taxa, lon, lat, year, start, end, doy, hex_id, ends_with("summer")) %>%
      drop_na()
  }  else if (phe_type == "GSL_winter") {
    phe_metric <- "doy_diff"
    diff <- tree_climate_phe %>%
      filter(thres == 0.5) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = tree_climate_phe %>%
      filter(direction == "up" & thres == 0.5) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, canonicalname = taxa, lon, lat, year, doy_diff, hex_id, ends_with("winter")) %>%
      drop_na()
  } else if (phe_type == "GSL_summer") {
    phe_metric <- "doy_diff"
    diff <- tree_climate_phe %>%
      filter(thres == 0.5) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = tree_climate_phe %>%
      filter(direction == "down" & thres == 0.5) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, canonicalname = taxa, lon, lat, year, doy_diff, hex_id, ends_with("summer")) %>%
      drop_na()
  } else if (phe_type == "GSL_range_season") {
    phe_metric <- "doy_diff"
    diff <- tree_climate_phe %>%
      filter(thres == 0.5) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = tree_climate_phe %>%
      filter(direction == "down" & thres == 0.5) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, canonicalname = taxa, lon, lat, year, doy_diff, hex_id, ends_with("growing")) %>%
      drop_na()
  } else if (phe_type == "GreenUp") {
    phe_metric <- "doy_diff"
    direction <- "up"
    diff <- tree_climate_phe %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = tree_climate_phe %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, canonicalname = taxa, lon, lat, year, doy_diff, hex_id, ends_with("spring")) %>%
      drop_na()
  } else {
    direction <- "down"
    diff <- tree_climate_phe %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = tree_climate_phe %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, canonicalname = taxa, lon, lat, year, doy_diff, hex_id, ends_with("fall")) %>%
      drop_na()
  }
  return(model_input)
}
