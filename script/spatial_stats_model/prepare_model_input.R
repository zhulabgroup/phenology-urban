# city <- "NewYork"
# source <- "Grid" # "WU_buffer" or "Grid" or "WU_surface"
# taxa_col <- "canonicalname" # "canonicalname" or "genus"
# phe_type <- "SOS" # "SOS", "EOS"

prepare_data <- function(city, source, taxa_col, phe_type) {
  if (source == "WU_buffer") {
    input_dir <- list.files(
      path = "~/phenology-urban/data/proc/tree_climate/WU_buffer/",
      pattern = paste0("^", city, "_.*\\.rds$"),
      full.names = TRUE)
    tree_climate <- read_rds(input_dir) |>
      st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F) |>
      st_transform(crs = 5070)
    } else if (source == "WU_surface") {
      input_dir <- list.files(
        # path = "~/phenology-urban/data/proc/tree_climate/WU_surface/",
        path = paste0("~/phenology-urban/data/proc/urban/", city, "/WU_surface/"), # For example NY
        pattern = paste0("^", city, ".*\\.rds$"),
        full.names = TRUE)
      tree_climate <- read_rds(input_dir) |>
        st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F) |>
        st_transform(crs = 5070)
    } else if (source == "Grid") { # UHAT and Daymet
       input_dir <- list.files(,
        path = paste0("~/phenology-urban/data/proc/tree_climate/", source, "/"), # For example NY
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
                              path = paste0("~/phenology-urban/data/proc/tree_climate/", source, "/"), # For example NY
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
    } else {
      stop("Invalid source. Must be 'WU' or 'Grid'.")
    }
    
    # Assign each point to a hexagon
    xrange = diff(range(st_coordinates(tree_climate)[,1]))
    yrange = diff(range(st_coordinates(tree_climate)[,2]))
    diag_len = sqrt(xrange^2 +yrange^2)
    hex_size <- diag_len/15
    bbox <- st_bbox(tree_climate)
    
    hex_grid_geom <- st_make_grid(bbox, cellsize = hex_size, square = FALSE)
    hex_grid <- st_sf(geometry = hex_grid_geom) |>
      mutate(hex_id = row_number()) |>
      mutate(centroid = st_centroid(geometry)) |>
      mutate(
        hex_lon = st_coordinates(centroid)[, 1],
        hex_lat = st_coordinates(centroid)[, 2]
      ) |>
      select(-centroid)
    tree_climate_hex <- st_join(tree_climate, hex_grid, join = st_within) %>%
      mutate(hex_id = as.factor(hex_id)) |>
      select(-hex_lon, -hex_lat) |>
      st_drop_geometry()
    
    #Join to phenology data
    file_list <- list.files(path = "~/lab-data/datasets/vegetation/PS/CONUS_urban/doy/", pattern = paste0("^doy_",city,"_.*\\.rds$"), full.names = TRUE)
    data_list <- file_list %>%
      purrr::map(~ read_rds(.x))
    all_doy <- bind_rows(data_list)
    tree_climate_phe <- tree_climate_hex |>
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
      select(id, canonicalname = taxa, lon, lat, year, start, end, doy, hex_id, ends_with("winter")) %>%
      drop_na()
  } else if (phe_type == "EOS") {
    phe_metric <- "doy"
    direction <- "down"
    thres <- 0.5
    model_input <- tree_climate_phe %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, canonicalname = taxa, lon, lat, year, start, end, doy, hex_id, ends_with("summer")) %>%
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

prepare_only_phe <- function(city, phe_type, taxon_col, summary = TRUE){
  metadata <- read_csv("~/lab-data/datasets/vegetation/PS/urban/metadata.csv",  show_col_types = FALSE) %>%
    filter(site == !!city)
  file_list <- list.files(path = "~/lab-data/datasets/vegetation/PS/CONUS_urban/doy", pattern = paste0("^doy_",city,"_.*\\.rds$"), full.names = TRUE)
  data_list <- file_list %>%
    purrr::map(~ readRDS(.x))
  all_doy <- bind_rows(data_list)
  tree_location <- metadata %>%
    mutate(genus = str_extract(taxa, "^[^ ]+")) %>%
    left_join(all_doy, by = "id") %>%
    filter(!is.na(doy)) %>%
    mutate(geometry = paste0("c(", lon, ", ", lat, ")"))
  if (phe_type == "SOS") {
    phe_metric <- "doy"
    direction <- "up"
    thres <- 0.5
    model_input <- tree_location %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, genus, taxa, lon, lat, year, doy, ends_with("winter")) %>%
      drop_na()
  } else if (phe_type == "EOS") {
    phe_metric <- "doy"
    direction <- "down"
    thres <- 0.5
    model_input <- tree_location %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, genus, taxa, lon, lat, year, doy, ends_with("summer")) %>%
      drop_na()
  }
  gbif_taxon_keys <- model_input$taxa %>% # get the first column of the data.frame
    unique()  %>% # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta") %>%
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES") %>% # get only accepted names and remove anything that might have matched to a non-species
    select(canonicalName, verbatim_name) # check for missing values
  phe_clean_species <- model_input %>%
    left_join(gbif_taxon_keys, by = c("taxa" = "verbatim_name")) %>% # join the GBIF taxon keys
    filter(canonicalName %in% (model_input %>%
                                 left_join(gbif_taxon_keys, by = c("taxa" = "verbatim_name")) %>%
                                 drop_na(canonicalName) %>% 
                                 group_by(canonicalName) %>%
                                 summarise(n = n(), .groups = "drop") %>%
                                 filter(n >= 30) %>%
                                 # arrange(desc(n)) %>%
                                 # head(10) %>% # get the top 10 species, otherwise all sp
                                 pull(canonicalName))) |> # filter to the taxa in the original data)
    rename(canonicalname = canonicalName)
  
  if (summary) {
    phe_clean_species <- phe_clean_species %>%
      group_by(.data[[taxon_col]], year) %>%
      summarise(
        n = n(),
        median_doy = median(doy, na.rm = TRUE),
        IQR_d = quantile(doy, 0.75, na.rm = TRUE) - quantile(doy, 0.25, na.rm = TRUE),
        lower_whisker_d = max(quantile(doy, 0.25, na.rm = TRUE) - 1.5 * IQR_d, min(doy, na.rm = TRUE)),
        upper_whisker_d = min(quantile(doy, 0.75, na.rm = TRUE) + 1.5 * IQR_d, max(doy, na.rm = TRUE)),
        whisker_d = (upper_whisker_d - lower_whisker_d),
        sd = sd(doy, na.rm = TRUE),
        .groups = "drop"
      )
  }
  return(phe_clean_species)
}
