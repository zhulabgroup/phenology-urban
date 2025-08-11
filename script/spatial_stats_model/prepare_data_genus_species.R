

# city <- "NY"
# source <- "WU" # "WU" or "Grid"
# taxa_col <- "canonicalname" # "canonicalname" or "genus"
# phe_type <- "SOS" # "SOS", "EOS", "GreenUp", "Dormancy"

prepare_data <- function(city, source, taxa_col, phe_type) {
  if (source == "WU") {
    points_in_buffer <- read_csv(paste0("~/phenology-urban/data/proc/urban/",city,"/tree_WU_500_buffer_PS.csv"), show_col_types = FALSE)
    
    points_sf <- st_as_sf(points_in_buffer, coords = c("lon", "lat"), crs = 4326, remove = F) |>
      st_transform(crs = 9311)
    hex_size <- 5000
    bbox <- st_bbox(points_sf)
    hex_grid_geom <- st_make_grid(bbox, cellsize = hex_size, square = FALSE)
    hex_grid <- st_sf(geometry = hex_grid_geom) |>
      mutate(hex_id = row_number()) |>
      mutate(centroid = st_centroid(geometry)) |>
      mutate(
        hex_lon = st_coordinates(centroid)[, 1],
        hex_lat = st_coordinates(centroid)[, 2]
      ) |>
      select(-centroid)
    # Assign each point to a hexagon
    points_in_buffer <- st_join(points_sf, hex_grid, join = st_within) %>%
      mutate(hex_id = as.factor(hex_id)) |>
      select(-hex_lon, -hex_lat) |>
      st_drop_geometry()
    
    all_sites_temp <- readRDS(paste0("~/urban-cooling/data/raw/WU/",city,"/",city,"_wu.rds")) %>%
      dplyr::select(c(Date, HighTemp, AvgTemp, LowTemp, Sum_mm, name)) %>%
      mutate(season = case_when(
        month(Date) %in% 3:5 ~ "spring",
        month(Date) %in% 6:8 ~ "summer",
        month(Date) %in% 9:11 ~ "fall",
        month(Date) %in% c(12, 1, 2) ~ "winter"
      )) %>%
      mutate(year = year(Date),
             year = if_else(month(Date) == 12, year + 1, year)) %>%
      group_by(name, year, season) %>%
      summarise(
        Avg_AvgTemp = mean(AvgTemp, na.rm = TRUE),
        Sum_Sum_mm = sum(Sum_mm, na.rm = TRUE),
        AvgTemp_count = sum(!is.na(AvgTemp)),
        Sum_mm_count = sum(!is.na(Sum_mm)),
        .groups = "drop"
      )
    
    srad_daily <- readRDS(paste0("~/phenology-urban/data/raw/",city,"/Daymet/daily2016-2024/srad_daily_2016-2023.rds")) %>%
      select(id, value, date)  %>%
      mutate(season = case_when(
        month(date) %in% 3:5 ~ "spring",
        month(date) %in% 6:8 ~ "summer",
        month(date) %in% 9:11 ~ "fall",
        month(date) %in% c(12, 1, 2) ~ "winter"
      )) %>%
      mutate(year = year(date),
             year = if_else(month(date) == 12, year + 1, year)) %>%
      group_by(id, year, season) %>%
      summarise(
        Avg_srad = mean(value, na.rm = TRUE),
        srad_count = sum(!is.na(value)),
        .groups = "drop"
      )
    if (taxa_col == "canonicalname") {
      model_input_all <- prepare_species_data(phe_type, points_in_buffer, all_sites_temp, srad_daily, source = "WU")
    } else if (taxa_col == "genus") {
      model_input_all <- prepare_genus_data(phe_type, points_in_buffer, all_sites_temp, srad_daily, source = "WU")
    } else {
      stop("Invalid taxa_col. Must be 'canonicalname' or 'genus'.")
    }
    
  } else if (source == "Grid") {
    metadata <- read_csv("~/lab-data/datasets/vegetation/PS/urban/metadata.csv") %>%
      filter(site == !!city)
    file_list <- list.files(path = "~/lab-data/datasets/vegetation/PS/urban/doy", pattern = paste0("^doy_",city,"_.*\\.rds$"), full.names = TRUE)
    data_list <- file_list %>%
      map(~ readRDS(.x))
    all_doy <- bind_rows(data_list)
    tree_location <- metadata %>%
      mutate(genus = str_extract(taxa, "^[^ ]+")) %>%
      left_join(all_doy, by = "id") %>%
      filter(!is.na(doy)) %>%
      mutate(geometry = paste0("c(", lon, ", ", lat, ")"))
    
    points_sf <- st_as_sf(tree_location |> distinct(lon, lat), coords = c("lon", "lat"), crs = 4326, remove = F) |>
      st_transform(crs = 9311)
    hex_size <- 5000
    bbox <- st_bbox(points_sf)
    hex_grid_geom <- st_make_grid(bbox, cellsize = hex_size, square = FALSE)
    hex_grid <- st_sf(geometry = hex_grid_geom) |>
      mutate(hex_id = row_number()) |>
      mutate(centroid = st_centroid(geometry)) |>
      mutate(
        hex_lon = st_coordinates(centroid)[, 1],
        hex_lat = st_coordinates(centroid)[, 2]
      ) |>
      select(-centroid)
    # Assign each point to a hexagon
    tree_location_hex <- st_join(points_sf, hex_grid, join = st_within) %>%
      mutate(hex_id = as.factor(hex_id)) |>
      select(-hex_lon, -hex_lat) |>
      st_drop_geometry() |>
      right_join(tree_location, by = c("lon", "lat"))
    
    # List of Humid variable names and their corresponding filenames
    files <- c("temp_monthly_all_tree.rds",
               "srad_monthly_all_tree.rds")
    base_path <- "~/phenology-urban/data/raw/"
    processed_humid_list <- lapply(files, function(file_name) {
      # Read the file
      var_name <- sub("_monthly_all_tree.rds", "", file_name)
      data <- readRDS(paste0(base_path, city, "/HUMID/", file_name)) %>%
        select(id, value, date) %>%
        mutate(
          season = case_when(
            substr(date, 6, 7) %in% c("03", "04", "05") ~ "spring",
            substr(date, 6, 7) %in% c("06", "07", "08") ~ "summer",
            substr(date, 6, 7) %in% c("09", "10", "11") ~ "fall",
            substr(date, 6, 7) %in% c("12", "01", "02") ~ "winter"
          )) %>%
        mutate(year = as.integer(substr(date, 1, 4)),
               year = if_else(substr(date, 6, 7) == "12", year + 1, year)) %>%
        filter(year < 2019) %>% 
        group_by(id, year, season) %>%
        summarise(
          # Use sum for prcp, and mean for temp and srad
          Avg_value = mean(as.numeric(value), na.rm = TRUE),
          .groups = "drop"
        ) %>%
        # Rename the Avg_value column to the variable name (e.g., tmax, tmin, prcp)
        rename(!!var_name := Avg_value)
      return(data)
    })
    
    files <- c("prcp_monthly_all_tree.rds")
    processed_daymet_list <- lapply(files, function(file_name) {
      # Read the file
      var_name <- sub("_monthly_all_tree.rds", "", file_name)
      data <- readRDS(paste0(base_path, city, "/Daymet/monthly2016-2024/", file_name)) %>%
        select(id, value, date) %>%
        mutate(
          season = case_when(
            substr(date, 6, 7) %in% c("03", "04", "05") ~ "spring",
            substr(date, 6, 7) %in% c("06", "07", "08") ~ "summer",
            substr(date, 6, 7) %in% c("09", "10", "11") ~ "fall",
            substr(date, 6, 7) %in% c("12", "01", "02") ~ "winter"
          )) %>%
        mutate(year = as.integer(substr(date, 1, 4)),
               year = if_else(substr(date, 6, 7) == "12", year + 1, year)) %>%
        group_by(id, year, season) %>%
        summarise(
          # Use sum for prcp, and mean for temp and srad
          Avg_value = sum(as.numeric(value), na.rm = TRUE),
          .groups = "drop"
        ) %>%
        # Rename the Avg_value column to the variable name (e.g., tmax, tmin, prcp)
        rename(!!var_name := Avg_value)
      
      return(data)
    })
    
    all_sites_climate <- processed_humid_list[[1]] |>
      left_join(processed_daymet_list[[1]], by = c("id", "season", "year")) 
    
    srad_daily <- processed_humid_list[[2]] |>
      rename(Avg_srad = srad)
    if (taxa_col == "canonicalname") {
      model_input_all <- prepare_species_data(phe_type, tree_location_hex, all_sites_climate, srad_daily, source = "Grid")
    } else if (taxa_col == "genus") {
      model_input_all <- prepare_genus_data(phe_type, tree_location_hex, all_sites_climate, srad_daily, source = "Grid")
    } else {
      stop("Invalid taxa_col. Must be 'canonicalname' or 'genus'.")
    }
    
  } else {
    stop("Invalid source. Must be 'WU' or 'Grid'.")
  }
  return(model_input_all)
}

  

prepare_species_data <- function(phe_type, points_in_buffer, all_sites_temp, srad_daily, source) {
  if (source == "WU") {
    points_in_buffer_joined <- points_in_buffer %>%
      left_join((srad_daily %>%
                   select(id, year, season, Avg_srad) %>%
                   pivot_wider(names_from = season, values_from = Avg_srad, names_prefix = "srad_")),
                by = c("id", "year")) %>%
      left_join(all_sites_temp %>%
                  select(name, year, season, Avg_AvgTemp) %>%
                  pivot_wider(names_from = season, values_from = Avg_AvgTemp, names_prefix = "temp_"),
                by = c("Site" = "name", "year")) %>%
      left_join(all_sites_temp %>%
                  select(name, year, season, Sum_Sum_mm) %>%
                  pivot_wider(names_from = season, values_from = Sum_Sum_mm, names_prefix = "prcp_"),
                by = c("Site" = "name", "year"))
  } else if (source == "Grid") {
    ## All Daymet
    points_in_buffer_joined <- points_in_buffer %>%
      left_join((srad_daily %>%
                   select(id, year, season, Avg_srad) %>%
                   pivot_wider(names_from = season, values_from = Avg_srad, names_prefix = "srad_")),
                by = c("id", "year")) %>%
      left_join(all_sites_temp %>%
                  mutate(Avg_AvgTemp = temp) %>%
                  select(id, year, season, Avg_AvgTemp) %>%
                  pivot_wider(names_from = season, values_from = Avg_AvgTemp, names_prefix = "temp_"),
                by = c("id", "year")) %>%
      left_join(all_sites_temp %>%
                  select(id, year, season, prcp) %>%
                  pivot_wider(names_from = season, values_from = prcp, names_prefix = "prcp_"),
                by = c("id", "year"))
  }
  
  if (phe_type == "SOS") {
    phe_metric <- "doy"
    direction <- "up"
    thres <- 0.5
    model_input <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, taxa, lon, lat, year, doy, hex_id, ends_with("winter"), Site) %>%
      drop_na()
  } else if (phe_type == "EOS") {
    phe_metric <- "doy"
    direction <- "down"
    thres <- 0.5
    model_input <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, taxa, lon, lat, year, doy, hex_id, ends_with("summer"), Site) %>%
      drop_na()
  } else if (phe_type == "GreenUp") {
    phe_metric <- "doy_diff"
    direction <- "up"
    diff <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = points_in_buffer_joined %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, taxa, lon, lat, year, doy_diff, hex_id, ends_with("spring")) %>%
      drop_na()
    
  } else {
    direction <- "down"
    diff <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = points_in_buffer_joined %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, taxa, lon, lat, year, doy_diff, hex_id, ends_with("fall")) %>%
      drop_na()
  }
  
  gbif_taxon_keys <- model_input$taxa %>% # get the first column of the data.frame
    unique()  %>% # remove duplicates
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta") %>%
    filter(kingdom == "Plantae", phylum == "Tracheophyta") |>
    filter(status %in% c("ACCEPTED", "SYNONYM") & rank == "SPECIES") %>% # get only accepted names and remove anything that might have matched to a non-species
    select(canonicalName, verbatim_name) # check for missing values
  
  model_input_top_10_species <- model_input %>%
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
  return(model_input_top_10_species)
}

prepare_genus_data <- function(phe_type, points_in_buffer, all_sites_temp, srad_daily, source) {
  if (source == "WU") {
    ### WU sites
    points_in_buffer_joined <- points_in_buffer %>%
      left_join((srad_daily %>%
                   select(id, year, season, Avg_srad) %>%
                   pivot_wider(names_from = season, values_from = Avg_srad, names_prefix = "srad_")),
                by = c("id", "year")) %>%
      left_join(all_sites_temp %>%
                  select(name, year, season, Avg_AvgTemp) %>%
                  pivot_wider(names_from = season, values_from = Avg_AvgTemp, names_prefix = "temp_"),
                by = c("Site" = "name", "year")) %>%
      left_join(all_sites_temp %>%
                  select(name, year, season, Sum_Sum_mm) %>%
                  pivot_wider(names_from = season, values_from = Sum_Sum_mm, names_prefix = "prcp_"),
                by = c("Site" = "name", "year"))
  } else if (source == "Grid") {
    ## All Daymet
    points_in_buffer_joined <- points_in_buffer %>%
      left_join((srad_daily %>%
                   select(id, year, season, Avg_srad) %>%
                   pivot_wider(names_from = season, values_from = Avg_srad, names_prefix = "srad_")),
                by = c("id", "year")) %>%
      left_join(all_sites_temp %>%
                  mutate(Avg_AvgTemp = temp) %>%
                  select(id, year, season, Avg_AvgTemp) %>%
                  pivot_wider(names_from = season, values_from = Avg_AvgTemp, names_prefix = "temp_"),
                by = c("id", "year")) %>%
      left_join(all_sites_temp %>%
                  select(id, year, season, prcp) %>%
                  pivot_wider(names_from = season, values_from = prcp, names_prefix = "prcp_"),
                by = c("id", "year"))
  } else {
    stop("Invalid source. Must be 'WU' or 'Daymet'.")
  }
  
  if (phe_type == "SOS") {
    phe_metric <- "doy"
    direction <- "up"
    thres <- 0.5
    model_input <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, taxa,genus, lon, lat, year, doy, hex_id,  ends_with("winter"), Site) %>%
      drop_na()
  } else if (phe_type == "EOS") {
    phe_metric <- "doy"
    direction <- "down"
    thres <- 0.5
    model_input <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, taxa, genus, lon, lat, year, doy, hex_id,  ends_with("summer"), Site) %>%
      drop_na()
  } else if (phe_type == "GreenUp") {
    phe_metric <- "doy_diff"
    direction <- "up"
    diff <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = points_in_buffer_joined %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, taxa, genus, lon, lat, year, doy_diff, hex_id, ends_with("spring")) %>%
      drop_na()
    
  } else {
    direction <- "down"
    diff <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = points_in_buffer_joined %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, taxa, genus, lon, lat, year, doy_diff, hex_id, ends_with("fall")) %>%
      drop_na()
  }
  
  return(model_input)
}

prepare_only_phe <- function(city, phe_type, taxon_col, summary = TRUE){
  metadata <- read_csv("~/lab-data/datasets/vegetation/PS/urban/metadata.csv",  show_col_types = FALSE) %>%
    filter(site == !!city)
  file_list <- list.files(path = "~/lab-data/datasets/vegetation/PS/urban/doy", pattern = paste0("^doy_",city,"_.*\\.rds$"), full.names = TRUE)
  data_list <- file_list %>%
    map(~ readRDS(.x))
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
