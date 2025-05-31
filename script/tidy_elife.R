library(dplyr)
library(rgbif)

##################### read in the data #####################

cleaning <- FALSE

if (cleaning) {
  city_polygons <- st_read(here("data/raw/CONUS/boundary/USCities/tl_2024_us_uac20.shp")) |>
    separate(NAME20, into = c("city_raw", "state"), sep = ", ", remove = FALSE) |>
    # Split rows if city contains '--'
    mutate(city_list = str_split(city_raw, "--")) |>
    unnest(city_list) |>
    # Clean up column names
    rename(city = city_list) |>
    # remove space in the name
    mutate(city = str_replace_all(city, "\\s+", "")) |>
    select(NAME20, state, city, ALAND20,AWATER20, INTPTLAT20, INTPTLON20) |>
    mutate(city = if_else(NAME20 == "Washington--Arlington, DC--VA--MD" & city == "Washington",
                          "WashingtonDC",
                          city))
  
  
  file_path <- here("data/raw/CONUS/StreetTrees/elife/v2")
  tr_inv_files <- list.files(path = file_path, pattern = "\\.csv$", full.names = TRUE)
  
  tr_inv <- do.call(rbind,
                    Filter(Negate(is.null),
                           lapply(tr_inv_files, function(file) {
                             data <- read.csv(file)
                             
                             required_columns <- c("city", "state", "scientific_name", "common_name", 
                                                   "longitude_coordinate", "latitude_coordinate",
                                                   "condition", "native", "most_recent_observation")
                             
                             if (!all(required_columns %in% colnames(data))) {
                               message("Skipping file due to missing columns: ", file)
                               return(NULL)
                             }
                             
                             dplyr::select(data, all_of(required_columns)) |>
                               mutate(source_file = basename(file)) 
                           })))
  
  tr_inv <- tr_inv |>
    mutate(scientific_name = stri_replace_all_regex(scientific_name, "[^\\x20-\\x7E]", "")) |>
    mutate(scientific_name = str_to_lower(scientific_name)) |> 
    mutate(scientific_name = str_replace(scientific_name, "^(\\w)(\\w*)", function(m) {
      paste0(str_to_upper(substr(m, 1, 1)), substr(m, 2, nchar(m)))
    }))
  
  ##################### Filter out records with condition, coordinates species name #####################
  
  elife_sp_list <- tr_inv |>
    distinct(scientific_name) |>
    arrange(scientific_name) |>
    pull(scientific_name)
  # 1st divide the genus name/ species name
  
  elife_sp_backbone <- elife_sp_list |>
    name_backbone_checklist(kingdom = "Plantae", phylum = "Tracheophyta")
  
  elife_sp_valid <- elife_sp_backbone |>
    filter(rank %in% c("SPECIES", "GENUS"))
  
  replacement_table <- elife_sp_valid |>
    mutate(
      replacement = verbatim_name,
      replacement = if_else(matchType == "FUZZY", canonicalName, replacement),
      replacement = if_else(status == "SYNONYM" & rank == "GENUS", genus, replacement), 
      replacement = if_else(status == "SYNONYM" & rank == "SPECIES", species, replacement) 
    ) |>
    dplyr::select(verbatim_name, replacement, rank)
  
  elife_trees_with_dup <- tr_inv |>
    left_join(replacement_table, by = c("scientific_name" = "verbatim_name")) |>
    mutate(scientific_name = ifelse(is.na(replacement), scientific_name, replacement)) |>
    dplyr::select(-replacement, -rank)
  
  saveRDS(elife_trees_with_dup, file = here("data/raw/CONUS/StreetTrees/elife/elife_trees_with_dup.rds"))
  
  ##################### Filter out records with valid coordinates #####################
  
  single_coords <- elife_trees_with_dup |>
    mutate(most_recent_observation = as.Date(most_recent_observation, format = "%m/%d/%Y")) |>
    group_by(longitude_coordinate, latitude_coordinate, scientific_name) |>
    filter(n() == 1)
  
  # exactly same information
  duplicate_coords <- elife_trees_with_dup |>
    mutate(condition = factor(condition, levels = c("excellent", "good", "fair", "poor", "dead/dying", "dead"), ordered = TRUE)) |>
    mutate(most_recent_observation = as.Date(most_recent_observation, format = "%m/%d/%Y")) |>
    group_by(longitude_coordinate, latitude_coordinate, scientific_name) |>
    filter(n() > 1) |>
    group_by(longitude_coordinate, latitude_coordinate, scientific_name) |>
    arrange(desc(most_recent_observation), condition) |> #(NA in time, ignore)
    distinct(longitude_coordinate, latitude_coordinate, scientific_name, .keep_all = TRUE)
  
  elife_trees_clean <- rbind(single_coords, duplicate_coords) |>
    filter(is.na(longitude_coordinate) == FALSE & is.na(latitude_coordinate) == FALSE) |> 
    mutate(city_combine = str_split(source_file, "_", simplify = TRUE)[,1]) |>
    mutate(city_combine = if_else(city_combine == "AuroraCO", "Aurora", city_combine),
           city_combine = if_else(city_combine == "StLouis", "St.Louis", city_combine)) |>
    mutate(state_abb = state.abb[match(state,state.name)],
           state_abb = replace_na(state_abb, "DC")) 
  
  elife_trees_sf <- elife_trees_clean %>%
    st_as_sf(coords = c("longitude_coordinate", "latitude_coordinate"), crs = 4326)
  
  elife_city <- elife_trees_clean %>%
    group_by(city_combine, state_abb) |>
    summarise(
      TreeCount = n(),
      .groups = "drop"
    ) 
  
  possible_polygons <- elife_city |>
    tidyr::crossing(city_polygons) |>
    filter(
      str_detect(city, fixed(city_combine, ignore_case = TRUE)),
      str_detect(state, fixed(state_abb, ignore_case = TRUE))
    ) |>
    st_as_sf(crs = 4326)
  
  elife_trees_clean <- elife_trees_sf |>
    st_transform(3857) |>
    group_by(city_combine) |>
    group_split() |>  # split into list of data frames by city_combine
    map_dfr(function(tree_group) {
      city_id <- unique(tree_group$city_combine)[1]
      
      # Get corresponding polygon
      poly_match <- possible_polygons |>
        st_transform(3857) |>
        filter(city_combine == city_id)
      
      if (nrow(poly_match) == 0) return(NULL)  # skip if no polygon match
      message("Matching trees in ", city_id)
      st_join(tree_group, poly_match, join = st_within)
    })
  
  saveRDS(elife_trees_clean, file = here("data/raw/CONUS/StreetTrees/elife/elife_trees_clean.rds"))
}

tidy_elife <- function(city_polygons) {
  elife_trees_clean <- readRDS(here("data/raw/CONUS/StreetTrees/elife/elife_trees_clean.rds")) |>
    filter(!is.na(NAME20)) |>
    st_transform(4326)
  
  elife_city <- elife_trees_clean |>
    group_by(city_combine.x) |>
    summarise(
      TreeCount = n(),
      geometry = st_convex_hull(st_union(geometry)),
      .groups = "drop"
    )
  
  return(elife_city)
}






