library(zhulabtools)
load_packages(c("tidyverse", "ggplot2", "sf", "spBayes", "INLA", "lubridate"))
select <- dplyr::select


prepare_tree_WU_climate <- function(city, type, buffer_dist = NULL){
  #type %in% c(WU_buffer, WU_surface)
  if (type == "WU_buffer" && is.null(buffer_dist)) {
    stop("Argument `buffer_dist` must be provided when type = 'WU_buffer'.")
  }
  
  ################## Calculate seasonal climate ##################
  message(paste0("Processing ", city, "..."))
  message("Calculating seasonal climate...")
  
  ###### Daymet shortwave radiation #######
  daymet_dir <- paste0("~/phenology-urban/data/raw/CONUS/Daymet/all_city_trees/", city, "/srad_monthly_2016_2024.csv")
  daymet_climate <- read_csv(daymet_dir, show_col_types = FALSE)
  
  ###### WU temperature + precipitation #######
  WU_dir <- paste0("~/phenology-urban/data/raw/CONUS/WU_clean/", city, "_wu.rds")
  WU_loc_dir <- paste0("~/phenology-urban/data/raw/CONUS/WU/", city, "/location.csv")
  
  WU_location <- read_csv(WU_loc_dir, show_col_types = FALSE) |>
    filter(!is.na(lon) & !is.na(lat)) |>
    distinct(lon, lat) |>
    mutate(new_id = paste0(str_to_upper(city), row_number())) |>
    right_join(read_csv(WU_loc_dir, show_col_types = FALSE)|>
                 filter(!is.na(lon) & !is.na(lat)), 
               by = c("lon", "lat")) |>
    select(name = station, lon, lat, new_id)
  
  WU_climate <- read_rds(WU_dir) |>
    dplyr::select(c(Date, HighTemp, AvgTemp, LowTemp, Sum_mm, name)) |>
    mutate(season = case_when(
      month(Date) %in% 3:5 ~ "spring",
      month(Date) %in% 6:8 ~ "summer",
      month(Date) %in% 9:11 ~ "fall",
      month(Date) %in% c(12, 1, 2) ~ "winter"
    )) |>
    mutate(year = year(Date),
           year = if_else(month(Date) == 12, year + 1, year)) |>
    left_join(WU_location, by = "name") |>
    group_by(name, new_id, year, season) |>
    summarise(
      station_AvgTemp = mean(AvgTemp, na.rm = TRUE),
      station_Sum_mm = sum(Sum_mm, na.rm = TRUE),
      station_AvgTemp_count = sum(!is.na(AvgTemp)),
      station_Sum_mm_count = sum(!is.na(Sum_mm)),
      .groups = "drop"
    ) |>
    group_by(new_id, year, season) |>
    summarise(
      Avg_AvgTemp = mean(station_AvgTemp, na.rm = TRUE),
      Sum_Sum_mm = mean(station_Sum_mm, na.rm = TRUE),
      AvgTemp_count = sum(station_AvgTemp_count),
      Sum_mm_count = sum(station_AvgTemp_count),
      .groups = "drop"
    )
    
  
  ################## Attribute climate to tree ##################
  message("Attributing climate to trees...")
  
  if (type == "WU_buffer"){
    out_path <- paste0("~/phenology-urban/data/proc/tree_climate/WU_buffer/", city, "_", buffer_dist, ".rds")
    tree_dir <- paste0("~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/", city, ".rds")
    
    tree_location <- read_rds(tree_dir) |>
      st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F)
    WU_buffer <- WU_location |>
      distinct(new_id, lon, lat) |>
      st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F) |>
      st_buffer(WU_location, dist = 500)
    
    ggplot() +
      # geom_sf(data = boundary, fill = "grey") +
      geom_sf(data = WU_location |>
                distinct(new_id, lon, lat) |>
                st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F), color = "blue", size = 1) +
      geom_sf(data = WU_buffer, fill = "blue", color = NA, alpha = 0.2) +
      geom_sf(data = tree_location, color = "red", size = 0.01) +
      theme_minimal()
    
    
    intersects <- st_intersects(tree_location, WU_buffer)
    idx_tree   <- seq_len(nrow(tree_location))
    
    # pick one buffer per tree: if 1 hit -> that one; if >1 -> nearest among those; if 0 -> NA
    picked_buffer_idx <- map2_int(intersects, idx_tree, function(hits, i) {
      if (length(hits) == 0) return(NA_integer_)
      if (length(hits) == 1) return(hits[1])
      d <- st_distance(st_geometry(tree_location)[i], st_geometry(WU_buffer)[hits])
      hits[which.min(as.numeric(d))]
    })
    
    id_buffer <- tibble(
      id        = tree_location$id,
      buffer_id = picked_buffer_idx
    ) |>
      left_join(
        WU_buffer |> st_drop_geometry() |> mutate(buffer_id = row_number()),
        by = "buffer_id"
      ) |>
      filter(!is.na(buffer_id)) |>
      transmute(id, new_id)
    
    points_in_buffer <- tree_location |>
      left_join(id_buffer, by = "id") |>
      left_join(daymet_climate |>
                  select(id, year, season, Avg_srad) |>
                  pivot_wider(names_from = season, values_from = Avg_srad, names_prefix = "srad_"), 
                by = c("id")) |>
      left_join(WU_climate |>
                  select(new_id, year, season, Avg_AvgTemp) |>
                  pivot_wider(names_from = season, values_from = Avg_AvgTemp, names_prefix = "temp_"),
                by = c("new_id", "year")) |>
      left_join(WU_climate |>
                  select(new_id, year, season, Sum_Sum_mm) |>
                  pivot_wider(names_from = season, values_from = Sum_Sum_mm, names_prefix = "prcp_"),
                by = c("new_id", "year")) |>
      filter(!is.na(new_id))
      
    write_rds(points_in_buffer, out_path)
    
  } else if (source == "WU_surface"){

    
    
    
  } else {
    stop("source must be either 'WU' or 'WU_surface'")
  }
}



##################################################################
# Run for WU_buffer
file_path <- path.expand("~/phenology-urban/data/raw/CONUS/WU")
sub_folders <- list.dirs(file_path, recursive = FALSE, full.names = FALSE) 
# example of running the function
for (city in sub_folders){
  tryCatch(
    {
      prepare_tree_WU_climate(city, type = "WU_buffer", buffer_dist = 500)
    },
    error = function(e) {
      message("Skipping ", city, ": ", e$message)
    }
  )
}
  