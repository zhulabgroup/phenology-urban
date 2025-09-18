library(zhulabtools)
load_packages(c("tidyverse", "ggplot2", "sf", "spBayes", "INLA", "lubridate", "terra"))
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
  tree_dir <- paste0("~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/", city, ".rds")
  tree_location <- read_rds(tree_dir) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F)
  
  ################## Attribute climate to tree ##################
  message("Attributing climate to trees...")
  
  if (type == "WU_buffer"){
    out_path <- paste0("~/phenology-urban/data/proc/tree_climate/WU_buffer/", city, "_", buffer_dist, ".rds")
    WU_buffer <- WU_location |>
      distinct(new_id, lon, lat) |>
      st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F) |>
      st_buffer(WU_location, dist = 500)
  
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
      
    write_rds(points_in_buffer, out_path, compress = "gz")
    
  } else if (type == "WU_surface"){
    out_path <- paste0("~/phenology-urban/data/proc/tree_climate/WU_surface/", city, ".rds")
    
    seasons  <- c("summer","winter")
    var_map  <- c(Avg_AvgTemp = "temp", Sum_Sum_mm = "prcp")  
    proj_crs <- 5070                                           # NAD83 / Conus Albers
    log_prcp <- TRUE 

    clim_sf <- WU_climate |>
      inner_join(
        WU_location |> group_by(new_id) |>
          summarise(lon = median(lon, na.rm = TRUE), lat = median(lat, na.rm = TRUE), .groups = "drop"),
        by = "new_id"
      ) |>
      filter(year < 2025) |>
      st_as_sf(coords = c("lon","lat"), crs = 4326) |>
      st_transform(proj_crs)
    
    trees_sf <- st_transform(tree_location, proj_crs)
    
    spde_predict_one <- function(obs_sf, pred_sf, y_col,
                                 max_edge = c(5000, 20000), cutoff = 1000) {
      mesh <- inla.mesh.2d(
        loc     = rbind(st_coordinates(obs_sf), st_coordinates(pred_sf)),
        max.edge = max_edge, cutoff = cutoff
      )
      spde <- inla.spde2.matern(mesh, alpha = 2)
      Aobs  <- inla.spde.make.A(mesh, st_coordinates(obs_sf))
      Apred <- inla.spde.make.A(mesh, st_coordinates(pred_sf))
      Afix_obs  <- Diagonal(nrow(obs_sf), 1)
      Afix_pred <- Diagonal(nrow(pred_sf), 1)
      sidx <- inla.spde.make.index("spatial", spde$n.spde)
      stk_obs <- inla.stack(
        data    = list(y = obs_sf[[y_col]]),
        A       = list(Aobs, Afix_obs),
        effects = list(sidx, data.frame(Intercept = rep(1, nrow(obs_sf)))),
        tag     = "obs"
      )
      stk_prd <- inla.stack(
        data    = list(y = NA),
        A       = list(Apred, Afix_pred),
        effects = list(sidx, data.frame(Intercept = rep(1, nrow(pred_sf)))),
        tag     = "pred"
      )
      stk_all <- inla.stack(stk_obs, stk_prd)
      
      fit <- inla(y ~ -1 + Intercept + f(spatial, model = spde),
                  family = "gaussian",
                  data   = inla.stack.data(stk_all, spde = spde),
                  control.predictor = list(A = inla.stack.A(stk_all), compute = TRUE),
                  control.compute   = list(dic = FALSE, waic = FALSE, cpo = FALSE),
                  control.inla = list(int.strategy = "eb"),
                  verbose = FALSE)
      
      ipred <- inla.stack.index(stk_all, "pred")$data
      tibble(id = pred_sf$id,
             mean = fit$summary.fitted.values[ipred, "mean"],
             sd   = fit$summary.fitted.values[ipred, "sd"])
    }
    
    interpo_WU_climate <- expand_grid(year = sort(unique(clim_sf$year)), season = seasons, vcol = names(var_map)) |>
      mutate(
        out = pmap(
          list(year, season, vcol),
          function(yy, ssn, vc) {
            message(paste(yy, ssn, vc))
            obs <- clim_sf |> filter(year == yy, season == ssn)
            if (vc == "Sum_Sum_mm" && log_prcp) {
              obs[[vc]] <- log1p(obs[[vc]])
              spde_predict_one(obs_sf = obs, pred_sf = trees_sf, y_col = vc) |>
                mutate(value = pmax(expm1(mean), 0))
            } else {
              spde_predict_one(obs_sf = obs, pred_sf = trees_sf, y_col = vc) |>
                mutate(value = mean)
            }
          }
        )
      ) |>
      select(year, season, vcol, out) |>
      unnest(out) |>
      transmute(id, year,
                col = paste0(unname(var_map[vcol]), "_", season),  # e.g., temp_fall
                val = value) |>
      pivot_wider(names_from = col, values_from = val) |>
      arrange(id, year)
    
    points_surface <- tree_location |>
      left_join(daymet_climate |>
                  select(id, year, season, Avg_srad) |>
                  pivot_wider(names_from = season, values_from = Avg_srad, names_prefix = "srad_"), 
                by = c("id")) |>
      left_join(interpo_WU_climate, by = c("id", "year"))
    
    write_rds(points_surface, out_path, compress = "gz")
    
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
  

# Run for WU_surface
file_path <- path.expand("~/phenology-urban/data/raw/CONUS/WU")
sub_folders <- list.dirs(file_path, recursive = FALSE, full.names = FALSE) 
# example of running the function
for (city in sub_folders[54:78]){
  tryCatch(
    {
      prepare_tree_WU_climate(city, type = "WU_surface")
    },
    error = function(e) {
      message("Skipping ", city, ": ", e$message)
    }
  )}


##################################################################
# Test WU-surface for existing city

prepare_tree_Grid_climate <- function(city){
  ################## Calculate seasonal climate ##################
  message(paste0("Processing example city: ", city, "..."))
  message("Calculating seasonal climate...")
  
  ###### Daymet shortwave radiation #######
  daymet_srad <- read_csv(paste0("~/phenology-urban/data/raw/CONUS/Daymet/all_city_trees/", city, "/srad_monthly_2016_2024.csv"), show_col_types = FALSE)
  daymet_prcp <- read_csv(paste0("~/phenology-urban/data/raw/CONUS/Daymet/all_city_trees/", city, "/prcp_monthly_2016_2024.csv"), show_col_types = FALSE)
    
  ###### UHAT temperature #######
    
  uhat_files <- list.files(
    path = "~/lab-data/datasets/climate/UHAT/clean/",
    pattern = paste0("^", city, ".*\\.tif$"), 
    full.names = TRUE
  )
  missing_files <- uhat_files[!file.exists(uhat_files)]
  if (length(missing_files) > 0) {
    stop("The following CHELSA TIF file(s) do not exist: {str_c(missing_files, collapse = ', ')}")
  }
  uhat_climate <- map(uhat_files, function(file) {
    r <- rast(file)
    vals <- terra::extract(r, tree_location, ID = FALSE)
    
    fname <- basename(file)
    fname <- str_remove(fname, "\\.tif$")
    season <- str_extract(fname, "summer|winter")
    year   <- str_extract(fname, "\\d{4}")
    colnames(vals) <- "value"
    vals <- bind_cols(
      tibble(id = tree_location$id, year = year, season = season),
      vals
    )
    vals
  }) |> 
    bind_rows() |>
    pivot_wider(
      names_from  = season,
      values_from = value,
      names_glue  = "UHAT_temp_{season}"
    )
  
  points_surface <- tree_location |>
    left_join(daymet_srad |>
                select(id, year, season, Avg_srad) |>
                pivot_wider(names_from = season, values_from = Avg_srad, names_prefix = "Daymet_srad_"),
              by = c("id")) |>
    left_join(daymet_prcp |>
                select(id, year, season, Sum_Sum_mm) |>
                pivot_wider(names_from = season, values_from = Sum_Sum_mm, names_prefix = "Daymet_prcp_"),
              by = c("id", "year")) |>
    left_join(uhat_climate |>
                mutate(year = as.double(year)), by = c("id", "year"))
  
    write_rds(points_surface, out_path, compress = "gz")
}

prepare_tree_Grid_climate(city = "NewYork")
