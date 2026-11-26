library(zhulabtools)
load_packages(c("here", "tidyverse", "ggplot2", "sf", "lubridate", "terra", "patchwork", "scales", "viridis", "leaflet"))
select <- dplyr::select
pacman::p_load(terra)
source(here(
  "script/prepare_tree_Grid_climate_adaptive.R"
))

summarise_Daymet_seasons <- function(daymet, season_type = "full_season") {
  season_type <- match.arg(season_type, c("full_season", "early_season", "later_season"))
  winter_months <- switch(season_type,
                          full_season = c(12, 1, 2), early_season = c(12, 1), later_season = c(1, 2)
  )
  summer_months <- switch(season_type,
                          full_season = c(6, 7, 8), early_season = c(6, 7), later_season = c(7, 8)
  )
  if (!all(c("id", "date", "srad", "prcp") %in% names(daymet))) {
    stop("Daily Daymet data must contain id, date, srad, and prcp")
  }
  daymet <- daymet |>
    mutate(date = as.Date(date))
  if (anyNA(daymet$id) || anyNA(daymet$date)) stop("Daily Daymet data contains missing IDs or dates")
  if (anyDuplicated(daymet[c("id", "date")])) stop("Daily Daymet data contains duplicate tree dates")
  
  daymet |>
    mutate(month = lubridate::month(date),
           year = lubridate::year(date) + as.integer(month == 12),
           season = case_when(
             month %in% winter_months ~ "winter",
             month %in% summer_months ~ "summer",
             season_type == "full_season" & month %in% 3:5 ~ "spring",
             season_type == "full_season" & month %in% 9:11 ~ "fall",
             TRUE ~ NA_character_
           )) |>
    filter(!is.na(season)) |>
    group_by(id, year, season) |>
    summarise(
      Avg_srad = if (all(is.na(srad))) NA_real_ else mean(srad, na.rm = TRUE),
      Sum_Sum_mm = if (all(is.na(prcp))) NA_real_ else sum(prcp, na.rm = TRUE),
      .groups = "drop"
    )
}

prepare_tree_Grid_climate <- function(city, raw_uhat, daymet_raw, fixed_season = TRUE, season_type = "full_season", preseason_length = 30, phenophase = "SOS", ...) {
  if (!is.logical(fixed_season) || length(fixed_season) != 1L || is.na(fixed_season)) {
    stop("fixed_season must be TRUE or FALSE")
  }
  if (!fixed_season) {
    return(
      prepare_tree_Grid_climate_adaptive(
        city = city, raw_uhat, daymet_raw,
        preseason_length = preseason_length,
        phenophase = phenophase,
        ...
      )
    )
  }
  season_type <- match.arg(season_type, c("full_season", "early_season", "later_season"))
  output_dir <- switch(season_type, full_season = "Grid", early_season = "Grid_early_season",
                       later_season = "Grid_later_season")
  out_path <- here("data/proc/tree_climate", output_dir, paste0(city, ".rds"))
  
  ################## Calculate seasonal climate ##################
  
  tree_dir <- paste0("~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/", city, ".rds")
  tree_location <- read_rds(tree_dir) |>
    filter(!(lon == "NaN") & !(lat == "NaN")) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F)
  message(paste0("Daymet: Processing city: ", city))
  daymet <- daymet_raw |>
    summarise_Daymet_seasons(season_type)
  if (season_type != "full_season") {
    source(here("script/process_UHAT.R"), local = TRUE)
    process_UHAT_city(city, raw_uhat, season_type)
  }
  
  ###### UHAT temperature #######
  
  uhat_files <- list.files(
    path = paste0("~/lab-data/datasets/climate/UHAT/", switch(season_type,
                                                  full_season = "clean", early_season = "clean_early_season", later_season = "clean_later_season")),
    pattern = paste0("^", city, ".*tavg\\.tif$"),
    full.names = TRUE
  )
  uhat_files <- uhat_files[startsWith(basename(uhat_files), paste0(city, "_")) &
                             grepl("_(winter|summer)_[0-9]{4}_tavg\\.tif$", basename(uhat_files))]
  if (length(uhat_files) == 0L) stop("No winter/summer UHAT rasters found for: ", city)
  missing_files <- uhat_files[!file.exists(uhat_files)]
  if (length(missing_files) > 0) {
    stop("The following CHELSA TIF file(s) do not exist: {str_c(missing_files, collapse = ', ')}")
  }
  uhat_climate <- map(uhat_files, function(file) {
    r <- rast(file)
    dir_path <- ("~/lab-data/datasets/climate/UHAT/grid/")
    tifs <- list.files(
      dir_path,
      pattern = "^MYD11A1\\..*LST_Day_1km\\.tif$",
      full.names = TRUE
    )
    r_template <- terra::rast(tifs[1])
    crs(r) <- crs(r_template)  # ensure correct CRS
    tree_location <- tree_location|>
      st_transform(crs(r))
    
    vals <- terra::extract(r, tree_location, ID = FALSE)
    fname <- basename(file)
    fname <- str_remove(fname, "\\.tif$")
    season <- str_extract(fname, "winter|summer")
    year   <- str_extract(fname, "\\d{4}")
    colnames(vals) <- "value"
    vals <- bind_cols(
      tibble(id = tree_location$id, year = year, season = season),
      vals
    )
    vals}) |> 
    bind_rows() |>
    pivot_wider(
      names_from  = season,
      values_from = value,
      names_glue  = "UHAT_temp_{season}"
    )
  
  points_surface <- tree_location |>
    mutate(id = as.character(id)) |>
    left_join(daymet |>
                pivot_wider(names_from = season, values_from = c(Avg_srad, Sum_Sum_mm),
                            names_glue = "{.value}_{season}") |>
                rename_with(~ sub("^Avg_srad_", "Daymet_srad_", .x), starts_with("Avg_srad_")) |>
                rename_with(~ sub("^Sum_Sum_mm_", "Daymet_prcp_", .x), starts_with("Sum_Sum_mm_")),
              by = "id") |>
    left_join(uhat_climate |>
                mutate(year = as.double(year), id = as.character(id)), by = c("id", "year"))
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write_rds(points_surface, out_path, compress = "gz")
}


tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

if (sys.nframe() == 0L) {
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "6"))
  
  results <- parallel::mclapply(cities[23:76], function(city) {
    tryCatch({
      uhat_file <- paste0("~/lab-data/datasets/climate/UHAT/raw/", city, ".csv")
      if (!file.exists(uhat_file)) stop("No UHAT found")
      raw_uhat <- readr::read_csv(uhat_file, show_col_types = FALSE, num_threads = 1)
      
      daymet_files <- list.files(here("data/raw/CONUS/Daymet/all_city_trees", city),
                                 pattern = "^daymet_daily_[0-9]{4}_[0-9]{4}\\.csv\\.gz$", full.names = TRUE)
      if (length(daymet_files) != 1L) stop("Expected one Daymet file, found ", length(daymet_files))
      
      daymet_raw <- readr::read_csv(daymet_files[1],
                                    col_select = all_of(c("id", "date", "prcp", "srad")),
                                    col_types = cols(id = col_character(), date = col_character(),
                                                     prcp = col_double(), srad = col_double()),
                                    num_threads = 1, show_col_types = FALSE)
      
      for (season in c("early_season", "later_season"))
        prepare_tree_Grid_climate(city, raw_uhat, daymet_raw, fixed_season = TRUE, season_type = season)
      
      for (phe in c("SOS", "EOS")) {
        prepare_tree_Grid_climate(city, raw_uhat, daymet_raw,
                                  fixed_season = FALSE, preseason_length = 30, phenophase = phe)
        prepare_tree_Grid_climate(city, raw_uhat, daymet_raw,
                                  fixed_season = FALSE, preseason_length = 60, phenophase = phe)
      }
      
      message("Finished: ", city)
      data.frame(city, status = "success")
    }, error = function(e) {
      message("Failed: ", city, " - ", conditionMessage(e))
      data.frame(city, status = conditionMessage(e))
    })
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  
  print(dplyr::bind_rows(results))
}