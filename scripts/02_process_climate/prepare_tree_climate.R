pacman::p_load(terra)

prepare_tree_Grid_climate <- function(city, fixed_season = TRUE, season_type = "full_season") {
  if (!is.logical(fixed_season) || length(fixed_season) != 1L || is.na(fixed_season)) {
    stop("fixed_season must be TRUE or FALSE")
  }
  if (!fixed_season) stop("fixed_season = FALSE is not implemented")
  season_type <- match.arg(season_type, c("full_season", "early_season", "later_season"))
  output_dir <- switch(season_type, full_season = "Grid", early_season = "Grid_early_season",
                       later_season = "Grid_later_season")
  out_path <- here("data/input/tree_climate", output_dir, paste0(city, ".rds"))
  
  ################## Calculate seasonal climate ##################
  message(paste0("Processing example city: ", city, "..."))
  
  tree_dir <- here(paste0("data/input/tree_list/", city, ".rds"))
  tree_location <- read_rds(tree_dir) |>
    filter(!(lon == "NaN") & !(lat == "NaN")) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F)
  
  ###### Daymet shortwave radiation #######
  daymet_srad <- read_csv(paste0(here("data/input/climate/Daymet/"), city, "/srad_monthly_2016_2024.csv"), show_col_types = FALSE)
  daymet_prcp <- read_csv(paste0(here("data/input/climate/Daymet/"), city, "/prcp_monthly_2016_2024.csv"), show_col_types = FALSE)
  if (season_type != "full_season") {
    source(here("scripts/02_process_climate/process_UHAT.R"), local = TRUE)
    process_UHAT_city(city, season_type)
  }
    
  ###### UHAT temperature #######
    
  uhat_files <- list.files(
    path = here("data/input/climate/UHAT", switch(season_type,
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
    dir_path <- here("data/input/climate/UHAT/grid/")
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
  
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    write_rds(points_surface, out_path, compress = "gz")
}


tree_dir <- here("data/input/tree_list")
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))
if (sys.nframe() == 0L) for (city in cities){
  tryCatch(
    {
      prepare_tree_Grid_climate(city)
    },
    error = function(e) {
      message("Skipping ", city, ": ", e$message)
    }
  )}

