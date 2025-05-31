extract_chelsa <- function(chelsa_ras, occ_data, lon_col = "longitude", lat_col = "latitude") {
  if (!all(c(lon_col, lat_col) %in% colnames(occ_data))) {
    cli_alert_warning("The data lacks required longitude/latitude columns.")
    return(NULL)
  }
  
  chelsa_values <- terra::extract(chelsa_ras, occ_data[, c(lon_col, lat_col)], ID = FALSE)
  colnames(chelsa_values) <- names(chelsa_ras)
  # Manually apply scale/offset per variable
  combined_data <- bind_cols(occ_data, chelsa_values) |>
    mutate(bio1 = bio1 * 0.1 -273.15,
           bio12 = bio12 * 0.1)
  return(combined_data)
}

load_chelsa <- function(chelsa_vars = c("bio1"), chelsa_dir) {
  chelsa_files <- chelsa_vars |>
    map(~ file.path(chelsa_dir, str_c("CHELSA_", .x, "_1981-2010_V.2.1.tif"))) |>
    unlist()
  
  missing_files <- chelsa_files[!file.exists(chelsa_files)]
  if (length(missing_files) > 0) {
    stop("The following CHELSA TIF file(s) do not exist: {str_c(missing_files, collapse = ', ')}")
  }
  
  chelsa_ras <- terra::rast(chelsa_files, raw = TRUE)
  names(chelsa_ras) <- chelsa_vars
  return(chelsa_ras)
}