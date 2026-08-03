run_freq_species <- function(model_input_all, r_output_dir, phe_type, city, fitting_type = "doy"){
  species_list <- unique(model_input_all$canonicalname)
  failed_sp <- c()
  for (species in species_list) {
    tryCatch({
      folder_path <- paste0(r_output_dir, "/",city,"/canonicalname/")
      if (!dir.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE)
      }
      species_clean <- gsub(" ", "_", species)
      model_path <- paste0(folder_path, species_clean, "_", phe_type, ".rds")
      
      cat(sprintf("\n Running freq model for species: %s, phe_type: %s, city: %s\n", 
                  species, phe_type, city))
      
      phe_metric <- "doy"
      model_input_sf <- model_input_all %>%
        filter(canonicalname == species)
      
      x <- model_input_sf[[phe_metric]]
      q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
      iqr <- diff(q)
      lo  <- q[1] - 1.5 * iqr
      hi  <- q[2] + 1.5 * iqr
      model_input_sf <- model_input_sf %>%
        filter(is.na(.data[[phe_metric]]) |
                 dplyr::between(.data[[phe_metric]], lo, hi))
      model_input_sf <- model_input_sf %>%
        st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
        st_transform(crs = 5070) %>%
        mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2],
               year = as.factor(year)) 
      model_input <- model_input_sf %>%
        st_drop_geometry()
      if (nrow(model_input) < 30) {
        message(paste("Not enough unique coordinates for taxon:", species, "with", nrow(model_input), "unique points. Skipping..."))
        next  # Skip the current iteration of the loop
      }
      if (file.exists(model_path)) {
        cat(sprintf("\n File already exists for species: %s, phe_type: %s, city: %s. Skipping freq run.\n", 
                    species, phe_type, city))
        next  # Skip the current iteration of the loop
      } # if the file exists, skip to the next iteration
      
      ######################################## Run/Save the model ########################################
      ### Run the model
      phe_metric <- "doy"
      temp_col <- model_input %>% select(starts_with("temp")) %>% names()
      prcp_col <- model_input %>% select(starts_with("prcp")) %>% names()
      srad_col <- model_input %>% select(starts_with("srad")) %>% names()
      
      form <- reformulate(
        termlabels = c(temp_col, prcp_col, srad_col),
        response = phe_metric
      )
      
      res <- splm(
        form,
        model_input_sf,
        spcov_type = "matern",
        random = ~ year
      )
      # res$waic$waic
      # summary(res)
      
      fix_effect <- as.data.frame(tidy(res))%>%
        mutate(
          Lower = estimate + qnorm(0.025) * std.error,
          Upper = estimate + qnorm(0.975) * std.error
        ) %>%
        rename(Estimate = estimate, Factor = term) %>%
        mutate(Sig = ifelse(p.value<0.05, "*", "")) %>%
        filter(stringr::str_starts(Factor, "temp")) %>%
        mutate(Factor = ifelse(stringr::str_starts(Factor, "temp"), species, Factor))
      
      list(res = res,
           fix_effect = fix_effect,
           species = species) %>%
        write_rds(model_path, compress = "gz")
  }, error = function(e){
    message("Something wrong", species)
    failed_sp <- c(failed_sp, species)
  })
  }
message("Failed species, ", failed_sp)
}