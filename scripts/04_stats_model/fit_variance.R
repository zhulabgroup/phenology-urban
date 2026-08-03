pacman::p_load(tidytext, lme4)

fit_variance_partition <- function(df, species_name) {
  d <- df |> 
    filter(canonicalname == species_name) |>
    mutate(tree_uid = interaction(city, id, drop = TRUE))
  
  if (n_distinct(d$city) < 2 || n_distinct(d$tree_uid) < 10) {
    return(tibble(
      species = species_name,
      sigma_city = NA_real_,
      sigma_within = NA_real_,
      sigma_residual = NA_real_,
      prop_city = NA_real_,
      prop_within = NA_real_
    ))
  }
  
  model <- tryCatch(
    lmer(doy ~ 1 +  (1 | year) + (1 | city) + (1 | city:id), data = d),
    error = function(e) NULL
  )
  
  if (is.null(model)) {
    return(tibble(
      species = species_name,
      sigma_city = NA_real_,
      sigma_within = NA_real_,
      sigma_residual = NA_real_,
      prop_city = NA_real_,
      prop_within = NA_real_,
      compared_ratio = NA_real_,
      phe_type = phe_type,
      sample_city = n_distinct(d$city)
    ))
  }
  
  vc <- as.data.frame(VarCorr(model))
  sigma_city     <- vc$vcov[vc$grp == "city"]
  sigma_within  <- vc$vcov[vc$grp == "city:id"]
  sigma_residual <- vc$vcov[vc$grp == "Residual"]
  sigma_year <- vc$vcov[vc$grp == "year"]
  total_var <- sigma_city + sigma_within + sigma_residual + sigma_year
  tibble(
    species = species_name,
    sigma_city = sigma_city,
    sigma_within = sigma_within,
    sigma_residual = sigma_residual,
    prop_city = sigma_city / total_var,
    prop_within = sigma_within / total_var,
    prop_residual = sigma_residual / total_var,
    prop_city2 = prop_city / (sum(prop_city, prop_within)),
    prop_within2 = prop_within / (sum(prop_city, prop_within)),
    phe_type = phe_type,
    sample_city = n_distinct(d$city)
  )
}

for (phe_type in c("SOS", "EOS")) {
  all_city_input <- read_rds(here::here(paste0("release/input/phenology/no_outlier/nooutlier_", phe_type, "_input.rds"))) |>
    filter(year<2024) |>
    mutate(year = as.factor(year))
  species_list <- unique(all_city_input$canonicalname)
  variance_partition_all <- map_dfr(
    species_list,
    ~ fit_variance_partition(all_city_input, .x)
    )
  write_rds(variance_partition_all, here::here(paste0("release/intermediate/variance_decompose/variance_partition_", phe_type, ".rds")))
}
