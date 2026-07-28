####################### INLA model results tidy #######################
tree_dir <- here("data/raw/CONUS/StreetTrees/cleaned_ps_input/")
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))
phe_types <- c("SOS", "EOS")

prepare_coef_summary <- function(r_output_dir, phe_type){
  pattern = paste0(phe_type,"\\.rds$")
  files <- list.files(r_output_dir, pattern = pattern, full.names = TRUE)
  coef_summary <- data.frame()
  # for one genus
  for (file in files) {
    message(paste0("Processing file: ", file))
    model_data <- readRDS(file)
    species <- model_data$species
    # beta_names <- c("beta[0]", "beta[AvgTemp_mean]", "beta[Sum_mm_sum]", "beta[srad_mean]")
    # df_list <- lapply(seq_along(model_data$res$marginals.fixed), function(i) {
    #   data.frame(x = model_data$res$marginals.fixed[[i]][,1], 
    #              density = model_data$res$marginals.fixed[[i]][,2], 
    #              beta = beta_names[i])
    # })
    # df <- bind_rows(df_list)
    beta_avg_temp <- model_data$res$marginals.fixed[[2]]
    
    coef_summary <- coef_summary %>%
      rbind(model_data$fix_effect)
  }
  
  coef_summary <- coef_summary %>%
    mutate(SigAlpha = ifelse(Sig == "*", 1, 0.2)) %>%
    mutate(Color = ifelse(Estimate >0, "#0D0887FF" , "#FCA636FF"))
  write_rds(coef_summary, paste0(r_output_dir, phe_type, "_summary.rds"))
}

for (city in cities) { 
  message(paste0("Processing city: ", city))
  r_output_dir <- paste0(here("data/proc/model_output/final_models_clean_outlier/spatial_spde/"), city, "/canonicalname/")
  for (phe_type in phe_types){
    prepare_coef_summary(r_output_dir, phe_type)
  }
}

############################### MCMC model result #########################

prepare_coef_summary_mcmc <- function(r_output_dir, phe_type){
  files <- list.files(r_output_dir, pattern = paste0(phe_type, "\\.rds$"), full.names = TRUE)
  combined_param_data <- purrr::map(files, function(file) {
    model_data <- readRDS(file)
    genus <- sub(paste0("(.*?)_", phe_type, "\\.rds"), "\\1", basename(file))
    ggmcmc::ggs(model_data$samples) |>
      dplyr::filter(stringr::str_starts(Parameter, stringr::fixed("temp"))) |>
      dplyr::rename(ParameterOriginal = Parameter) |>
      dplyr::mutate(
        Parameter = genus,
        converge  = model_data$converge
        # waic = model_data$WAIC$WAIC,
        # milk = model_data$mean_logml
      )
  }) |>
    bind_rows()
  
  # Summarize
  summary_data <- combined_param_data |>
    group_by(Parameter) |>
    summarise(
      Estimate = median(value, na.rm = TRUE),
      Lower = quantile(value, 0.025, na.rm = TRUE),
      Upper = quantile(value, 0.975, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      # waic = first(waic),
      converge = first(converge),
      # milk = first(milk),
      .groups = "drop"
    ) |>
    mutate(Factor= Parameter,
      Sig = ifelse(Upper * Lower > 0, "*", ""),
           SigAlpha = ifelse(Sig == "*", 1, 0.2),
           Color = ifelse(Estimate >0, "#0D0887FF" , "#FCA636FF")) 
  write_rds(summary_data, paste0(r_output_dir, phe_type, "_summary.rds"))
}
  
for (city in cities) { 
  message(paste0("Processing city: ", city))
  r_output_dir <- paste0(here("data/proc/model_output/final_models_clean_outlier/mcmc/spatial_mat/"), city, "/canonicalname/")
  for (phe_type in phe_types){
    prepare_coef_summary_mcmc(r_output_dir, phe_type)
  }
}
