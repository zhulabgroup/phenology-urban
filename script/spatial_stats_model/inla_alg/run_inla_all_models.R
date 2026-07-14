################################## package loading ##################################
library(zhulabtools)

load_packages(c("here", "rgbif", "ggplot2", "sf","sp",
                "glue", "tidyverse", "broom","lme4", "viridis",
                "broom.mixed", "coda", "patchwork", "terra", "stats",
                "foreach", "doParallel"))
select <- dplyr::select
## Rscript.r
options(inla.call = "inla")
library(INLA)

source(here("script/spatial_stats_model/prepare_model_input.R"))
source(here("script/spatial_stats_model/inla_alg/run_sp_inla.R"))

################################## parameter loading ##################################

# Read file paths from environment variables
# r_input_dir <- Sys.getenv("R_INPUT_DIR")
# r_output_dir <- Sys.getenv("R_OUTPUT_DIR")
# city <- Sys.getenv("CITY")
# 
r_input_dir <- "~/phenology-urban/data/"
r_output_dir <- "~/phenology-urban/data/proc/model_output/final_models_clean_outlier/"

tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

# r_output_dir <- paste0(r_output_dir, "inla_start/")
# r_output_dir <- paste0(r_output_dir, "inla_end/")
# r_output_dir_min <- paste0(r_output_dir, "inla_min/")
# r_output_dir_max <- paste0(r_output_dir, "inla_max/")
# r_output_dir <- paste0(r_output_dir, "inla_standardize/")
# r_output_dir <- paste0(r_output_dir, "inla_gsl/")
r_output_dir <- paste0(r_output_dir, "inla_sensitive/")

standardize <- FALSE

################################## Preparing standardization function ##################################
# remove_outliers_iqr <- function(df, value_col) {
#   df %>%
#     group_by(canonicalname) %>% 
#     mutate(
#       q1  = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
#       q3  = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
#       iqr = q3 - q1,
#       lo  = q1 - 1.5 * iqr,
#       hi  = q3 + 1.5 * iqr
#     ) %>%
#     ungroup() %>%
#     filter(
#       is.na(.data[[value_col]]) | dplyr::between(.data[[value_col]], lo, hi)
#     ) %>%
#     select(-q1, -q3, -iqr, -lo, -hi)
# }
# 
# 
# all_data_by_phe <- lapply(c("SOS", "EOS"), function(phe_type) {
#   d <- purrr::map_dfr(
#     cities,
#     ~ {
#       df <- prepare_data(.x,
#                          source = "Grid",
#                          taxa_col = "canonicalname",
#                          phe_type = phe_type)
#       df$city <- .x
#       df <- remove_outliers_iqr(df, value_col = "doy")
#       df
#     }
#   )
#   d$phe_type <- phe_type
#   d
# })
# names(all_data_by_phe) <- c("SOS", "EOS")
# 
# get_vars_to_std <- function(d) {
#   grep("^(temp_|prcp_|srad_)", names(d), value = TRUE)
# }
# std_params <- lapply(all_data_by_phe, function(d) {
#   vars_to_std <- get_vars_to_std(d)
#   mu <- sapply(vars_to_std, function(v) mean(d[[v]], na.rm = TRUE))
#   sd <- sapply(vars_to_std, function(v) stats::sd(d[[v]], na.rm = TRUE))
#   list(mu = mu, sd = sd)
# })
# write_rds(std_params, ("~/phenology-urban/data/proc/model_output/final_models_clean_outlier/inla_standardize/inla_standardization_params.rds"))

# std_params <- read_rds("~/phenology-urban/data/proc/model_output/final_models_clean_outlier/inla_standardize/inla_standardization_params.rds")
# standardize_with_params <- function(df, params, vars) {
#   for (v in vars) {
#     df[[v]] <- (df[[v]] - params$mu[[v]]) / params$sd[[v]]
#   }
#   df
# }

################################## run spatial model ##################################

##### Supplementary analysis for SOS/EOS
# for (city in cities) {
#   for (phe_type in c("GSL_range_season")) {
#     model_input_all <- prepare_data(city, source = "Grid_GSL", taxa_col= "canonicalname", phe_type = phe_type)
#     
#     if (standardize) {
#       model_input_all <- standardize_with_params(
#         model_input_all,
#         params = std_params[[phe_type]],
#         vars = get_vars_to_std(model_input_all)
#       )
#     }
#     run_inla_species(model_input_all, r_output_dir, phe_type, city)
#   }
# }

##### Sensitivity analysis for SOS/EOS threshold

for (city in cities[25:76]) {
  for (phe_type in c("SOS_04", "SOS_06", "EOS_04", "EOS_06")) {
    model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
    batch_status <- run_inla_species(model_input_all, r_output_dir, phe_type, city, n_cores = 10)
    print(batch_status)
  }
}

##### Formal analysis for SOS/EOS 

# for (phe_type in c("SOS", "EOS")) {
#   model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
#   run_inla_species(model_input_all, r_output_dir, phe_type, city)
#   # model_input_all <- prepare_data(city, source = "WU", taxa_col= "genus", phe_type = phe_type)
#   # run_inla_genus(model_input_all, r_output_dir, phe_type, city)
# }

################################## run non-spatial model ##################################
# source("script/spatial_stats_model/inla_alg/run_nonsp_inla.R")
# for (phe_type in c("SOS", "EOS")) {
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "canonicalname", phe_type = phe_type)
#   run_nonsp_inla(model_input_all, r_output_dir, phe_type, city, 
#                  taxon_col = "canonicalname", model_type = "no spatial")
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "genus", phe_type = phe_type)
#   run_nonsp_inla(model_input_all, r_output_dir, phe_type, city, 
#                  taxon_col = "genus", model_type = "no spatial")
# }
# 
# for (phe_type in c("SOS", "EOS")) {
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "canonicalname", phe_type = phe_type)
#   run_nonsp_inla(model_input_all, r_output_dir, phe_type, city, 
#                  taxon_col = "canonicalname", model_type = "random spatial")
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "genus", phe_type = phe_type)
#   run_nonsp_inla(model_input_all, r_output_dir, phe_type, city, 
#                  taxon_col = "genus", model_type = "random spatial")
# }
# 
# for (phe_type in c("SOS", "EOS")) {
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "canonicalname", phe_type = phe_type)
#   run_nonsp_inla(model_input_all, r_output_dir, phe_type, city, 
#                  taxon_col = "canonicalname", model_type = "fix spatial")
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "genus", phe_type = phe_type)
#   run_nonsp_inla(model_input_all, r_output_dir, phe_type, city, 
#                  taxon_col = "genus", model_type = "fix spatial")
# }

