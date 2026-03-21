################################## package loading ##################################
library(zhulabtools)

load_packages(c("here", "tidyverse", "sf","sp", "spmodel"))
select <- dplyr::select

################################## parameter loading ##################################

r_input_dir <- "~/phenology-urban/data/"
r_output_dir <- "~/phenology-urban/data/proc/model_output/final_models_clean_outlier/"

tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

source(here("script/spatial_stats_model/prepare_model_input.R"))

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
source(here("script/spatial_stats_model/freq_alg/run_sp_freq.R"))

r_output_dir <- paste0(r_output_dir, "freq/")
if (!dir.exists(r_output_dir)) {
  dir.create(r_output_dir, recursive = TRUE)
}
standardize <- FALSE

for (city in cities) {
  for (phe_type in c("SOS", "EOS")) {
    model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
    # if (standardize) {
    #   model_input_all <- standardize_with_params(
    #     model_input_all,
    #     params = std_params[[phe_type]],
    #     vars = get_vars_to_std(model_input_all)
    #   )
    # }
    run_freq_species(model_input_all, r_output_dir, phe_type, city)
  }
}

