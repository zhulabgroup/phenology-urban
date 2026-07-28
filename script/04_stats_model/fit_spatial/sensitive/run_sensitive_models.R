################################## package loading ##################################
library(zhulabtools)

load_packages(c("here", "tidyverse", "sf","sp", "coda", "terra", "nimble", "spmodel", "spBayes", "foreach", "doParallel", "RANN"))

################################## parameter loading ##################################
source(here("script/04_stats_model/fit_spatial/prepare_model_input.R"))
source(here("script/04_stats_model/fit_spatial/sensitive/alg_sp_freq.R"))
source(here("script/04_stats_model/fit_spatial/sensitive/alg_spbayes_mat.R"))

r_input_dir <- here("data/")
r_output_dir <- here("data/proc/model_output/final_models_clean_outlier/")

tree_dir <- here("data/raw/CONUS/StreetTrees/cleaned_ps_input/")
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

################################## run freq model ##################################

r_output_dir_freq <- paste0(r_output_dir, "freq/")
if (!dir.exists(r_output_dir_freq)) {
  dir.create(r_output_dir_freq, recursive = TRUE)
}

for (city in cities) {
  for (phe_type in c("SOS", "EOS")) {
    model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
    run_freq_species(model_input_all, r_output_dir, phe_type, city)
  }
}
################################## run mcmc model ##################################

r_output_dir_mcmc <- paste0(r_output_dir, "mcmc/")
if (!dir.exists(r_output_dir_mcmc)) {
  dir.create(r_output_dir_mcmc, recursive = TRUE)
}

for (city in cities) {
  for (phe_type in c("SOS", "EOS")) {
    model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
    run_mcmc_spatial_mat(model_input_all, r_output_dir_mcmc, phe_type, city, taxon_col = "canonicalname")
  }
}