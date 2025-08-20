################################## package loading ##################################
library(zhulabtools)

load_packages(c("here", "rgbif", "dplyr","lubridate", "nimble", "sf","sp", "taxize","bridgesampling",
                "tidyverse", "coda", "terra", "nimble", "spBayes", "fields",
                "ggplot2", "RANN", "foreach", "doParallel"))


# Read file paths from environment variables
r_input_dir <- Sys.getenv("R_INPUT_DIR")
r_output_dir <- Sys.getenv("R_OUTPUT_DIR")
city <- Sys.getenv("CITY")

# r_input_dir <- "~/phenology-urban/data/"
# r_output_dir <- "~/phenology-urban/data/proc/model_output/final_models/"
# city = "NY"
# phe_type = "SOS"

select <- dplyr::select
r_output_dir <- paste0(r_output_dir, "mcmc/")
source(here("script/spatial_stats_model/prepare_data_genus_species.R"))

################################## run spatial model ##################################
# source("script/spatial_stats_model/mcmc_alg/run_spbayes_exp.R")
# for (phe_type in c("SOS", "EOS")) {
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "canonicalname", phe_type = phe_type)
#   run_mcmc_spatial_exp(model_input_all, r_output_dir, phe_type, city, taxon_col = "canonicalname")
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "genus", phe_type = phe_type)
#   run_mcmc_spatial_exp(model_input_all, r_output_dir, phe_type, city, taxon_col = "genus")
# }


source("script/spatial_stats_model/mcmc_alg/run_spbayes_mat.R")
for (phe_type in c("SOS", "EOS")) {
  model_input_all <- prepare_data(city, source = "WU", taxa_col= "genus", phe_type = phe_type)
  run_mcmc_spatial_mat(model_input_all, r_output_dir, phe_type, city, taxon_col = "genus")
  model_input_all <- prepare_data(city, source = "WU", taxa_col= "canonicalname", phe_type = phe_type)
  run_mcmc_spatial_mat(model_input_all, r_output_dir, phe_type, city, taxon_col = "canonicalname")
}

# 
# ################################## run non-spatial model ##################################
# 
# source("script/spatial_stats_model/mcmc_alg/run_bayes_ols.R")
# for (phe_type in c("SOS", "EOS")) {
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "canonicalname", phe_type = phe_type)
#   run_mcmc_nosp(model_input_all, r_output_dir, phe_type, city, taxon_col = "canonicalname")
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "genus", phe_type = phe_type)
#   run_mcmc_nosp(model_input_all, r_output_dir, phe_type, city, taxon_col = "genus")
# }
# 
# ################################## run random-group model ##################################
# source("script/spatial_stats_model/mcmc_alg/run_bayes_randomsp.R")
# for (phe_type in c("SOS", "EOS")) {
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "canonicalname", phe_type = phe_type)
#   run_mcmc_randomsp(model_input_all, r_output_dir, phe_type, city, taxon_col = "canonicalname")
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "genus", phe_type = phe_type)
#   run_mcmc_randomsp(model_input_all, r_output_dir, phe_type, city, taxon_col = "genus")
# }
# 
# ################################## run fixed-group model ##################################
# source("script/spatial_stats_model/mcmc_alg/run_bayes_fixedsp.R")
# for (phe_type in c("SOS", "EOS")) {
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "canonicalname", phe_type = phe_type)
#   run_mcmc_fixedsp(model_input_all, r_output_dir, phe_type, city, taxon_col = "canonicalname")
#   model_input_all <- prepare_data(city, source = "WU", taxa_col= "genus", phe_type = phe_type)
#   run_mcmc_fixedsp(model_input_all, r_output_dir, phe_type, city, taxon_col = "genus")
# }
