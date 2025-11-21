################################## package loading ##################################
library(zhulabtools)

load_packages(c("here", "rgbif", "dplyr","lubridate", "ggplot2", "sf","sp", 
                "readr", "purrr", "stringr","glue", "tidyr", "broom","lme4", 
                "broom.mixed", "coda", "patchwork", "terra", "stats", "viridis"))
select <- dplyr::select
## Rscript.r
options(inla.call = "inla")
library(INLA)
################################## parameter loading ##################################

# Read file paths from environment variables
r_input_dir <- Sys.getenv("R_INPUT_DIR")
r_output_dir <- Sys.getenv("R_OUTPUT_DIR")
city <- Sys.getenv("CITY")

# r_input_dir <- "~/phenology-urban/data/"
# r_output_dir <- "~/phenology-urban/data/proc/model_output/final_models_clean_outlier/"
# city = "Albuquerque"

tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

source(here("script/spatial_stats_model/prepare_model_input.R"))


################################## run spatial model ##################################
source(here("script/spatial_stats_model/inla_alg/run_sp_inla.R"))

r_output_dir <- paste0(r_output_dir, "inla_start/")
r_output_dir <- paste0(r_output_dir, "inla_end/")

for (city in cities) {
  for (phe_type in c("SOS", "EOS")) {
    model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
    run_inla_species(model_input_all, r_output_dir, phe_type, city, fitting_type = "end")
  }
}

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

