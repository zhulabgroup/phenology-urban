################################## package loading ##################################
pacman::p_load(rgbif, foreach, doParallel, terra)

select <- dplyr::select
## Rscript.r
options(inla.call = "inla")
library(INLA)

source(here("scripts/04_stats_model/fit_spatial/prepare_model_input.R"))
source(here("scripts/04_stats_model/fit_spatial/inla/alg_sp_inla.R"))

################################## parameter loading ##################################

r_output_dir <- here("data/intermediate/spatial/")

tree_dir <- here("data/raw/CONUS/StreetTrees/cleaned_ps_input/")
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

r_output_dir <- paste0(r_output_dir, "inla/")

##### Formal analysis for SOS/EOS 

for (city in cities) {
  for (phe_type in c("SOS", "EOS")) {
    model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
    batch_status <- run_inla_species(model_input_all, r_output_dir, phe_type, city, n_cores = 1)
    print(batch_status)
  }
}



