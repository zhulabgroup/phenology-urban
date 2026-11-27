################################## package loading ##################################
pacman::p_load(coda, terra, nimble, spmodel, spBayes, foreach, doParallel, RANN)

################################## parameter loading ##################################
source(here("scripts/04_stats_model/fit_spatial/prepare_model_input.R"))
source(here("scripts/04_stats_model/fit_spatial/sensitive/alg_sp_freq.R"))
source(here("scripts/04_stats_model/fit_spatial/sensitive/alg_spbayes_mat.R"))
source(here("scripts/04_stats_model/fit_spatial/inla/alg_sp_inla.R"))
options(inla.call = "inla")
library(INLA)


output_root <- here("data/intermediate/spatial/")
tree_dir <- here("data/input/tree_list/")
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

################################## run different preseason model ##################################

source_output_dirs <- c(
  Grid_early_season = "inla_early_season",
  Grid_later_season = "inla_later_season",
  Grid_pre60 = "inla_pre60",
  Grid_pre30 = "inla_pre30"
)

for (source_name in names(source_output_dirs)) {
  r_output_dir <- paste0(file.path(output_root, source_output_dirs[[source_name]]), "/")
  for (city in cities) {
    for (phe_type in c("SOS", "EOS")) {
      model_input_all <- prepare_data(city, source = source_name, taxa_col = "canonicalname", phe_type = phe_type)
      batch_status <- run_inla_species(model_input_all, r_output_dir, phe_type, city, n_cores = 1)
      print(batch_status)
    }
  }
}

################################## run freq model ##################################

r_output_dir_freq <- paste0(output_root, "freq/")
if (!dir.exists(r_output_dir_freq)) {
  dir.create(r_output_dir_freq, recursive = TRUE)
}

for (city in cities) {
  for (phe_type in c("SOS", "EOS")) {
    model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
    run_freq_species(model_input_all, r_output_dir_freq, phe_type, city)
  }
}
################################## run mcmc model ##################################

r_output_dir_mcmc <- paste0(output_root, "mcmc/")
if (!dir.exists(r_output_dir_mcmc)) {
  dir.create(r_output_dir_mcmc, recursive = TRUE)
}

for (city in cities) {
  for (phe_type in c("SOS", "EOS")) {
    model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
    run_mcmc_spatial_mat(model_input_all, r_output_dir_mcmc, phe_type, city, taxon_col = "canonicalname")
  }
}

################################## run GSL model ##################################

r_output_dir_gsl <- paste0(output_root, "inla_gsl/")
if (!dir.exists(r_output_dir_gsl)) {
  dir.create(r_output_dir_gsl, recursive = TRUE)
}

for (city in cities) {
  for (phe_type in c("GSL_winter", "GSL_summer")) {
    model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
    run_mcmc_spatial_mat(model_input_all, r_output_dir_gsl, phe_type, city, taxon_col = "canonicalname")
  }
}