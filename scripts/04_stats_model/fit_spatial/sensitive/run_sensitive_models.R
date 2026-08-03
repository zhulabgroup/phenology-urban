################################## package loading ##################################
pacman::p_load(coda, terra, nimble, spmodel, spBayes, foreach, doParallel, RANN)

################################## parameter loading ##################################
source(here("scripts/04_stats_model/fit_spatial/prepare_model_input.R"))
source(here("scripts/04_stats_model/fit_spatial/sensitive/alg_sp_freq.R"))
source(here("scripts/04_stats_model/fit_spatial/sensitive/alg_spbayes_mat.R"))

r_output_dir <- here("release/intermediate/spatial/")
cities <- names(city_label)

################################## run freq model ##################################

r_output_dir_freq <- paste0(r_output_dir, "freq/")
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

################################## run GSL model ##################################

r_output_dir_gsl <- paste0(r_output_dir, "inla_gsl/")
if (!dir.exists(r_output_dir_gsl)) {
  dir.create(r_output_dir_gsl, recursive = TRUE)
}

for (city in cities) {
  for (phe_type in c("GSL_winter", "GSL_summer")) {
    model_input_all <- prepare_data(city, source = "Grid", taxa_col= "canonicalname", phe_type = phe_type)
    run_mcmc_spatial_mat(model_input_all, r_output_dir_gsl, phe_type, city, taxon_col = "canonicalname")
  }
}