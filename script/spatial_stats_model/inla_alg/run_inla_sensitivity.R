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

################################## parameter loading ##################################
source(here("script/spatial_stats_model/prepare_model_input.R"))
source(here("script/spatial_stats_model/inla_alg/run_sp_inla.R"))

output_root <- here("data/proc/model_output/final_models_clean_outlier")
tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
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
      batch_status <- run_inla_species(model_input_all, r_output_dir, phe_type, city, n_cores = 15)
      print(batch_status)
    }
  }
}