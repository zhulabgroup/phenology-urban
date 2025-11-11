install.packages('BatchPlanet', repos = c('https://yiluansong.r-universe.dev', 'https://cloud.r-project.org'))
library(BatchPlanet)

tree_dir <- "~/lab-data/datasets/vegetation/PS/CONUS_urban/Tree_point"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
tree_df <- rds_files |>
  map_dfr(read_rds) |>
  mutate(taxa = if_else(genus == "Validation", "Validation", gsub(" ", "", taxa))) 

dir_data <- "~/lab-data/datasets/vegetation/PS/CONUS_urban/"

site_list <- tree_df |>
  distinct(site) |>
  pull()

for (site in c("StLouis", "StGeorge")){
  df_thres <- set_thresholds(thres_up = 0.5, thres_down = 0.5)
  calculate_phenological_metrics_batch(dir = dir_data, v_year = 2017:2023, v_site = site,
                                        df_thres = df_thres, var_index = "evi",
                                       min_days = 80, check_seasonality = T, extend_to_previous_year = 275,
                                       extend_to_next_year = 90, num_cores = 1)
}


# dir <- "~/lab-data/datasets/vegetation/PS/CONUS_urban/doy"
# files <- list.files(dir, pattern = "St\\.Louis", full.names = T)
# new_files <- sub("St\\.Louis", "StLouis", files)
# file.rename(files, new_files)
