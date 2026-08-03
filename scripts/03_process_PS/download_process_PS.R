if (!require("BatchPlanet")){
  install.packages('BatchPlanet', repos = c('https://yiluansong.r-universe.dev', 'https://cloud.r-project.org'))
}
pacman::p_load(BatchPlanet)

tree_dir <- here("data/input/tree_list")
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
tree_df <- rds_files |> #25
  map_dfr(read_rds) |>
  mutate(taxa = if_else(genus == "Validation", "Validation", gsub(" ", "", taxa))) 

dir_data <- here("data/input/phenology")
setting <- set_planetscope_parameters(
  api_key = set_api_key(),
  item_name = "PSScene",
  asset = "ortho_analytic_4b_sr",
  product_bundle = "analytic_sr_udm2",
  cloud_lim = 1,
  harmonized = T
)

order_planetscope_imagery_batch(dir = dir_data, df_coordinates = tree_df, v_year = 2016:2024, setting = setting)
site_list <- tree_df |>
  distinct(site) |>
  pull()

# n_split = 4
# year_splits <- split(2016:2024, cut(seq_along(2016:2024), n_split, labels = FALSE))
for (site in c(site_list)){
  tree_df_site <- tree_df |> filter(site == !!site)
  big_genera <- names(which(table(tree_df_site$genus) > 2000))
  tree_df_site_group <- tree_df_site |>
    mutate(group = if_else(genus %in% big_genera & !is.na(taxa) & nzchar(taxa),
                           taxa,
                           genus
    )) |>
    dplyr::select(site, id, lon, lat, group)
  update_group <- unique(tree_df_site_group$group)

  message("Downloading images...")
  download_planetscope_imagery_batch(dir = dir_data,  v_site = site, v_year = 2016:2024, setting = setting, num_cores = 30, overwrite = F)
  
  message("Retrieving time series...")
  retrieve_planetscope_time_series_batch(dir = dir_data, df_coordinates = tree_df_site_group, v_site = site, v_group = update_group, num_cores = 30)
  
  message("Cleaning time series...")
  clean_planetscope_time_series_batch(dir = dir_data, v_site = site, num_cores = 12, calculate_evi = T) #

  message("Processing time series...")
  df_thres <- set_thresholds(thres_up = c(0.2, 0.4, 0.5, 0.6, 0.8), thres_down = c(0.2, 0.4, 0.5, 0.6, 0.8))
  calculate_season_metrics_batch(dir = dir_data, v_year = 2017:2023, v_site = site, df_thres = df_thres, var_index = "evi",
                                       min_days = 1, check_seasonality = T, extend_to_previous_year = 275,
                                       extend_to_next_year = 90, num_cores = 15)
}