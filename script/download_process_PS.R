# remotes::install_local("~/batchplanet", dependencies = TRUE, build_vignettes = FALSE, force = TRUE)

library(batchplanet)

# matching_npn <- readRDS("C:/Users/feixi/Desktop/PlanetScope/City_tree/Tree_point/summary/matching_npn.rds")
# npn_points <- readRDS("C:/Users/feixi/Desktop/PlanetScope/CONUS/npn_points.rds") |>
#   left_join(matching_npn, by = c("site" = "site.x")) |>
#   filter(!is.na(city)) |>
#   mutate(city = gsub("\\.", "_", city))

## Add the validation point to consistent city_list

# npn_points_updated <- npn_points |>
#   rowwise() |>
#   mutate(
#     new_id = {
#       file_path <- file.path("C:/Users/feixi/Desktop/PlanetScope/City_tree/Tree_point",
#                              paste0(city, ".rds"))
#       
#       if (file.exists(file_path)) {
#         city_df <- readRDS(file_path)
#       } else {
#         city_df <- tibble(id = integer(),
#                           lon = double(),
#                           lat = double(),
#                           genus = character())
#         message(paste0("Not found"))
#       }
#       match_idx <- which(abs(city_df$lon - lon) < 1e-5 & abs(city_df$lat - lat) < 1e-5)
#       
#       if (length(match_idx) > 0) {
#         matched_id <- city_df$id[match_idx[1]]
#         message("Existed")
#       } else {
#         new_id <- ifelse(nrow(city_df) == 0, 1, max(city_df$id, na.rm = TRUE) + 1)
#         new_row <- tibble(id = new_id,
#                           site = city,
#                           lon = lon,
#                           lat = lat,
#                           genus = "Validation")
#         city_df <- bind_rows(city_df, new_row)
#         write_rds(city_df, file_path)
#         message("Add new row")
#         matched_id <- new_id
#       }
#       matched_id
#     }
#   ) |>
#   ungroup()
# write_rds(npn_points_updated, "C:/Users/feixi/Desktop/PlanetScope/City_tree/Tree_point/summary/matching_npn.rds")


tree_dir <- "~/lab-data/datasets/vegetation/PS/CONUS_urban/Tree_point/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
tree_df <- rds_files[1] |>
  map_dfr(read_rds) |>
  dplyr::select(site, id, lon, lat, group = genus)

dir_data <- "~/lab-data/datasets/vegetation/PS/CONUS_urban/"
setting <- set_planetscope_parameters(
  api_key = set_api_key(),
  item_name = "PSScene",
  asset = "ortho_analytic_4b_sr",
  product_bundle = "analytic_sr_udm2",
  cloud_lim = 1,
  harmonized = T
)

# order_planetscope_imagery_batch(dir = dir_data, df_coordinates = tree_df, v_year = 2016:2024, setting = setting)
# Then wait until all orders success

site_list <- tree_df |>
  distinct(site) |>
  pull()

for (site in site_list){
  message(paste0("Downloading data for city:", site))
  download_planetscope_imagery_batch(dir = dir_data,  v_site = site, v_year = 2016:2024, setting = setting, num_cores = 12, overwrite = F)
  retrieve_planetscope_time_series_batch(dir = dir_data, df_coordinates = tree_df, v_site = site, num_cores = 12)
  
  # message("Retrieval done. Cleaning raw imagery...")
  # raw_dir <- file.path(dir_data, "raw", site)
  # items <- list.files(raw_dir, full.names = TRUE)
  # items_to_delete <- items[basename(items) != "orders"]
  # unlink(items_to_delete, recursive = TRUE, force = TRUE)
  
  message("Cleaning done. Processing time series...")
  clean_planetscope_time_series_batch(dir = dir_data, v_site = site, num_cores = 12, calculate_evi = T)
  df_thres <- set_thresholds(thres_up = 0.5, thres_down = 0.5)
  calculate_phenological_metrics_batch(dir = dir_data, v_site = site, df_thres = df_thres, var_index = "evi", 
                                       min_days = 80, check_seasonality = T, extend_to_previous_year = 275, 
                                       extend_to_next_year = 90, num_cores = 12)
  # check_seasonality?
}
# 
# 
# clean_28117_allGroup <- clean_28117_allGroup |>
#   filter(year!=2016)
# saveRDS(clean_28117_allGroup,"C:/Users/feixi/Desktop/PS_NPN/CONUS/clean/clean_28117_allGroup.rds" )

