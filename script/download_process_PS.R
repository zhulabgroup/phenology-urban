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
tree_df <- rds_files[67] |>
  map_dfr(read_rds) |>
  mutate(taxa = if_else(genus == "Validation", "Validation", gsub(" ", "", taxa))) 

dir_data <- "/scratch/zhukai_root/zhukai0/jializ/PS/CONUS_urban/"
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

n_split = 4
year_splits <- split(2016:2024, cut(seq_along(2016:2024), n_split, labels = FALSE))
for (site in site_list){
  tree_df_site <- tree_df |> filter(site == !!site)
  big_genera <- names(which(table(tree_df_site$genus) > 2000))
  tree_df2 <- tree_df_site |>
    mutate(group = if_else(genus %in% big_genera & !is.na(taxa) & nzchar(taxa),
                           taxa,
                           genus
    )) |>
    dplyr::select(site, id, lon, lat, group)
  update_group <- unique(tree_df2$group)
  
  message("Downloading time series...")
  for (j in 4:length(year_splits)){
    message(paste0("Downloading Part ", j," data for city: ", site))
    download_planetscope_imagery_batch(dir = dir_data,  v_site = site, v_year = year_splits[[j]], setting = setting, num_cores = 30, overwrite = F)
    retrieve_planetscope_time_series_batch(dir = dir_data, df_coordinates = tree_df2, v_site = site, num_cores = 30)

    message("Retrieval done. Cleaning Part ", j," raw imagery...")
    raw_dir <- file.path(dir_data, "raw", site)
    items <- list.files(raw_dir, full.names = TRUE)
    items_to_delete <- items[basename(items) != "orders"]
    unlink(items_to_delete, recursive = TRUE, force = TRUE)

    ts_dir <- paste0(dir_data, "ts/")
    files <- list.files(ts_dir, pattern = paste0("^ts_", site, "_.*\\.rds$"), full.names = TRUE)
    new_files <- file.path(dirname(files),paste0("part_", j, "_", basename(files)))
    file.rename(files, new_files)
  }
  files <- list.files(ts_dir, pattern = paste0("^part_[0-9]+_ts_", site, "_.*\\.rds$"), full.names = TRUE)
  keys <- sub("^part_[0-9]+_", "", basename(files))
  keys <- sub("\\.rds$", "", keys)
  split_files <- split(files, keys)
  for (key in names(split_files)) {
    message("Merging group: ", key)
    file_list <- split_files[[key]]
    data_list <- lapply(file_list, readRDS)
    merged <- do.call(rbind, data_list)
    out_file <- file.path(ts_dir, paste0(key, ".rds"))
    write_rds(merged, out_file, compress = "gz")
    file.remove(file_list)
  }
  
  # message("Cleaning time series...")
  # clean_planetscope_time_series_batch(dir = dir_data, v_site = site, num_cores = 12, calculate_evi = T) #

  # message("Processing time series...")
  # df_thres <- set_thresholds(thres_up = 0.5, thres_down = 0.5)
  # calculate_phenological_metrics_batch(dir = dir_data, v_year = 2017:2023, v_site = site, df_thres = df_thres, var_index = "evi",
  #                                      min_days = 1, check_seasonality = T, extend_to_previous_year = 275,
  #                                      extend_to_next_year = 90, num_cores = 30)
}

# > ?retrieve_planetscope_time_series_batch
# Warning messages:
#   1: In browseURL(paste0("http://127.0.0.1:", port, "/library/", pkgname,  :
#                            closing unused connection 16 (<-localhost:11679)
#                          2: In browseURL(paste0("http://127.0.0.1:", port, "/library/", pkgname,  :
#                                                   closing unused connection 15 (<-localhost:11679)

# 
# 
# clean_28117_allGroup <- clean_28117_allGroup |>
#   filter(year!=2016)
# saveRDS(clean_28117_allGroup,"C:/Users/feixi/Desktop/PS_NPN/CONUS/clean/clean_28117_allGroup.rds" )

