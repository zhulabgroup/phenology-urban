tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))

source(("~/phenology-urban/script/spatial_stats_model/prepare_model_input.R"))

phe_type <- "SOS"
taxon_col <- "canonicalname"
phe_metric <- "doy"

all_city_input <- lapply(cities, function(city) {
  message(paste0("Reading data in ", city))
  prepare_data(city, source = "Grid", taxa_col= taxon_col, phe_type = phe_type) |>
    mutate(city = city) |> 
    group_by(.data[[taxon_col]]) |>
    mutate(
      q1 = stats::quantile(.data[[phe_metric]], probs = 0.25, na.rm = TRUE),
      q3 = stats::quantile(.data[[phe_metric]], probs = 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      lower_bound = q1 - 1.5 * iqr,
      upper_bound = q3 + 1.5 * iqr
    ) |>
    filter(is.na(.data[[phe_metric]]) |
             dplyr::between(.data[[phe_metric]], lower_bound, upper_bound)) |>
    select(-q1, -q3, -iqr, -lower_bound, -upper_bound)
  })|>
  bind_rows()|>
  mutate(genus = str_extract(.data[[taxon_col]], "^[^ ]+"))


write_rds(all_city_input, "~/phenology-urban/data/proc/model_output/nooutlier_SOS_input.rds", compress = "gz")