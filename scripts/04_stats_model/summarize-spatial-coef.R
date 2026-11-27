### City and species level data
cities <- names(city_label)
inventory_centroids <- read_rds(here("data/input/tree_list/summary/inventory_centroids.rds")) |>
  mutate(city_label = unname(city_label[site]))
boundary_us <- st_read(here("data/input/boundary/US/us_boundarydissolve.shp"), quiet = T) |>
  st_transform(4326)

load(here("data/input/tree_list/summary/tree_species.RData"))
city_inven_backbone <- read_rds(here("data/input/tree_list/summary/city_inven_backbone.rds"))|>
  filter(
    kingdom == "Plantae", phylum == "Tracheophyta",
    rank %in% c("GENUS", "SPECIES"),
    status %in% c("ACCEPTED", "SYNONYM")
  )|>
  mutate(usageKey = if_else(status == "SYNONYM", acceptedUsageKey, usageKey))

city_feature <- read_rds(here("data/input/tree_list/summary/city_feature.rds")) |>
  select(-`ID.1`) |>
  left_join(read_rds(here("data/input/tree_list/summary/city_bioclimate.rds")) |>
              select(site, bio1, bio12), by = "site") |>
  rename(city = site) |>
  mutate(city = ifelse(city == "St_George", "StGeorge", city)) |>
  mutate(city = ifelse(city == "St.Louis", "StLouis", city))|>
  mutate(city_label = unname(city_label[city]))

### Coefficients results
summarize_spatial_coef <- function(model_type = "inla"){
  save_summary_path <- paste0(here("data/output/spatial_coef_summary_"), model_type, ".rds")
  phe_types <- c("SOS", "EOS")s
  if (model_type %in% c("inla", "inla_early_season", "inla_later_season", "inla_pre30", "inla_pre60")) {
    model_dir <- here("data/intermediate/spatial/", model_type, "spatial_spde/")
  } else if (model_type == "mcmc") {
    model_dir <- here("data/intermediate/spatial/mcmc/spatial_mat/")
  } else if (model_type == "freq") {
    model_dir <- here("data/intermediate/spatial/freq/")
  }  else if (model_type == "inla_gsl") {
    model_dir <- here("data/intermediate/spatial/inla_gsl/spatial_spde/")
    phe_types <- c("GSL_winter", "GSL_summer")
  } else {
    stop("Unknown model_type: ", model_type)
  }

  if (file.exists(save_summary_path)) {
    coef_city_species <- readr::read_rds(save_summary_path)
  } else {
    coef_summary_all <- data.frame()
    read_errors <- list()
    source(here("scripts/04_stats_model/fit_spatial/tidy_model_results.R"))
    for (city in cities) {
      r_output_dir <- paste0(model_dir, city, "/canonicalname/")
      for (phe_type in phe_types) {
        file_path <- file.path(r_output_dir, paste0(phe_type, "_summary.rds"))
        record_error <- function(path, error) {
          read_errors[[length(read_errors) + 1L]] <<- data.frame(
            city = city, phe_type = phe_type, file = path,
            error = conditionMessage(error), stringsAsFactors = FALSE
          )
        }
        coef_summary <- NULL
        if (file.exists(file_path)) {
          coef_summary <- tryCatch({
            cached <- readr::read_rds(file_path)
            if (!is.data.frame(cached) || !"Factor" %in% names(cached)) {
              stop("Invalid coefficient summary")
            }
            cached
          }, error = function(e) {
            record_error(file_path, e)
            NULL
          })
        }
        if (is.null(coef_summary)) {
          message(paste0("Processing ", city, " Phenophase: ", phe_type))
          coef_summary <- tryCatch({
            if (!dir.exists(r_output_dir)) stop("Model directory does not exist")
            if (model_type != "mcmc") {
              prepare_coef_summary(r_output_dir, phe_type, on_error = record_error)
            } else {
              prepare_coef_summary_mcmc(r_output_dir, phe_type, on_error = record_error)
            }
          }, error = function(e) {
            record_error(r_output_dir, e)
            NULL
          })
        }
        if (!is.null(coef_summary) && nrow(coef_summary) > 0) {
          coef_summary_all <- tryCatch({
            coef_summary <- dplyr::mutate(coef_summary, phe_type = phe_type, city = city)
            dplyr::bind_rows(coef_summary_all, coef_summary)
          }, error = function(e) {
            record_error(if (file.exists(file_path)) file_path else r_output_dir, e)
            coef_summary_all
          })
        }
      }
    }
    errors <- dplyr::bind_rows(read_errors)
    if (nrow(errors) > 0) {
      error_path <- tempfile(paste0("spatial_coef_summary_", model_type, "_errors_"),
                             tmpdir = dirname(save_summary_path), fileext = ".csv")
      readr::write_csv(errors, error_path)
      warning(nrow(errors), " coefficient read error(s); details saved to ", error_path)
    }
    if (!nrow(coef_summary_all)) {
      stop("No coefficient summaries could be loaded")
    }
    coef_summary_all <- coef_summary_all |> 
      mutate(phe_type = factor(phe_type, levels = phe_types),
             genus = sub(" .*", "", Factor),
             species = Factor) |>
      mutate(city_label = unname(city_label[city]))
    
    coef_city_species <- coef_summary_all |>
      mutate(Factor := gsub("_", " ", species)) |>
      left_join(city_inven_backbone |>
                  select(usageKey, verbatim_name), by = c("Factor" = "verbatim_name")) |>
      left_join(city_feature |> 
                  select(-city_label) |>
                  st_drop_geometry(), 
                by = "city")
    if (!nrow(errors)) write_rds(coef_city_species, save_summary_path)
    attr(coef_city_species, "read_errors") <- errors
  }
  return(coef_city_species)
}
