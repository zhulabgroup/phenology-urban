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
  phe_types <- c("SOS", "EOS")
  if (model_type %in% c("inla", "inla_early_season", "inla_later_season", "inla_pre30", "inla_pre60")) {
    model_dir <- here("data/intermediate/spatial/", model_type, "spatial_spde/")
  } else if (model_type == "mcmc") {
    model_dir <- here("data/intermediate/spatial/mcmc/spatial_mat/")
  } else if (model_type == "freq") {
    model_dir <- here("data/intermediate/spatial/freq/")
  }  else if (model_type == "inla_gsl") {
    model_dir <- here("data/intermediate/spatial/inla_gsl/spatial_spde/")
    phe_types <- c("GSL_winter", "GSL_summer")
  }
  
  if (file.exists(save_summary_path)) {
    coef_city_species <- readr::read_rds(save_summary_path)
  } else {  
    coef_summary_all <- data.frame()
    for (city in cities) {
      r_output_dir <- paste0(model_dir, city, "/canonicalname/")
      for (phe_type in phe_types) {
        file_path <- file.path(r_output_dir, paste0(phe_type, "_summary.rds"))
        if (file.exists(file_path)) {
          coef_summary <- readr::read_rds(file_path) |>
            dplyr::mutate(phe_type = phe_type, city = city)
          coef_summary_all <- dplyr::bind_rows(coef_summary_all, coef_summary)
        }
      }
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
    write_rds(coef_city_species, save_summary_path)
  }
  return(coef_city_species)
}