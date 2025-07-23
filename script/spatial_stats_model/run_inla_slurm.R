################################## package loading ##################################

library(zhulabtools)

load_packages(c("rgbif", "dplyr","lubridate", "ggplot2", "sf","sp", "readr", "purrr", "stringr","glue", "tidyr", "broom","lme4", "broom.mixed", "coda", "patchwork", "terra", "stats", "viridis"))
library(INLA)
################################## parameter loading ##################################

# Read file paths from environment variables
r_input_dir <- Sys.getenv("R_INPUT_DIR")
r_output_dir <- Sys.getenv("R_OUTPUT_DIR")
city <- Sys.getenv("CITY")
phe_type <- Sys.getenv("PHE_TYPE")
taxon_col <- Sys.getenv("TAXON_COL")  # column name for species

# r_input_dir <- "~/phenology-urban/data/"
# r_output_dir <- "~/phenology-urban/"
# city <- "NY"
# phe_type <- "EOS"
# taxon_col <- "canonicalname"  # column name for species

################################## data input ##################################

points_in_buffer <- read_csv(paste0("~/phenology-urban/data/proc/urban/",city,"/tree_WU_500_buffer_PS.csv"), show_col_types = FALSE)

all_sites_temp <- readRDS(paste0("~/urban-cooling/data_gl/raw/WU/",city,"/",city,"_wu.rds")) %>%
  dplyr::select(c(Date, HighTemp, AvgTemp, LowTemp, Sum_mm, name)) %>%
  mutate(season = case_when(
    month(Date) %in% 3:5 ~ "spring",
    month(Date) %in% 6:8 ~ "summer",
    month(Date) %in% 9:11 ~ "fall",
    month(Date) %in% c(12, 1, 2) ~ "winter"
  )) %>%
  mutate(year = year(Date),
         year = if_else(month(Date) == 12, year + 1, year)) %>%
  group_by(name, year, season) %>%
  summarise(
    Avg_AvgTemp = mean(AvgTemp, na.rm = TRUE),
    Sum_Sum_mm = sum(Sum_mm, na.rm = TRUE),
    AvgTemp_count = sum(!is.na(AvgTemp)),
    Sum_mm_count = sum(!is.na(Sum_mm)),
    .groups = "drop"
  )

srad_daily <- readRDS(paste0("~/phenology-urban/data/raw/",city,"/Daymet/daily2016-2024/srad_daily_2016-2023.rds")) %>%
  select(id, value, date)  %>%
  mutate(season = case_when(
    month(date) %in% 3:5 ~ "spring",
    month(date) %in% 6:8 ~ "summer",
    month(date) %in% 9:11 ~ "fall",
    month(date) %in% c(12, 1, 2) ~ "winter"
  )) %>%
  mutate(year = year(date),
         year = if_else(month(date) == 12, year + 1, year)) %>%
  group_by(id, year, season) %>%
  summarise(
    Avg_srad = mean(value, na.rm = TRUE),
    srad_count = sum(!is.na(value)),
    .groups = "drop"
  )

r_output_dir <- paste0(r_output_dir, city, "/")

if (taxon_col  == "canonicalname"){
  source("script/spatial_stats_model/run_species_inla.R")
  model_input_all <- prepare_species_data(phe_type, points_in_buffer, all_sites_temp, srad_daily, source = "WU")
  run_inla_species(model_input_all, r_output_dir, phe_type, city)
} else if (taxon_col == "genus") {
  source("script/spatial_stats_model/run_genus_inla_fixpre.R")
  model_input_all <- prepare_genus_data(phe_type, points_in_buffer, all_sites_temp, srad_daily, source = "WU")
  run_inla_genus(model_input_all, r_output_dir, phe_type, city)
} else {
  stop("Invalid taxon_col. Must be 'canonicalname' or 'species'.")
}

