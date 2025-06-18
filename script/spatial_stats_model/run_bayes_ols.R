library(nimble)
library(sp)
library(sf)
library(dplyr)
library(stringr)
library(readr)
library(coda)
library(lubridate)
library(tidyr)
library(taxize)
library(purrr)
library(terra)

# Read file paths from environment variables
r_input_dir <- Sys.getenv("R_INPUT_DIR")
r_output_dir <- Sys.getenv("R_OUTPUT_DIR")
city <- Sys.getenv("CITY")
phe_type <- Sys.getenv("PHE_TYPE")
taxon_col <- Sys.getenv("TAXON_COL")  # column name for species

# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) < 3) {
#   stop("5 arguments needed: r_output_file, city, phe_type")
# }
# 
# r_input_dir <- args[1]
# r_output_dir <- args[2]
# city <- args[3]
# phe_type <- args[4]

# r_input_dir <- "~/phenology-urban/data/"
# r_output_dir <- "~/phenology-urban/data/proc/urban/ols_model/bayes_species/"
# city <- "NY"
# phe_type <- "SOS"
# taxon_col <- "canonicalname"

#################### Read data ####################

points_in_buffer <- read_csv(paste0(r_input_dir, "proc/urban/",city,"/tree_WU_500_buffer_PS.csv"), show_col_types = FALSE)

all_sites_temp <- readRDS(paste0(r_input_dir, "proc/urban/",city,"/",city,"_wu.rds")) %>%
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

srad_daily <- readRDS(paste0(r_input_dir, "raw/",city,"/Daymet/daily2016-2024/srad_daily_2016-2023.rds")) %>%
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
source("script/spatial_stats_model/run_species_inla.R")
model_input_all <- prepare_species_data(phe_type, points_in_buffer, all_sites_temp, srad_daily)

#################### Find the optimal cov model ####################

taxon_list <- unique(model_input_all[[taxon_col]])
# boundary <- switch(
#   city,
#   "NY" = st_read("~/urban-cooling/data/raw/WU/NY/Untitled/nybb_dissolved.shp"),
#   "DT" = st_read("~/urban-cooling/data/raw/WU/DT/City_of_Detroit_Boundary/City_of_Detroit_Boundary.shp"),
#   "ST" = st_read("~/urban-cooling/data/raw/WU/ST/seattle_boundary/CityUGA.shp") %>%
#     filter(CITY_NM == "Seattle"),
#   "HT" = st_read("~/urban-cooling/data/raw/WU/HT/HGAC_City_Boundaries/HGAC_City_Boundaries.shp") %>%
#     filter(NAME == "Houston"),  # change to some distance
#   "TP" = st_read("~/urban-cooling/data/raw/WU/TP/Municipal_Boundary/Municipal_Boundary.shp") %>%
#     filter(NAME == "Tampa"),
#   "DV" = st_read("~/urban-cooling/data/raw/WU/DV/DRCOG_Municipalities/DRCOG_Municipalities.shp") %>%
#     filter(city == "Denver"),
#   "AT" = st_read("~/urban-cooling/data/raw/WU/AT/TxDOT_City_Boundaries/Cities.shp") %>%
#     filter(CITY_NM == "Austin"),
#   stop("Invalid city abbreviation.")
# )

for (taxon in taxon_list) {
  taxon_clean <- gsub(" ", "_", taxon)
  out_dir <- paste0(r_output_dir, city, "/")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  model_path <- paste0(out_dir, taxon_clean, "_", phe_type, ".rds")
  
  if (file.exists(model_path)) {
    cat(sprintf("\n File already exists for taxon: %s, phe_type: %s, city: %s. Skipping INLA run.\n", 
                taxon, phe_type, city))
    next  # Skip the current iteration of the loop
  } # if the file exists, skip to the next iteration  
  
  cat(sprintf("\n Running INLA for taxon: %s, phe_type: %s, city: %s\n", 
              taxon, phe_type, city))
  model_input_sf <- model_input_all %>%
    filter(taxon == !!taxon) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    # st_transform(crs = st_crs(boundary)) %>%
    mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2],
           year = as.factor(year)) 
  if (phe_type %in% c("SOS","EOS")) {
    phe_metric <- "doy"
  } else {
    phe_metric <- "doy_diff"
  }
  
  tavg_phe_df <- model_input_sf %>% 
    st_drop_geometry()
  
  #################### Run the spatial model ####################
  temp_col <- names(tavg_phe_df)[grepl("^temp", names(tavg_phe_df))]
  prcp_col <- names(tavg_phe_df)[grepl("^prcp", names(tavg_phe_df))]
  srad_col <- names(tavg_phe_df)[grepl("^srad", names(tavg_phe_df))]
  
  sp_formula <- as.formula(paste0(phe_metric, " ~ ", temp_col, "+", prcp_col, "+", srad_col))
  
  # Define Nimble model code
  model_code <- nimbleCode({
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- inprod(X[i, 1:P], beta[]) + year_effect[year[i]]
    }
    for (q in 1:nyears) {
      year_effect[q] ~ dnorm(0, sd = sqrt(tau.year))
    }
    for (j in 1:P) {
      beta[j] ~ dnorm(0, sd = 1000)
    }
    
    tau ~ dgamma(0.001, 0.001)
    tau.year ~ dinvgamma(2.001, 0.667)
  })
  
  # Data inputs
  
  N <- nrow(tavg_phe_df)
  nyears <- length(unique(tavg_phe_df$year))
  X <- model.matrix(sp_formula, data = tavg_phe_df)
  y <- tavg_phe_df[[all.vars(sp_formula)[1]]]
  year_index <- as.integer(factor(tavg_phe_df$year))  # Create fixed year indices
  data <- list(y = y, X = X)  # Remove 'year' from data
  inits <- list(
    beta = rep(0, ncol(X)),
    tau = 1, 
    tau.year = 1,
    year_effect = rep(0, nyears)
  )
  
  ### Change based on semi-variogram fitting
  constants <- list(N = N, P = ncol(X), 
                    nyears = nyears, 
                    year = year_index, 
                    I = diag(1,N))
  monitors <- c("beta", "tau.year", "year_effect")
  
  ni <- 20000   # nb iterations 
  nt <- 1       # thinning interval
  nb <- 5000    # nb iterations as burn-in 
  nc <- 1         # nb chains
  
  
  #################### Run/Save the model ####################
  
  modelMCMCrun <- nimbleMCMC(code = model_code, data = data, constants = constants, inits = inits, 
                             monitors = monitors,
                             niter = ni, nburnin = nb, thin = nt, nchains = nc, setSeed = 2025,
                             progressBar = TRUE, 
                             samplesAsCodaMCMC = TRUE, 
                             summary = TRUE, 
                             WAIC = TRUE)
  
  # beta_label <- data.frame(Parameter = "beta[2]", Label = "AvgTemp")
  # param_data <- ggs(modelMCMCrun$samples, par_labels = beta_label) %>%
  #   filter(Parameter == "AvgTemp")
  # beta_mean <- mean(param_data$value)
  # param_data %>%
  #   ggs_density() +
  #   theme_bw() +
  #   ggtitle(paste("taxon:", taxon)) +
  #   theme(plot.title = element_text(hjust = 0.5)) +
  #   annotate("text", x = beta_mean, y = 0, label = paste0("Mean: ", round(beta_mean, 2)),vjust = -1.5, color = "red", size = 4)
  # 
  saveRDS(modelMCMCrun, file = model_path)
}


