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
library(rgbif)
library(bridgesampling)

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
# r_output_dir <- "~/phenology-urban/data/proc/urban/sp_model/fixed_preseason/random_sp_humid/"
# city <- "NY"
# phe_type <- "SOS"
# taxon_col <- "canonicalname"

#################### Read data ####################

# Given the temporal variation required, we will use the WU data
# Try the HUMID one

## prepare input data

metadata <- read_csv("~/lab-data/datasets/vegetation/PS/urban/metadata.csv") %>%
  filter(site == !!city)
file_list <- list.files(path = "~/lab-data/datasets/vegetation/PS/urban/doy", pattern = paste0("^doy_",city,"_.*\\.rds$"), full.names = TRUE)
data_list <- file_list %>%
  map(~ readRDS(.x))
all_doy <- bind_rows(data_list)
tree_location <- metadata %>%
  mutate(genus = str_extract(taxa, "^[^ ]+")) %>%
  left_join(all_doy, by = "id") %>%
  filter(!is.na(doy)) %>%
  mutate(geometry = paste0("c(", lon, ", ", lat, ")"))

points_sf <- st_as_sf(tree_location |> distinct(lon, lat), coords = c("lon", "lat"), crs = 4326, remove = F) |>
  st_transform(crs = 9311)
hex_size <- 5000
bbox <- st_bbox(points_sf)
hex_grid_geom <- st_make_grid(bbox, cellsize = hex_size, square = FALSE)
hex_grid <- st_sf(geometry = hex_grid_geom) |>
  mutate(hex_id = row_number()) |>
  mutate(centroid = st_centroid(geometry)) |>
  mutate(
    hex_lon = st_coordinates(centroid)[, 1],
    hex_lat = st_coordinates(centroid)[, 2]
  ) |>
  select(-centroid)
# Assign each point to a hexagon
tree_location_hex <- st_join(points_sf, hex_grid, join = st_within) %>%
  mutate(hex_id = as.factor(hex_id)) |>
  select(-hex_lon, -hex_lat) |>
  st_drop_geometry() |>
  right_join(tree_location, by = c("lon", "lat"))



# List of Humid variable names and their corresponding filenames
files <- c("temp_monthly_all_tree.rds",
           "srad_monthly_all_tree.rds")
base_path <- "~/phenology-urban/data/raw/"
processed_humid_list <- lapply(files, function(file_name) {
  # Read the file
  var_name <- sub("_monthly_all_tree.rds", "", file_name)
  data <- readRDS(paste0(base_path, city, "/HUMID/", file_name)) %>%
    select(id, value, date) %>%
    mutate(
      season = case_when(
        substr(date, 6, 7) %in% c("03", "04", "05") ~ "spring",
        substr(date, 6, 7) %in% c("06", "07", "08") ~ "summer",
        substr(date, 6, 7) %in% c("09", "10", "11") ~ "fall",
        substr(date, 6, 7) %in% c("12", "01", "02") ~ "winter"
      )) %>%
    mutate(year = as.integer(substr(date, 1, 4)),
           year = if_else(substr(date, 6, 7) == "12", year + 1, year)) %>%
    filter(year < 2019) %>% 
    group_by(id, year, season) %>%
    summarise(
      # Use sum for prcp, and mean for temp and srad
      Avg_value = mean(as.numeric(value), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Rename the Avg_value column to the variable name (e.g., tmax, tmin, prcp)
    rename(!!var_name := Avg_value)
  return(data)
})

files <- c("prcp_monthly_all_tree.rds")
processed_daymet_list <- lapply(files, function(file_name) {
  # Read the file
  var_name <- sub("_monthly_all_tree.rds", "", file_name)
  data <- readRDS(paste0(base_path, city, "/Daymet/monthly2016-2024/", file_name)) %>%
    select(id, value, date) %>%
    mutate(
      season = case_when(
        substr(date, 6, 7) %in% c("03", "04", "05") ~ "spring",
        substr(date, 6, 7) %in% c("06", "07", "08") ~ "summer",
        substr(date, 6, 7) %in% c("09", "10", "11") ~ "fall",
        substr(date, 6, 7) %in% c("12", "01", "02") ~ "winter"
      )) %>%
    mutate(year = as.integer(substr(date, 1, 4)),
           year = if_else(substr(date, 6, 7) == "12", year + 1, year)) %>%
    group_by(id, year, season) %>%
    summarise(
      # Use sum for prcp, and mean for temp and srad
      Avg_value = sum(as.numeric(value), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Rename the Avg_value column to the variable name (e.g., tmax, tmin, prcp)
    rename(!!var_name := Avg_value)
  
  return(data)
})

all_sites_climate <- processed_humid_list[[1]] |>
  left_join(processed_daymet_list[[1]], by = c("id", "season", "year")) 

srad_daily <- processed_humid_list[[2]] |>
  rename(Avg_srad = srad)

source("~/phenology-urban/script/spatial_stats_model/prepare_species_data.R")

if (taxon_col  == "canonicalname"){
  model_input_all <- prepare_species_data(phe_type, tree_location_hex, all_sites_climate, srad_daily, source = "Grid")
} else if (taxon_col == "genus") {
  model_input_all <- prepare_genus_data(phe_type, tree_location_hex, all_sites_climate, srad_daily, source = "Grid")
} else {
  stop("Invalid taxon_col. Must be 'canonicalname' or 'species'.")
}


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
    cat(sprintf("\n File already exists for taxon: %s, phe_type: %s, city: %s. Skipping model run.\n", 
                taxon, phe_type, city))
    next  # Skip the current iteration of the loop
  } # if the file exists, skip to the next iteration  
  
  cat(sprintf("\n Running model for taxon: %s, phe_type: %s, city: %s\n", 
              taxon, phe_type, city))
  model_input_sf <- model_input_all %>%
    filter(!!sym(taxon_col) == !!taxon) %>%
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
  ols_model <- lm(sp_formula, data = tavg_phe_df)
  # Define Nimble model code
  model_code <- nimbleCode({
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- inprod(X[i, 1:P], beta[]) + year_effect[year[i]] + sp_effect[hex_id[i]]
      y_pred[i] ~ dnorm(mu[i], tau)
    }
    for (q in 1:nyears) {
      z_year[q] ~ dnorm(0, 1)                     
      year_effect[q] <- z_year[q] * sigma_year   
    }
    sigma_year ~ dexp(1) 
    for (s in 1:nhex) {
      z_hex[s] ~ dnorm(0, 1)  
      sp_effect[s] <- z_hex[s] * sigma_hex   
    }
    sigma_hex ~ dexp(1) 
    beta[1] ~ dnorm(180, sd = 40)    
    beta[2] ~ dnorm(0, sd = 5) 
    beta[3] ~ dnorm(0, sd = 1) 
    beta[4] ~ dnorm(0, sd = 1) 
    tau ~ dgamma(0.001, 0.001)
  })
  
  # Data inputs
  
  N <- nrow(tavg_phe_df)
  nyears <- length(unique(tavg_phe_df$year))
  nhex <- length(unique(tavg_phe_df$hex_id))
  X <- model.matrix(sp_formula, data = tavg_phe_df)
  y <- tavg_phe_df[[all.vars(sp_formula)[1]]]
  year_index <- as.integer(factor(tavg_phe_df$year))  # Create fixed year indices
  hex_index <- as.integer(factor(tavg_phe_df$hex_id))  # Create fixed year indices
  data <- list(y = y, X = X)  # Remove 'year' from data
  inits <- list(
    beta = as.numeric(ols_model$coefficients),
    tau = 1, 
    sigma_year = 1,
    sigma_hex = 1,
    year_effect = rep(0, nyears),
    sp_effect = rep(0, nhex),
    mu = rep(180, N),
    y_pred = rep(180, N),  # Initialize y_pred
    z_year = rep(0, nyears),  # Initialize year effects
    z_hex = rep(0, nhex)  # Initialize spatial effects
  )
  
  ### Change based on semi-variogram fitting
  constants <- list(N = N, P = ncol(X), 
                    nyears = nyears, 
                    year = year_index,
                    nhex = nhex, 
                    hex_id = hex_index,
                    I = diag(1,N))
  monitors <- c("beta", "year_effect", "sp_effect")
  
  ni <- 300000   # nb iterations 
  nt <- 10       # thinning interval
  nb <- 250000    # nb iterations as burn-in 
  nc <- 4         # nb chains

  
  #################### Run/Save the model ####################
  
  my_model <- nimbleModel(code = model_code, data = data, inits = inits, constants = constants)
  compiled_model <- compileNimble(my_model)
  
  conf <- configureMCMC(my_model, enableWAIC = TRUE)
  conf$removeSampler(paste0("beta[", 1:4, "]"))
  conf$addSampler(target = paste0("beta[", 1:4, "]"), type = "RW_block")
  conf$removeSampler("tau")
  conf$addSampler(target = "tau", type = "slice")
  conf$addMonitors(monitors)
  # conf$addMonitors2("y_pred")  # Add y_pred to monitors
  
  mymcmc <- buildMCMC(conf)
  compiled_mcmc <- compileNimble(mymcmc, project = compiled_model)
  modelMCMCrun <- runMCMC(compiled_mcmc,
                          niter = ni, nburnin = nb, thin = nt, nchains = nc,
                          samplesAsCodaMCMC = T, summary = TRUE, WAIC = TRUE)
  
  beta2_mcmc <- lapply(modelMCMCrun$samples, function(chain) {
    as.mcmc(chain[, "beta[2]", drop = FALSE])
  }) |> as.mcmc.list()
  converge <- gelman.diag(beta2_mcmc, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
              multivariate=TRUE)
  
  # library(ggmcmc)
  # beta_label <- data.frame(Parameter = "beta[2]", Label = "AvgTemp")
  # param_data <- ggs(modelMCMCrun$samples, par_labels = beta_label) %>%
  #   filter(Parameter == "AvgTemp")
  # beta_mean <- mean(param_data$value)
  # param_data %>%
  #   ggs_density() + # traceplot
  #   theme_bw() +
  #   ggtitle(paste("taxon:", taxon)) +
  #   theme(plot.title = element_text(hjust = 0.5)) +
  #   annotate("text", x = beta_mean, y = 0, label = paste0("Mean: ", round(beta_mean, 2)),vjust = -1.5, color = "red", size = 4)
  # 
  
  
  
  log_posterior_nimble <- function(samples_row, data) {
    constants <- data$constants
    param_names <- names(samples_row)
    if (is.null(param_names)) stop("samples_row does not have named parameters")
    
    # unpack parameter values
    beta <- samples_row[grep("^beta\\[", param_names)]
    tau <- samples_row["tau"]
    sigma_year <- samples_row["sigma_year"]
    sigma_hex <- samples_row["sigma_hex"]
    
    z_year <- samples_row[grep("^z_year\\[", param_names)]
    z_hex  <- samples_row[grep("^z_hex\\[",  param_names)]
    
    # log-prior
    log_prior <- sum(
      dnorm(beta[1], mean = 180, sd = 40, log = TRUE),
      dnorm(beta[2], 0, 5, log = TRUE),
      dnorm(beta[3], 0, 1, log = TRUE),
      dnorm(beta[4], 0, 1, log = TRUE),
      dexp(sigma_year, rate = 1, log = TRUE),
      dexp(sigma_hex, rate = 1, log = TRUE),
      dgamma(tau, shape = 0.001, rate = 0.001, log = TRUE),
      sum(dnorm(z_year, 0, 1, log = TRUE)),
      sum(dnorm(z_hex, 0, 1, log = TRUE))
    )
    
    # log-likelihood
    X <- data$X
    y <- data$y
    year <- constants$year
    hex_id <- constants$hex_id
    N <- constants$N
    
    year_effect <- z_year * sigma_year
    sp_effect   <- z_hex  * sigma_hex
    
    mu <- as.vector(X %*% beta) + year_effect[year] + sp_effect[hex_id]
    log_lik <- sum(dnorm(y, mean = mu, sd = 1 / sqrt(tau), log = TRUE))
    
    return(log_lik + log_prior)
  }
  
  posterior_samples <- as.matrix(modelMCMCrun$samples)
  param_names <- colnames(posterior_samples)
  model_bridge <- bridge_sampler(samples = posterior_samples,
                              log_posterior = log_posterior_nimble,
                              data = list(y = y, X = X, constants = constants),
                              lb=setNames(rep(-Inf, length(param_names)), param_names),
                              ub=setNames(rep( Inf, length(param_names)), param_names),
                              param_types = rep("real", ncol(posterior_samples)),
                              repetitions = 10,
                              method = "normal",
                              cores = 1,
                              use_neff = T,
                              silent = FALSE,
                              verbose = TRUE,
                              varlist = ls(),
                              envir = environment())
  mean_logml <- mean(model_bridge$logml)
  
  
  
  modelMCMCrun <- list(
    samples = modelMCMCrun$samples,
    summary = modelMCMCrun$summary,
    WAIC = modelMCMCrun$WAIC,
    converge = converge$psrf[1],
    mean_logml = mean_logml
  )
  write_rds(modelMCMCrun, file = model_path, compress = "gz")
}


