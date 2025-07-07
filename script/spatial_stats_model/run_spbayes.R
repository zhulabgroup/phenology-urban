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
# r_output_dir <- "~/phenology-urban/data/proc/urban/sp_model/fixed_preseason/spbayes_species_humid_all/"
# city <- "NY"
# phe_type <- "SOS"
# taxon_col <- "canonicalname"

#################### Read data ####################

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

source("script/spatial_stats_model/run_species_inla.R")
model_input_all <- prepare_species_data(phe_type, tree_location, all_sites_climate, srad_daily)

taxon_list <- unique(model_input_all[[taxon_col]])
#################### Find the optimal cov model ####################
for (taxon in taxon_list) {
  taxon_clean <- gsub(" ", "_", taxon)
  out_dir <- paste0(r_output_dir, city, "/")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  model_path <- paste0(out_dir, taxon_clean, "_", phe_type, ".rds")
  
  if (file.exists(model_path)) {
    next  # Skip the current iteration of the loop
  }
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
  tavg_phe <- model_input_sf |>
    st_drop_geometry()
  temp_col <- names(tavg_phe)[grepl("^temp", names(tavg_phe))]
  prcp_col <- names(tavg_phe)[grepl("^prcp", names(tavg_phe))]
  srad_col <- names(tavg_phe)[grepl("^srad", names(tavg_phe))]
  
  # OLS_formula <- as.formula(paste0(phe_metric, " ~ ", temp_col, "+", prcp_col, "+", srad_col , "+ (1|year)"))
  # random_OLS <- lmer(OLS_formula, data = tavg_phe)
  # tavg_phe$residual_random <- tavg_phe[[phe_metric]] - fitted(random_OLS)
  # coordinates(tavg_phe) <- ~ lon + lat
  # vgm.emt <- variogram(residual_random ~ 1, tavg_phe, cloud = FALSE)
  # vgm.las.init <- vgm(model = c("Exp", "Sph", "Gau", "Mat"))
  # vgm_exn_fit <- fit.variogram(vgm.emt, vgm.las.init)
  # cov_type <- as.character(vgm_exn_fit$model)
  cov_type = "Mat"
  
  #################### Run the spatial model ####################
  coords <- coordinates(tavg_phe)
  dist_matrix <- as.matrix(dist(coords))
  N_knots = 100
  N <- nrow(tavg_phe)
  N_unique = nrow(unique(coords))
  if (N_knots <= N_unique){
    knots <- kmeans(coords, centers = N_knots)$centers
    knots_dist_matrix <- as.matrix(dist(knots))
    data_knots_dist_matrix <- as.matrix(dist(rbind(coords, knots)))[1:N, (N + 1):(N + N_knots)]
  } else {
    knots <- coords
    knots_dist_matrix <- dist_matrix
    data_knots_dist_matrix <- dist_matrix
  }
    
  
  tavg_phe_df <- as.data.frame(tavg_phe)
  sp_formula <- as.formula(paste0(phe_metric, " ~ ", temp_col, "+", prcp_col, "+", srad_col))
  
  # Define covariance functions
  
  spherical_covariance <- nimbleFunction(
    run = function(d = double(2), range = double(0), sigma.sq = double(0)) {
      returnType(double(2))  # Return type is a matrix
      nrow_d <- dim(d)[1]
      ncol_d <- dim(d)[2]
      result <- matrix(0, nrow_d, ncol_d)  # Initialize result matrix
      for (i in 1:nrow_d) {
        for (j in 1:ncol_d) {
          if (d[i, j] <= range) {
            result[i, j] = sigma.sq * (1 - (3 * d[i, j] / (2 * range)) + (d[i, j]^3 / (2 * range^3)))
          } else {
            result[i, j] = 0.0
          }
        }
      }
      return(result)
    }
  )
  
  matern_covariance <- nimbleFunction(
    run = function(d = double(2), range = double(0), sigma.sq = double(0), nu = double(0)) {
      returnType(double(2))
      nrow_d <- dim(d)[1]
      ncol_d <- dim(d)[2]
      result <- matrix(0, nrow_d, ncol_d)  # Initialize result matrix
      for (i in 1:nrow_d) {
        for (j in 1:ncol_d) {
            if (d[i, j] == 0) {
              result[i, j] = 0.0
          } else {
              part <- (2^(1 - nu)) / gamma(nu)  * (sqrt(2 * nu) * d[i, j] / range)^nu
              result[i,j] = sigma.sq * part * besselK(d[i, j] / range, nu)
          }
        }
      }  
      return(result)
    }
  )
  
  # Define Nimble model code
  model_code <- nimbleCode({
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- inprod(X[i, 1:P], beta[]) + spatial_effect[i] + year_effect[year[i]]
    }
    for (q in 1:nyears) {
      year_effect[q] ~ dnorm(0, sd = sqrt(tau.year))
    }
    for (j in 1:P) {
      beta[j] ~ dnorm(0, sd = 1000)
    }
    spatial_effect[1:N] ~ dmnorm(Mu[1:N], cov = cov_matrix[1:N, 1:N])
    
    
    # Choose the covariance function based on variogram fitting
    if (N_knots <= N_unique){
      if (cov_type == "Mat") {
        knots_cov_matrix[1:N_knots, 1:N_knots] <- matern_covariance(knots_dist_matrix[1:N_knots, 1:N_knots], range, sigma.sq, nu) + I[1:N_knots, 1:N_knots]
        data_knots_cov_matrix[1:N, 1:N_knots] <- matern_covariance(data_knots_dist_matrix[1:N, 1:N_knots], range, sigma.sq, nu)
      } else if (cov_type == "Gau") {
        knots_cov_matrix[1:N_knots, 1:N_knots] <- sigma.sq * exp(-phi * knots_dist_matrix[1:N_knots, 1:N_knots]^2) + I[1:N_knots, 1:N_knots]
        data_knots_cov_matrix[1:N, 1:N_knots] <- sigma.sq * exp(-phi * data_knots_dist_matrix[1:N, 1:N_knots]^2)
      } else if (cov_type == "Exp") {
        knots_cov_matrix[1:N_knots, 1:N_knots] <- sigma.sq * exp(-phi * knots_dist_matrix[1:N_knots, 1:N_knots]) + I[1:N_knots, 1:N_knots]
        data_knots_cov_matrix[1:N, 1:N_knots] <- sigma.sq * exp(-phi * data_knots_dist_matrix[1:N, 1:N_knots])
      # } else if (cov_type == "Sph") {
      #   for (m in 1:N_knots) {
      #     for (n in 1:N_knots) {
      #       knots_cov_matrix[m, n] <- spherical_covariance(knots_dist_matrix[m, n], range, sigma.sq)
      #     }
      #   }
      #   for (i in 1:N) {
      #     for (j in 1:N_knots) {
      #       data_knots_cov_matrix[i, j] <- spherical_covariance(data_knots_dist_matrix[i, j], range, sigma.sq)
      #     }
      #   }
      #   cov_matrix_ori[1:N, 1:N] <- data_knots_cov_matrix[1:N, 1:N_knots] %*% 
      #     inverse(knots_cov_matrix[1:N_knots, 1:N_knots]) %*% 
      #     t(data_knots_cov_matrix[1:N, 1:N_knots])
      #   cov_matrix[1:N, 1:N] <- cov_matrix_ori[1:N, 1:N] + I[1:N, 1:N]
      # }
      } else if (cov_type == "Sph") {
        knots_cov_matrix[1:N_knots, 1:N_knots] <- spherical_covariance(knots_dist_matrix[1:N_knots, 1:N_knots], range, sigma.sq) + I[1:N_knots, 1:N_knots]
        data_knots_cov_matrix[1:N, 1:N_knots] <- spherical_covariance(data_knots_dist_matrix[1:N, 1:N_knots], range, sigma.sq)
      }
      knots_cov_chol_inv[1:N_knots, 1:N_knots] <- inverse(chol(knots_cov_matrix[1:N_knots, 1:N_knots]))
      cov_matrix[1:N, 1:N] <- data_knots_cov_matrix[1:N, 1:N_knots] %*% 
        knots_cov_chol_inv[1:N_knots, 1:N_knots] %*% t(knots_cov_chol_inv[1:N_knots, 1:N_knots]) %*% 
        t(data_knots_cov_matrix[1:N, 1:N_knots]) + I[1:N, 1:N]
    } else {
      if (cov_type == "Mat") {
        cov_matrix[1:N, 1:N] <-  matern_covariance(dist_matrix[1:N, 1:N], range, sigma.sq, nu) + I[1:N, 1:N]
      } else if (cov_type == "Gau") {
        cov_matrix[1:N, 1:N] <- sigma.sq * exp(-phi * dist_matrix[1:N, 1:N]^2) + I[1:N, 1:N]
      } else if (cov_type == "Exp") {
        cov_matrix[1:N, 1:N] <- sigma.sq * exp(-phi * dist_matrix[1:N, 1:N]) + I[1:N, 1:N]
      } else if (cov_type == "Sph") {
        cov_matrix[1:N, 1:N] <- spherical_covariance(dist_matrix[1:N, 1:N], range, sigma.sq) + I[1:N, 1:N]
      }
    }
      
    tau ~ dgamma(0.001, 0.001)
    sigma.sq ~ dinvgamma(2.5, 120)
    phi ~ dunif(5 / 1, 5 / 0.1)
    tau.year ~ dinvgamma(2.001, 0.667)
    range ~ dunif(0, 1)
    nu ~ dunif(0.5, 5)
  })
  
  # Data inputs
  
  nyears <- length(unique(tavg_phe_df$year))
  X <- model.matrix(sp_formula, data = tavg_phe_df)
  y <- tavg_phe_df[[all.vars(sp_formula)[1]]]
  year_index <- as.integer(factor(tavg_phe_df$year))  # Create fixed year indices
  data <- list(y = y, X = X)  # Remove 'year' from data
  inits <- list(
    beta = rep(0, ncol(X)),
    tau = 1, sigma.sq = 80, phi = 25, tau.year = 1,
    range = 0.5,
    spatial_effect = rep(0, N), 
    nu = 1,
    year_effect = rep(0, nyears)
  )
  
  ### Change based on semi-variogram fitting
  constants <- list(N = N, N_knots = N_knots, N_unique = N_unique, P = ncol(X), Mu = rep(0, N), nyears = nyears, 
                    dist_matrix = dist_matrix, knots_dist_matrix = knots_dist_matrix, 
                    data_knots_dist_matrix = data_knots_dist_matrix,
                    year = year_index, I = diag(1e-6,N), cov_type = cov_type)  # Add 'year' to constants
  monitors <- switch(
    cov_type,
    "Mat" = c("beta", "tau.year", "range", "sigma.sq", "nu", "spatial_effect", "year_effect"),
    "Gau" = c("beta", "tau", "sigma.sq", "phi", "tau.year", "spatial_effect", "year_effect"),
    "Exp" = c("beta", "tau", "sigma.sq", "phi", "tau.year", "spatial_effect", "year_effect"),
    "Sph" = c("beta", "tau.year", "range", "sigma.sq", "spatial_effect", "year_effect"),
    stop("Unknown cov_type!")
  )
  
  ni <- 10000   # nb iterations 
  nt <- 50       # thinning interval
  nb <- ni / 2    # nb iterations as burn-in 
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
  # beta_median <- median(param_data$value)
  # param_data %>%
  #   ggs_density() +
  #   theme_bw() +
  #   ggtitle(paste("Genus:", genus)) +
  #   theme(plot.title = element_text(hjust = 0.5)) +
  #   annotate("text", x = beta_median, y = 0, label = paste0("Median: ", round(beta_median, 2)),vjust = -1.5, color = "red", size = 4)
  # 
  
  saveRDS(modelMCMCrun, file = model_path)
}


