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
# r_output_dir <- "~/phenology-urban/data/proc/urban/sp_model/fixed_preseason/spbayes/"
# city <- "NY"
# phe_type <- "EOS"
# taxon_col <- "canonicalname"

#################### Read data ####################

# Given the temporal variation required, we will use the WU data

## prepare input data
points_in_buffer <- read_csv(paste0("~/phenology-urban/data/proc/urban/",city,"/tree_WU_500_buffer_PS.csv"), show_col_types = FALSE)

points_sf <- st_as_sf(points_in_buffer, coords = c("lon", "lat"), crs = 4326, remove = F) |>
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
points_in_buffer <- st_join(points_sf, hex_grid, join = st_within) %>%
  mutate(hex_id = as.factor(hex_id)) |>
  select(-hex_lon, -hex_lat) |>
  st_drop_geometry()


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

if (taxon_col  == "canonicalname"){
  source("script/spatial_stats_model/run_species_inla.R")
  model_input_all <- prepare_species_data(phe_type, points_in_buffer, all_sites_temp, srad_daily, source = "WU")
} else if (taxon_col == "genus") {
  source("script/spatial_stats_model/run_genus_inla_fixpre.R")
  model_input_all <- prepare_genus_data(phe_type, points_in_buffer, all_sites_temp, srad_daily, source = "WU")
} else {
  stop("Invalid taxon_col. Must be 'canonicalname' or 'species'.")
}


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
  tavg_phe <- model_input_sf
  temp_col <- names(tavg_phe)[grepl("^temp", names(tavg_phe))]
  prcp_col <- names(tavg_phe)[grepl("^prcp", names(tavg_phe))]
  srad_col <- names(tavg_phe)[grepl("^srad", names(tavg_phe))]
  
  #################### Run the spatial model ####################
  coords <- st_coordinates(tavg_phe)
  dist_matrix <- as.matrix(dist(coords))
  N_knots = 50
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
  ols_model <- lm(sp_formula, data = tavg_phe_df) 
  # Define Nimble model code
  model_code <- nimbleCode({
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- inprod(X[i, 1:P], beta[]) + spatial_effect[i] + year_effect[year[i]]
    }
    for (q in 1:nyears) {
      z_year[q] ~ dnorm(0, 1)                     
      year_effect[q] <- z_year[q] * sigma_year   
    }
    sigma_year ~ dexp(1) 
    beta[1] ~ dnorm(180, sd = 40)    
    beta[2] ~ dnorm(0, sd = 5) 
    beta[3] ~ dnorm(0, sd = 1) 
    beta[4] ~ dnorm(0, sd = 1) 
    spatial_effect[1:N] ~ dmnorm(Mu[1:N], cov = cov_matrix[1:N, 1:N])
    # spatial_effect[1:N] ~ dmnorm_cholesky(Mu[1:N], chol = chol_cov_matrix[1:N, 1:N])
    
    # Choose the covariance function based on variogram fitting
    if (N_knots <= N_unique){
      knots_cov_matrix[1:N_knots, 1:N_knots] <- sigma.sq * exp(-phi * knots_dist_matrix[1:N_knots, 1:N_knots]^2) + I[1:N_knots, 1:N_knots]
      data_knots_cov_matrix[1:N, 1:N_knots] <- sigma.sq * exp(-phi * data_knots_dist_matrix[1:N, 1:N_knots]^2)
      knots_cov_chol_inv[1:N_knots, 1:N_knots] <- inverse(chol(knots_cov_matrix[1:N_knots, 1:N_knots]))
      cov_matrix[1:N, 1:N] <- data_knots_cov_matrix[1:N, 1:N_knots] %*% 
        knots_cov_chol_inv[1:N_knots, 1:N_knots] %*% t(knots_cov_chol_inv[1:N_knots, 1:N_knots]) %*% 
        t(data_knots_cov_matrix[1:N, 1:N_knots]) + I[1:N, 1:N]
    } else {
      cov_matrix[1:N, 1:N] <- sigma.sq * exp(-phi * dist_matrix[1:N, 1:N]^2) + I[1:N, 1:N]
    }
      
    tau ~ dexp(1)
    sigma.sq ~ dexp(1)
    phi ~ dunif(5 / 1, 5 / 0.1)
    range ~ dunif(0.01, 1)
  })
  
  # Data inputs
  
  nyears <- length(unique(tavg_phe_df$year))
  X <- model.matrix(sp_formula, data = tavg_phe_df)
  y <- tavg_phe_df[[all.vars(sp_formula)[1]]]
  year_index <- as.integer(factor(tavg_phe_df$year))  # Create fixed year indices
  data <- list(y = y, X = X)  # Remove 'year' from data
  inits <- list(
    beta = rep(0, ncol(X)),
    tau = 1, 
    spatial_effect = rep(0, N), 
    sigma_year = 1,
    year_effect = rep(0, nyears),
    
    sigma.sq = 10, 
    range = 0.2,
    phi = 10
  )
  
  ### Change based on semi-variogram fitting
  constants <- list(N = N, N_knots = N_knots, N_unique = N_unique, P = ncol(X), Mu = rep(0, N), nyears = nyears, 
                    dist_matrix = dist_matrix, knots_dist_matrix = knots_dist_matrix, 
                    data_knots_dist_matrix = data_knots_dist_matrix,
                    year = year_index, I = diag(1e-6,N))  # Add 'year' to constants
  monitors <- c("beta", "tau", "sigma.sq", "phi", "spatial_effect", "year_effect")
  
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
  for (p in c("sigma.sq", "range", "tau", "phi")) {
    conf$removeSampler(p)
    conf$addSampler(target = p, type = "slice")
  }
  conf$addMonitors(monitors)
  
  mcmc <- buildMCMC(conf)
  compiled_mcmc <- compileNimble(mcmc, project = compiled_model)
  modelMCMCrun <- runMCMC(compiled_mcmc,
                          niter = ni, nburnin = nb, thin = nt, nchains = nc,
                          samplesAsCodaMCMC = TRUE, summary = TRUE, WAIC = TRUE)
  
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
    
    
    beta <- samples_row[grep("^beta\\[", param_names)]
    tau <- samples_row["tau"]
    sigma_year <- samples_row["sigma_year"]
    sigma_sq <- samples_row["sigma.sq"]
    phi <- samples_row["phi"]
    spatial_effect <- samples_row[grep("^spatial_effect\\[", param_names)]
    z_year <- samples_row[grep("^z_year\\[", param_names)]
    
    
    log_prior <- sum(
      dnorm(beta[1], mean = 180, sd = 40, log = TRUE),
      dnorm(beta[2], 0, 5, log = TRUE),
      dnorm(beta[3], 0, 1, log = TRUE),
      dnorm(beta[4], 0, 1, log = TRUE),
      if (length(beta) > 4) sum(dnorm(beta[5:length(beta)], 0, 10, log = TRUE)) else 0,
      dexp(sigma_year, rate = 1, log = TRUE),
      dexp(sigma_sq, rate = 1, log = TRUE),
      dunif(phi, 5/1, 5/0.1, log = TRUE),
      dgamma(tau, shape = 0.001, rate = 0.001, log = TRUE),
      sum(dnorm(z_year, 0, 1, log = TRUE))
    )
    
    
    coords <- data$coords
    N <- constants$N
    N_knots <- data$N_knots
    N_unique <- nrow(unique(coords))
    
    if (N_knots <= N_unique) {
      # --- Knots-based approximation ---
      knots <- data$knots
      knots_dist_matrix <- data$knots_dist_matrix
      data_knots_dist_matrix <- data$data_knots_dist_matrix
      
      
      K_kk <- sigma_sq * exp(-phi * knots_dist_matrix^2) + diag(1, N_knots)
      K_nk <- sigma_sq * exp(-phi * data_knots_dist_matrix^2)
      
      
      K_kk_chol_inv <- chol(K_kk)
      K_kk_inv <- chol2inv(K_kk_chol_inv)
      cov_approx <- K_nk %*% K_kk_inv %*% t(K_nk) + diag(1, N)
      
      log_prior <- log_prior + mvtnorm::dmvnorm(
        spatial_effect, mean = rep(0, N),
        sigma = cov_approx, log = TRUE
      )
    } else {
      # --- Full GP covariance ---
      dist_matrix <- data$dist_matrix
      cov_matrix <- sigma_sq * exp(-phi * dist_matrix^2) + diag(1, N)
      log_prior <- log_prior + mvtnorm::dmvnorm(
        spatial_effect, mean = rep(0, N),
        sigma = cov_matrix, log = TRUE
      )
    }
    
    # 4. log-likelihood
    X <- data$X
    y <- data$y
    year <- constants$year
    year_effect <- z_year * sigma_year
    
    mu <- as.vector(X %*% beta) + year_effect[year] + spatial_effect
    log_lik <- sum(dnorm(y, mean = mu, sd = 1 / sqrt(tau), log = TRUE))
    
    return(log_lik + log_prior)
  }
  
  if (N_knots <= N_unique) {
    knots <- kmeans(coords, centers = N_knots)$centers
    knots_dist_matrix <- as.matrix(dist(knots))
    data_knots_dist_matrix <- as.matrix(dist(rbind(coords, knots)))[1:N, (N + 1):(N + N_knots)]
  } else {
    knots <- coords
    knots_dist_matrix <- dist_matrix
    data_knots_dist_matrix <- dist_matrix
  }
  
  
  data_bridge <- list(
    y = y,
    X = X,
    coords = coords,
    dist_matrix = dist_matrix,
    knots = knots,
    knots_dist_matrix = knots_dist_matrix,
    data_knots_dist_matrix = data_knots_dist_matrix,
    N_knots = N_knots,
    constants = constants
  )
  
  posterior_samples <- as.matrix(modelMCMCrun$samples)
  param_names <- colnames(posterior_samples)
  
  model_bridge <- bridge_sampler(
    samples = posterior_samples,
    log_posterior = log_posterior_nimble,
    data = data_bridge,
    lb = setNames(rep(-Inf, length(param_names)), param_names),
    ub = setNames(rep( Inf, length(param_names)), param_names),
    param_types = rep("real", ncol(posterior_samples)),
    repetitions = 10,
    method = "normal",
    cores = 1,
    use_neff = TRUE,
    silent = FALSE,
    verbose = TRUE,
    varlist = ls(),
    envir = environment()
  )
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


