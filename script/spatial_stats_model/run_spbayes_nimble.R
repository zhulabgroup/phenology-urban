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
  
  # OLS_formula <- as.formula(paste0(phe_metric, " ~ ", temp_col, "+", prcp_col, "+", srad_col , "+ (1|year)"))
  # random_OLS <- lmer(OLS_formula, data = tavg_phe)
  # tavg_phe$residual_random <- tavg_phe[[phe_metric]] - fitted(random_OLS)
  # coordinates(tavg_phe) <- ~ lon + lat
  # vgm.emt <- variogram(residual_random ~ 1, tavg_phe, cloud = FALSE)
  # vgm.las.init <- vgm(model = c("Exp", "Sph", "Gau", "Mat"))
  # vgm_exn_fit <- fit.variogram(vgm.emt, vgm.las.init)
  # cov_type <- as.character(vgm_exn_fit$model)
  cov_type = "Gau"
  
  #################### Run the spatial model ####################
  coords <- st_coordinates(tavg_phe)
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
  ols_model <- lm(sp_formula, data = tavg_phe_df) 
  # Define covariance functions
  # 
  # spherical_covariance <- nimbleFunction(
  #   run = function(d = double(2), range = double(0), sigma.sq = double(0)) {
  #     returnType(double(2))  # Return type is a matrix
  #     nrow_d <- dim(d)[1]
  #     ncol_d <- dim(d)[2]
  #     result <- matrix(0, nrow_d, ncol_d)  # Initialize result matrix
  #     for (i in 1:nrow_d) {
  #       for (j in 1:ncol_d) {
  #         if (d[i, j] <= range) {
  #           result[i, j] = sigma.sq * (1 - (3 * d[i, j] / (2 * range)) + (d[i, j]^3 / (2 * range^3)))
  #         } else {
  #           result[i, j] = 0.0
  #         }
  #       }
  #     }
  #     return(result)
  #   }
  # )
  # 
  # matern_covariance <- nimbleFunction(
  #   run = function(d = double(2), range = double(0), sigma.sq = double(0), nu = double(0)) {
  #     returnType(double(2))
  #     nrow_d <- dim(d)[1]
  #     ncol_d <- dim(d)[2]
  #     result <- matrix(0, nrow_d, ncol_d)  # Initialize result matrix
  #     for (i in 1:nrow_d) {
  #       for (j in 1:ncol_d) {
  #           if (d[i, j] == 0) {
  #             result[i, j] = 0.0
  #         } else {
  #             part <- (2^(1 - nu)) / gamma(nu)  * (sqrt(2 * nu) * d[i, j] / range)^nu
  #             result[i,j] = sigma.sq * part * besselK(d[i, j] / range, nu)
  #         }
  #       }
  #     }  
  #     return(result)
  #   }
  # )
  
  # Define Nimble model code
  model_code <- nimbleCode({
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- inprod(X[i, 1:P], beta[]) + spatial_effect[i] + year_effect[year[i]]
      y_pred[i] ~ dnorm(mu[i], tau)
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
      
    tau ~ dexp(1)
    sigma.sq ~ dexp(1)
    phi ~ dunif(5 / 1, 5 / 0.1)
    range ~ dunif(0.01, 1)
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
    tau = 1, 
    spatial_effect = rep(0, N), 
    sigma_year = 1,
    year_effect = rep(0, nyears),
    
    sigma.sq = 10, 
    range = 0.2,
    phi = 10,
    nu = 1.5
  )
  
  ### Change based on semi-variogram fitting
  constants <- list(N = N, N_knots = N_knots, N_unique = N_unique, P = ncol(X), Mu = rep(0, N), nyears = nyears, 
                    dist_matrix = dist_matrix, knots_dist_matrix = knots_dist_matrix, 
                    data_knots_dist_matrix = data_knots_dist_matrix,
                    year = year_index, I = diag(1e-6,N), cov_type = cov_type)  # Add 'year' to constants
  monitors <- switch(
    cov_type,
    "Mat" = c("beta", "range", "sigma.sq", "nu", "spatial_effect", "year_effect"),
    "Gau" = c("beta", "tau", "sigma.sq", "phi", "spatial_effect", "year_effect"),
    "Exp" = c("beta", "tau", "sigma.sq", "phi", "spatial_effect", "year_effect"),
    "Sph" = c("beta", "range", "sigma.sq", "spatial_effect", "year_effect"),
    stop("Unknown cov_type!")
  )
  
  ni <- 2000   # nb iterations 
  nt <- 1       # thinning interval
  nb <- 500    # nb iterations as burn-in 
  nc <- 2         # nb chains
  
  
  #################### Run/Save the model ####################
  
  my_model <- nimbleModel(code = model_code, data = data, inits = inits, constants = constants)
  compiled_model <- compileNimble(my_model)
  
  conf <- configureMCMC(my_model, enableWAIC = TRUE)
  conf$removeSampler(paste0("beta[", 1:4, "]"))
  conf$addSampler(target = paste0("beta[", 1:4, "]"), type = "RW_block")
  conf$removeSampler(c("sigma.sq", "range", "tau", "nu"))
  conf$addSampler("sigma.sq", type = "slice")
  conf$addSampler("range", type = "slice")
  conf$addSampler("tau", type = "slice")
  conf$addSampler("nu", type = "slice")
  
  conf$addMonitors(monitors)
  conf$addMonitors2("y_pred")
  
  mcmc <- buildMCMC(conf)
  compiled_mcmc <- compileNimble(mcmc, project = compiled_model)
  modelMCMCrun <- runMCMC(compiled_mcmc,
                          niter = ni, nburnin = nb, thin = nt, nchains = nc,
                          samplesAsCodaMCMC = TRUE, summary = TRUE, WAIC = TRUE)
  
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
  
  samples_matrix <- as.matrix(modelMCMCrun$samples2)
  y_pred_cols <- grep("y_pred\\[", colnames(samples_matrix))
  y_pred_mean <- colMeans(samples_matrix[, y_pred_cols])
  rmse <- sqrt(mean((y - y_pred_mean)^2))
  
  modelMCMCrun <- list(
    samples = modelMCMCrun$samples,
    summary = modelMCMCrun$summary,
    WAIC = modelMCMCrun$WAIC,
    rmse = rmse
  )
  write_rds(modelMCMCrun, file = model_path, compress = "gz")
}


