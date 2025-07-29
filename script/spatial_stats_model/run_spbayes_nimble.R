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
select <- dplyr::select

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
  N <- nrow(coords)
  N_knots <- 50
  N_unique <- nrow(unique(coords))
  if (N_knots <= N_unique){
    knots <- kmeans(coords, centers = N_knots)$centers
    knots_dist_matrix <- as.matrix(dist(knots))
    data_knots_dist_matrix <- as.matrix(dist(rbind(coords, knots)))[1:N, (N + 1):(N + N_knots)]
    
  } else {
    knots <- coords
    knots_dist_matrix <- dist_matrix
    data_knots_dist_matrix <- dist_matrix
  }
    
  pre_exp_knots <- exp(-knots_dist_matrix^2)
  pre_exp_data_knots <- exp(-data_knots_dist_matrix^2)
  
  # 2. Design matrix
  tavg_phe_df <- as.data.frame(tavg_phe)
  sp_formula <- as.formula(paste0(phe_metric, " ~ ", temp_col, "+", prcp_col, "+", srad_col))
  X <- model.matrix(sp_formula, data = tavg_phe_df)
  y <- tavg_phe_df[[all.vars(sp_formula)[1]]]
  nyears <- length(unique(tavg_phe_df$year))
  year_index <- as.integer(factor(tavg_phe_df$year))

  update_cov_matrices <- nimbleFunction(
    setup = function(pre_exp_knots, pre_exp_data_knots, N_knots, N) {
      kd <- pre_exp_knots
      dkd <- pre_exp_data_knots
      m <- N_knots
      n <- N
      
      K_store <- matrix(0, m, m)
      B_store <- matrix(0, n, m)
      chol_store <- matrix(0, m, m)
    },
    run = function(phi = double(0), sigma_sq = double(0)) {
      returnType(void())
      for (i in 1:m) {
        for (j in 1:m) {
          K_store[i, j] <<- sigma_sq * (kd[i, j] ^ phi)
        }
        K_store[i, i] <<- K_store[i, i] + 1e-6
      }
      chol_store[1:m, 1:m] <<- chol(K_store[1:m, 1:m])
      
      for (i in 1:n) {
        for (j in 1:m) {
          B_store[i, j] <<- sigma_sq * (dkd[i, j] ^ phi)
        }
      }
      return()
    },
    methods = list(
      getB = function() {
        returnType(double(2))
        return(B_store)
      },
      getL = function() {
        returnType(double(2))
        return(chol_store)
      }
    )
  )
  
  update_cov <- update_cov_matrices(
    pre_exp_knots = pre_exp_knots,
    pre_exp_data_knots = pre_exp_data_knots,
    N_knots = N_knots,
    N = N
  )
  update_cov$run(phi = 5, sigma_sq = 1)
  B_init <- update_cov$getB()
  L_init <- update_cov$getL()
  
  #### 3. NIMBLE 模型 ####
  model_code <- nimbleCode({
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- inprod(X[i, 1:P], beta[]) + spatial_effect[i] + year_effect[year[i]]
    }
    for (q in 1:nyears) {
      z_year[q] ~ dnorm(0, 1)
      year_effect[q] <- z_year[q] * sigma_year
    }
    
    # 低秩空间效应：eta ~ N(0, I)
    for (k in 1:N_knots) {
      eta[k] ~ dnorm(0, 1)
    }
    # 空间效应： B %*% (L_eta %*% eta)
    for (i in 1:N) {
      spatial_effect[i] <- inprod(B[i, 1:N_knots], tmp_eta[1:N_knots])
    }
    tmp_eta[1:N_knots] <- inprod(L[1:N_knots, 1:N_knots], eta[1:N_knots])
    
    sigma_year ~ dexp(1)
    beta[1] ~ dnorm(180, sd = 40)
    beta[2] ~ dnorm(0, sd = 5)
    beta[3] ~ dnorm(0, sd = 1)
    beta[4] ~ dnorm(0, sd = 1)
    
    tau ~ dexp(1)
    sigma.sq ~ dexp(1)
    phi ~ dunif(5 / 1, 5 / 0.1)
  })
  
  data <- list(y = y, X = X)
  ### Change based on semi-variogram fitting
  constants <- list(N = N, N_knots = N_knots, P = ncol(X), nyears = nyears, 
                    B = B_init, L = L_init,
                    year = year_index)  # Add 'year' to constants
  ols_model <- lm(sp_formula, data = tavg_phe_df)
  inits <- list(
    beta = as.numeric(ols_model$coefficients),
    tau = 1,
    sigma_year = 1,
    z_year = rep(0, nyears),
    eta = rep(0, N_knots),
    sigma.sq = 1,
    phi = 10
  )
  
  
  #### 4. 定义 block sampler ####
  block_sampler_phi_sigma <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      phi_name <- control$phi_name
      sigma_name <- control$sigma_name
      update_cov <- control$update_cov
    },
    run = function() {
      current_phi <- model[[phi_name]]
      current_sigma <- model[[sigma_name]]
      
      # RW proposal
      proposed_phi <- rnorm(1, mean = current_phi, sd = 0.1)
      proposed_sigma <- rnorm(1, mean = current_sigma, sd = 0.1)
      
      if (proposed_phi <= 0 | proposed_sigma <= 0) return()
      
      logProb_current <- model$getLogProb(phi_name) + model$getLogProb(sigma_name)
      
      model[[phi_name]] <<- proposed_phi
      model[[sigma_name]] <<- proposed_sigma
      update_cov$run(proposed_phi, proposed_sigma)
      model$B[,] <<- update_cov$getB()
      model$L[,] <<- update_cov$getL()
      
      logProb_proposed <- model$calculate(phi_name) + model$calculate(sigma_name)
      log_MH_ratio <- logProb_proposed - logProb_current
      
      accept <- decide(log_MH_ratio)
      if (!accept) {
        model[[phi_name]] <<- current_phi
        model[[sigma_name]] <<- current_sigma
        update_cov$run(current_phi, current_sigma)
        model$B[,] <<- update_cov$getB()
        model$L[,] <<- update_cov$getL()
        model$calculate(phi_name)
        model$calculate(sigma_name)
      }
    },
    methods = list(
      reset = function() {
        return()
      }
    )
  )
  

  monitors <- c("beta", "tau", "sigma.sq", "phi", "eta", "year_effect")
  
  ni <- 3000   # nb iterations 
  nt <- 10       # thinning interval
  nb <- 2500    # nb iterations as burn-in 
  nc <- 2         # nb chains
  
  
  #################### Run/Save the model ####################
  
  my_model <- nimbleModel(code = model_code, data = data, inits = inits, constants = constants, buildDerivs = TRUE)
  
  compiled_model <- compileNimble(my_model)
  conf <- configureMCMC(my_model, enableWAIC = TRUE)

  conf$removeSampler(c("phi", "sigma.sq"))
  conf$addSampler(
    target = c("phi", "sigma.sq"),
    type = "RW_block",
    control = list(propCov = diag(c(0.05, 0.05)))  # 可调步长
  )
  # conf$addSampler(target = c("phi", "sigma.sq"),
  #                 type = block_sampler_phi_sigma,
  #                 control = list(
  #                   phi_name = "phi",
  #                   sigma_name = "sigma.sq",
  #                   update_cov = update_cov
  #                 ))
  
  
  conf$removeSampler(paste0("beta[", 1:4, "]"))
  conf$addSampler(target = "beta[1]", type = "slice")
  conf$addSampler(target = paste0("beta[", 2:4, "]"), type = "RW_block")
  
  conf$addMonitors(monitors)
  
  mcmc <- buildMCMC(conf)
  compiled_mcmc <- compileNimble(mcmc, project = compiled_model)
  system.time({
    modelMCMCrun <- runMCMC(compiled_mcmc,
                            niter = ni, nburnin = nb, thin = nt, nchains = nc,
                            samplesAsCodaMCMC = TRUE, summary = TRUE, WAIC = TRUE)
  })
  
  beta2_mcmc <- lapply(modelMCMCrun$samples, function(chain) {
    as.mcmc(chain[, "beta[2]", drop = FALSE])
  }) |> as.mcmc.list()
  converge <- gelman.diag(beta2_mcmc, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
                          multivariate=TRUE)
  library(ggmcmc)
  beta_label <- data.frame(Parameter = "beta[2]", Label = "AvgTemp")
  param_data <- ggs(modelMCMCrun$samples, par_labels = beta_label) %>%
    filter(Parameter == "AvgTemp")
  beta_mean <- mean(param_data$value)
  param_data %>%
    ggs_density() + # traceplot
    theme_bw() +
    ggtitle(paste("taxon:", taxon)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = beta_mean, y = 0, label = paste0("Mean: ", round(beta_mean, 2)),vjust = -1.5, color = "red", size = 4)
# 
#   
#   log_posterior_nimble <- function(samples_row, data) {
#     constants <- data$constants
#     param_names <- names(samples_row)
#     if (is.null(param_names)) stop("samples_row does not have named parameters")
#     
#     
#     beta <- samples_row[grep("^beta\\[", param_names)]
#     tau <- samples_row["tau"]
#     sigma_year <- samples_row["sigma_year"]
#     sigma_sq <- samples_row["sigma.sq"]
#     phi <- samples_row["phi"]
#     spatial_effect <- samples_row[grep("^spatial_effect\\[", param_names)]
#     z_year <- samples_row[grep("^z_year\\[", param_names)]
#     
#     if (tau <= 0 || sigma_sq <= 0 || phi <= 0) return(-Inf)
#     
#     log_prior <- sum(
#       dnorm(beta[1], mean = 180, sd = 40, log = TRUE),
#       dnorm(beta[2], 0, 5, log = TRUE),
#       dnorm(beta[3], 0, 1, log = TRUE),
#       dnorm(beta[4], 0, 1, log = TRUE),
#       if (length(beta) > 4) sum(dnorm(beta[5:length(beta)], 0, 10, log = TRUE)) else 0,
#       dexp(sigma_year, rate = 1, log = TRUE),
#       dexp(sigma_sq, rate = 1, log = TRUE),
#       dunif(phi, 5/1, 5/0.1, log = TRUE),
#       dgamma(tau, shape = 0.001, rate = 0.001, log = TRUE),
#       sum(dnorm(z_year, 0, 1, log = TRUE))
#     )
#     
#     N <- constants$N
#     N_knots <- data$N_knots
#     N_unique <- nrow(unique(data$coords))
#     spatial_mean <- rep(0, N)
#     
#     log_mvn_chol <- function(x, mean, Sigma) {
#       L <- tryCatch(chol(Sigma), error = function(e) return(NULL))
#       if (is.null(L)) return(-Inf)
#       diff <- x - mean
#       logdet <- sum(log(diag(L)))
#       quad <- backsolve(L, diff, transpose = TRUE)
#       quad <- sum(quad^2)
#       return(-0.5 * length(x) * log(2 * pi) - logdet - 0.5 * quad)
#     }
#     
#     if (N_knots <= N_unique) {
#       # --- Knots-based GP approximation ---
#       K_kk <- sigma_sq * exp(-phi * data$knots_dist_matrix_sq) + diag(1e-6, N_knots)
#       K_nk <- sigma_sq * exp(-phi * data$data_knots_dist_matrix_sq)
#       
#       # Cholesky
#       L_kk <- tryCatch(chol(K_kk), error = function(e) return(NULL))
#       if (is.null(L_kk)) return(-Inf)
#       K_kk_inv <- chol2inv(L_kk)
#       cov_approx <- K_nk %*% K_kk_inv %*% t(K_nk) + diag(1e-5, N)
#       
#       L_cov <- tryCatch(chol(cov_approx), error = function(e) return(NULL))
#       if (is.null(L_cov)) return(-Inf)
#       
#       log_prior <- log_prior + log_mvn_chol(spatial_effect, spatial_mean, cov_approx)
#     } else {
#       cov_matrix <- sigma_sq * exp(-phi * data$dist_matrix_sq) + diag(1e-6, N)
#       log_prior <- log_prior + log_mvn_chol(spatial_effect, spatial_mean, cov_matrix)
#     }
#     
#     # log-likelihood
#     X <- data$X
#     y <- data$y
#     year <- constants$year
#     year_effect <- z_year * sigma_year
#     
#     mu <- as.vector(X %*% beta) + year_effect[year] + spatial_effect
#     log_lik <- sum(dnorm(y, mean = mu, sd = 1 / sqrt(max(tau, 1e-8)), log = TRUE))
#     
#     return(log_lik + log_prior)
#   }
#   
#   if (N_knots <= N_unique) {
#     knots <- kmeans(coords, centers = N_knots)$centers
#     knots_dist_matrix <- as.matrix(dist(knots))
#     data_knots_dist_matrix <- as.matrix(dist(rbind(coords, knots)))[1:N, (N + 1):(N + N_knots)]
#   } else {
#     knots <- coords
#     knots_dist_matrix <- dist_matrix
#     data_knots_dist_matrix <- dist_matrix
#   }
#   
#   
#   data_bridge <- list(
#     y = y,
#     X = X,
#     coords = coords,
#     dist_matrix_sq = dist_matrix^2,
#     knots = knots,
#     knots_dist_matrix_sq = knots_dist_matrix^2,
#     data_knots_dist_matrix_sq = data_knots_dist_matrix^2,
#     N_knots = N_knots,
#     constants = constants
#   )
#   
#   posterior_samples <- as.matrix(modelMCMCrun$samples)
#   param_names <- colnames(posterior_samples)
#   
#   model_bridge <- bridge_sampler(
#     samples = posterior_samples,
#     log_posterior = log_posterior_nimble,
#     data = data_bridge,
#     lb = setNames(rep(-Inf, length(param_names)), param_names),
#     ub = setNames(rep( Inf, length(param_names)), param_names),
#     param_types = rep("real", ncol(posterior_samples)),
#     repetitions = 10,
#     method = "warp3",
#     cores = 1,
#     use_neff = TRUE,
#     silent = FALSE,
#     verbose = TRUE,
#     varlist = ls(),
#     envir = environment()
#   )
#   mean_logml <- mean(model_bridge$logml)
#   
  
  modelMCMCrun <- list(
    samples = modelMCMCrun$samples,
    summary = modelMCMCrun$summary,
    WAIC = modelMCMCrun$WAIC,
    converge = converge$psrf[1],
    # mean_logml = mean_logml
  )
  write_rds(modelMCMCrun, file = model_path, compress = "gz")
}


