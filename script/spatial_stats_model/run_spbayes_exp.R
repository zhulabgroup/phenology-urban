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
library(spBayes)
library(bridgesampling)
library(fields)
select <- dplyr::select


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
  out_dir <- paste0(r_output_dir, city, "/", taxon_col, "/")
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
  dist_mat <- as.matrix(dist(coords))
  max.d <- max(iDist(coords))
  N_knots <- 50
  N_unique <- nrow(unique(coords))
  if (N_unique < 30) {
    message(paste("Not enough unique coordinates for taxon:", taxon, "with", N_unique, "unique points. Skipping..."))
    next  # Skip the current iteration of the loop
  }
  if (N_knots <= N_unique){
    use_knots <- TRUE
    knots <- kmeans(coords, centers = N_knots)$centers
    dist_mat_nm <- as.matrix(rdist(coords, knots))
    dist_mat_mm <- as.matrix(rdist(knots, knots))
    
  } else {
    use_knots <- FALSE
    knots <- NULL
    dist_mat_nm <- dist_mat
    dist_mat_mm <- dist_mat
  }
  tavg_phe_df <- as.data.frame(tavg_phe)
  tavg_phe_df$year <- as.factor(tavg_phe_df$year) 
  
  sp_formula <- as.formula(paste0(
    phe_metric, " ~ ", 
    temp_col, " + ", 
    prcp_col, " + ", 
    srad_col, " + factor(year)"
  ))
  ols_model <- lm(sp_formula, data = tavg_phe_df)
  X <- model.matrix(sp_formula, data = tavg_phe_df)
  y <- tavg_phe_df[[phe_metric]]
  p <- ncol(X)
  
  tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
  
  priors <- list(
    "beta.Norm" = list(rep(0, p), diag(100, p)),
    "phi.Unif" = c(3/1, 3/0.1),  
    "sigma.sq.IG" = c(2, 2),
    "tau.sq.IG" = c(2, 0.1)
  )
  n.samples = 10000
  burn.in = 5000
  
  run_spLM_chains <- function(
    formula, data, coords, knots,
    starting_list, tuning, priors,
    cov.model = "matern", n.samples = n.samples,
    burn.in = burn.in,
    param_pattern = "temp_",  # 提取的参数匹配模式
    n.report = 2000, verbose = TRUE, use_knots
  ) {
    
    n.chains <- length(starting_list)
    if (is.null(burn.in)) burn.in <- floor(0.5 * n.samples)
    chain_results <- vector("list", n.chains)
    chain_recover <- vector("list", n.chains)
    
    for (i in seq_len(n.chains)) {
      if (verbose) message("Running ", i, " chain")
      
      args_list <- list(
        formula = formula,
        data = data,
        coords = coords,
        starting = starting_list[[i]],
        tuning = tuning,
        priors = priors,
        cov.model = cov.model,
        n.samples = n.samples,
        n.report = n.report
      )
      if (use_knots) args_list$knots <- knots
      
      # 运行 spLM
      sp_model <- do.call(spLM, args_list)
      
      chain_object <- spRecover(sp_model, start = burn.in, verbose = FALSE)
      chain_recover[[i]] <- chain_object
      # param_cols <- grep(param_pattern, colnames(chain_recover[[i]]$p.beta.recover.samples))
      chain <- chain_object$p.beta.recover.samples
      # varnames(chain) <- "beta_temp"
      chain_results[[i]] <- mcmc(chain)
    }
    beta_mcmc_list <- mcmc.list(chain_results)
    
    return(list(
      beta_mcmc_list = beta_mcmc_list,
      chains_recover = chain_recover
    ))
  }
  
  starting1 <- list("beta" = as.numeric(ols_model$coefficients),
                    "sigma.sq" = 1, "tau.sq" = 1, "phi" = 3/(0.5*max.d))
  starting2 <- list("beta" = rep(0, p), "sigma.sq" = 1, "tau.sq" = 1, "phi" = 3/0.5)
  starting_list <- list(chain1 = starting1, chain2 = starting2)
  res <- run_spLM_chains(
    formula = sp_formula,
    data = tavg_phe_df,
    coords = coords,
    knots = knots,
    starting_list = starting_list,
    tuning = tuning,
    priors = priors,
    cov.model = "exponential",
    n.samples = n.samples,
    burn.in = burn.in,
    param_pattern = "temp_",
    use_knots = use_knots
  )
  beta2_mcmc <- lapply(res$beta_mcmc_list, function(chain) {
    temp_cols <- grep("^temp_", colnames(chain), value = TRUE)  # 匹配以 temp_ 开头的列
    as.mcmc(chain[, temp_cols, drop = FALSE])
  }) |> as.mcmc.list()
  converge <- gelman.diag(beta2_mcmc, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
                          multivariate=TRUE)

  # library(ggmcmc)
  # beta_label <- data.frame(Parameter = "temp_summer", Label = "AvgTemp")
  # param_data <- ggs(res$beta_mcmc_list, par_labels = beta_label) %>%
  #   filter(Parameter == "AvgTemp")
  # beta_mean <- mean(param_data$value)
  # param_data %>%
  #   ggs_density() + # traceplot
  #   theme_bw() +
  #   ggtitle(paste("taxon:", taxon)) +
  #   theme(plot.title = element_text(hjust = 0.5)) +
  #   annotate("text", x = beta_mean, y = 0, label = paste0("Mean: ", round(beta_mean, 2)),vjust = -1.5, color = "red", size = 4)

  log_posterior_spLM <- function(params, data) {
    y <- data$y
    X <- data$X
    priors <- data$priors
    dist_mat_nm <- data$dist_mat_nm
    dist_mat_mm <- data$dist_mat_mm
    
    beta_names <- colnames(X)
    beta <- as.numeric(params[beta_names])
    sigma.sq <- params["sigma.sq"]
    tau.sq   <- params["tau.sq"]
    phi      <- params["phi"]
    
    if (any(is.na(c(beta, sigma.sq, tau.sq, phi))) || sigma.sq <= 0 || tau.sq <= 0) {
      return(-Inf)
    }
    n <- length(y)
    m <- nrow(dist_mat_mm)
    Rnm <- exp(-phi * dist_mat_nm)       # n x m
    Rmm <- exp(-phi * dist_mat_mm) + 1e-6 * diag(m)
    Lmm <- chol(Rmm)
    Q <- (sigma.sq^(0.5)) * forwardsolve(t(Lmm), t(Rnm))  # m x n
    mu <- X %*% beta
    
    logdet_Sigma <- n * log(tau.sq) +
      determinant(diag(m) + (1/tau.sq) * Q %*% t(Q), logarithm=TRUE)$modulus
    
    resid <- y - mu
    alpha <- resid / tau.sq -
      (1 / tau.sq^2) * t(Q) %*% solve(diag(m) + (1/tau.sq) * Q %*% t(Q), Q %*% resid)
    
    quadform <- sum(resid * alpha)
    loglik <- -0.5 * (quadform + logdet_Sigma + n * log(2 * pi))
    
    logprior <- mvtnorm::dmvnorm(beta,
                                 mean = priors$beta.Norm[[1]],
                                 sigma = priors$beta.Norm[[2]],
                                 log = TRUE)
    
    # 先验：sigma.sq 和 tau.sq 逆伽马
    logprior <- logprior + invgamma::dinvgamma(sigma.sq,
                                               shape = priors$sigma.sq.IG[1],
                                               scale = priors$sigma.sq.IG[2],
                                               log = TRUE)
    logprior <- logprior + invgamma::dinvgamma(tau.sq,
                                               shape = priors$tau.sq.IG[1],
                                               scale = priors$tau.sq.IG[2],
                                               log = TRUE)
    
    # 先验：phi uniform
    if (phi < priors$phi.Unif[1] || phi > priors$phi.Unif[2]) {
      return(-Inf)
    } else {
      logprior <- logprior + log(1 / (priors$phi.Unif[2] - priors$phi.Unif[1]))
    }
    
    return(loglik + logprior)
  }
  
  calc_marginal_likelihood <- function(chains_recover, log_posterior_fun, y, X, priors, dist_mat_nm, dist_mat_mm) {
    message("Calculating marginal likelihood...")
    sample_list <- lapply(chains_recover, function(chain) {
      cbind(chain$p.beta.recover.samples, chain$p.theta.recover.samples)
    })
    samples_all <- do.call(rbind, sample_list)
    samples <- as.matrix(samples_all)
    param_names <- colnames(samples)
    
    lb <- setNames(rep(-Inf, length(param_names)), param_names)
    ub <- setNames(rep( Inf, length(param_names)), param_names)
    if ("sigma.sq" %in% param_names) lb["sigma.sq"] <- 0
    if ("tau.sq"   %in% param_names) lb["tau.sq"]   <- 0
    if ("phi"      %in% param_names) {
      lb["phi"] <- priors$phi.Unif[1]
      ub["phi"] <- priors$phi.Unif[2]
    }
    
    bridge_result <- bridgesampling::bridge_sampler(
      samples = samples,
      log_posterior = log_posterior_fun,
      data = list(y = y, X = X, priors = priors, dist_mat_nm = dist_mat_nm, dist_mat_mm = dist_mat_mm),
      lb = lb,
      ub = ub,
      param_types = rep("real", ncol(samples)),
      repetitions = 10,
      method = "warp3",
      cores = 1,
      use_neff = TRUE,
      silent = FALSE,
      verbose = TRUE,
      varlist = ls(),
      envir = environment()
    )
    return((bridge_result$logml))
  }

  res_logml <- calc_marginal_likelihood(chains_recover = res$chains_recover, log_posterior_fun = log_posterior_spLM, 
                                        y = y, X = X, priors = priors, dist_mat_nm = dist_mat_nm, dist_mat_mm = dist_mat_mm)
  
  calc_WAIC_spLM <- function(chains_recover, formula, data) {
    message("Calculating WAIC...")
    X <- model.matrix(formula, data)
    y <- data[[all.vars(formula)[1]]]
    n <- length(y)
    
    loglik_list <- lapply(chains_recover, function(chain) {
      beta_samples <- chain$p.beta.recover.samples        # S x p
      tau_samples  <- chain$p.theta.recover.samples[, "tau.sq"]  # S x 1
      w_samples    <- t(chain$p.w.recover.samples)         # 转为 S x n
      S <- nrow(beta_samples)
      
      loglik_samples <- matrix(NA, nrow = S, ncol = n)
      for (s in 1:S) {
        mu_s <- X %*% beta_samples[s, ] + w_samples[s, ]
        
        loglik_samples[s, ] <- dnorm(y, mean = mu_s, sd = sqrt(tau_samples[s]), log = TRUE)
      }
      loglik_samples
    })
    
    loglik_samples_all <- do.call(rbind, loglik_list)
    lppd <- sum(log(colMeans(exp(loglik_samples_all))))
    p_waic <- sum(apply(loglik_samples_all, 2, var))
    waic <- -2 * (lppd - p_waic)
    return(waic)
  }
  
  res_waic <- calc_WAIC_spLM(chains_recover = res$chains_recover, formula = sp_formula, data = tavg_phe_df)
  
  modelMCMCrun <- list(
    samples = res$beta_mcmc_list,
    WAIC = list(WAIC = res_waic),
    converge = converge$psrf[1],
    mean_logml = mean(res_logml)
  )
  write_rds(modelMCMCrun, file = model_path, compress = "gz")
}


