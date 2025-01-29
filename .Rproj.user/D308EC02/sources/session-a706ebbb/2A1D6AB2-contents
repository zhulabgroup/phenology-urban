library(nimble)
library(sp)
library(sf)
library(dplyr)
library(stringr)
library(lme4)
library(coda)
library(car)
library(gstat)

# Read file paths from environment variables
r_input_dir <- Sys.getenv("R_INPUT_DIR")
r_output_dir <- Sys.getenv("R_OUTPUT_DIR")
city <- Sys.getenv("CITY")
phe_type <- Sys.getenv("PHE_TYPE")

# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) < 3) {
#   stop("5 arguments needed: r_output_file, city, phe_type")
# }
# 
# r_input_dir <- args[1]
# r_output_dir <- args[2]
# city <- args[3]
# phe_type <- args[4]

# r_input_dir <- "~/phenology-urban/data/proc/urban/"
# r_output_dir <- "~/phenology-urban/data/proc/urban/sp_model/"
# city <- "NY"
# phe_type <- "SOS"

#################### Read data ####################

tavg_allyear_var <- readRDS(paste0(r_input_dir, "tavg_allyear_var.rds"))

if (phe_type == "SOS") {
  phe_metric <- "doy"
  direction <- "up"
  thres <- 0.5
} else if (phe_type == "EOS") {
  phe_metric <- "doy"
  direction <- "down"
  thres <- 0.5
} else if (phe_type == "GreenUp") {
  phe_metric <- "doy_diff"
  direction <- "up"
  thres <- 0.2
  tavg_allyear_var$doy_diff <- as.numeric(tavg_allyear_var$doy_diff)
} else {
  phe_metric <- "doy_diff"
  direction <- "down"
  thres <- 0.2
  tavg_allyear_var$doy_diff <- as.numeric(tavg_allyear_var$doy_diff)
}


#################### Find the optimal cov model ####################

unique_genus <- unique(tavg_allyear_var$genus)

for (genus in unique_genus) {
  tavg_phe <- tavg_allyear_var %>%
    filter(direction == !!direction & thres == !!thres & genus == !!genus) %>%
    na.omit()
  print(nrow(tavg_phe))
  if (nrow(tavg_phe) == 0) next
  
  OLS_formula <- as.formula(paste0(phe_metric, " ~ AvgTemp_mean + Sum_mm_sum + srad_mean + (1|year)"))
  random_OLS <- lmer(OLS_formula, data = tavg_phe)
  tavg_phe$residual_random <- tavg_phe[[phe_metric]] - fitted(random_OLS)
  coordinates(tavg_phe) <- ~ lon + lat
  vgm.emt <- variogram(residual_random ~ 1, tavg_phe, cloud = FALSE)
  vgm.las.init <- vgm(model = c("Exp", "Sph", "Gau", "Mat"))
  vgm_exn_fit <- fit.variogram(vgm.emt, vgm.las.init)
  cov_type <- as.character(vgm_exn_fit$model)
  
  #################### Run the spatial model ####################
  
  coords <- coordinates(tavg_phe)
  tavg_phe_df <- as.data.frame(tavg_phe)
  sp_formula <- as.formula(paste0(phe_metric, " ~ AvgTemp_mean + Sum_mm_sum + srad_mean"))
  
  # Define covariance functions
  
  spherical_covariance <- nimbleFunction(
    run = function(d = double(0), range = double(0), sigma.sq = double(0)) {
      returnType(double(0)) 
      if (d <= range) {
        return (sigma.sq * (1 - (3 * d / (2 * range)) + (d^3 / (2 * range^3))))
      } else {
        return(0.0)
      }
    }
  )
  
  matern_covariance <- nimbleFunction(
    run = function(d = double(0), range = double(0), sigma.sq = double(0), nu = double(0)) {
      returnType(double(0))
        if (d == 0) {
          return(sigma.sq)
        } else {
          part <- (2^(1 - nu)) / gamma(nu)  * (sqrt(2 * nu) * d / range)^nu
          return(sigma.sq * part * besselK(d / range, nu))
        }
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
    
    if (cov_type == "Mat") {
      for (m in 1:N) {
        for (n in 1:N) {
          cov_matrix[m, n] <- matern_covariance(dist_matrix[m, n], range, sigma.sq, nu)
        }
      }
    } else if (cov_type == "Gau") {
      cov_matrix[1:N, 1:N] <- sigma.sq * exp(-phi * dist_matrix[1:N, 1:N]^2)
    } else if (cov_type == "Exp") {
      cov_matrix[1:N, 1:N] <- sigma.sq * exp(-phi * dist_matrix[1:N, 1:N])
    } else if (cov_type == "Sph") {
      for (m in 1:N) {
        for (n in 1:N) {
          cov_matrix[m, n] <- spherical_covariance(dist_matrix[m, n], range, sigma.sq)
        }
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
  
  N <- nrow(tavg_phe_df)
  nyears <- length(unique(tavg_phe_df$year))
  X <- model.matrix(sp_formula, data = tavg_phe_df)
  y <- tavg_phe_df[[all.vars(sp_formula)[1]]]
  dist_matrix <- as.matrix(dist(coords))
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
  constants <- list(N = N, P = ncol(X), Mu = rep(0, N), nyears = nyears, dist_matrix = dist_matrix,
                    year = year_index, I = diag(1,N), cov_type = cov_type)  # Add 'year' to constants
  monitors <- switch(
    cov_type,
    "Mat" = c("beta", "tau.year", "range", "sigma.sq", "nu", "spatial_effect", "year_effect"),
    "Gau" = c("beta", "tau", "sigma.sq", "phi", "tau.year", "spatial_effect", "year_effect"),
    "Exp" = c("beta", "tau", "sigma.sq", "phi", "tau.year", "spatial_effect", "year_effect"),
    "Sph" = c("beta", "tau.year", "range", "sigma.sq", "spatial_effect", "year_effect"),
    stop("Unknown cov_type!")
  )
  
  ni <- 20000   # nb iterations 
  nt <- 1       # thinning interval
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
  
  model_path <- paste0(r_output_dir, city, "_", genus, "_tavg_", phe_type,  ".rds")
  saveRDS(modelMCMCrun, file = model_path)
}


