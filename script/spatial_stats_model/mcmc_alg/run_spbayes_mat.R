# taxon = "Betula nigra"

run_mcmc_spatial_mat <- function(model_input_all, r_output_dir, phe_type, city, taxon_col, 
                                 n_cores = 6) {
  taxon_list <- unique(model_input_all[[taxon_col]])
  boundary <- switch(
    city,
    "NY" = st_read("~/urban-cooling/data/raw/WU/NY/Untitled/nybb_dissolved.shp"),
    "DT" = st_read("~/urban-cooling/data/raw/WU/DT/City_of_Detroit_Boundary/City_of_Detroit_Boundary.shp"),
    "ST" = st_read("~/urban-cooling/data/raw/WU/ST/seattle_boundary/CityUGA.shp") %>%
      filter(CITY_NM == "Seattle"),
    "HT" = st_read("~/urban-cooling/data/raw/WU/HT/HGAC_City_Boundaries/HGAC_City_Boundaries.shp") %>%
      filter(NAME == "Houston"),  # change to some distance
    "TP" = st_read("~/urban-cooling/data/raw/WU/TP/Municipal_Boundary/Municipal_Boundary.shp") %>%
      filter(NAME == "Tampa"),
    "DV" = st_read("~/urban-cooling/data/raw/WU/DV/DRCOG_Municipalities/DRCOG_Municipalities.shp") %>%
      filter(city == "Denver"),
    "AT" = st_read("~/urban-cooling/data/raw/WU/AT/TxDOT_City_Boundaries/Cities.shp") %>%
      filter(CITY_NM == "Austin"),
    stop("Invalid city abbreviation.")
  )
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
  
  cl <- parallel::makeCluster(n_cores, outfile = "")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, { library(sf); options(sf_use_s2 = FALSE) })
  parallel::clusterExport(cl, c("boundary"), envir = environment())
  
  `%dopar%` <- foreach::`%dopar%`
  res <- foreach(taxon = taxon_list, .packages = c(
    "spBayes", "dplyr", "data.table","here","lubridate", "sf",
    "tidyverse", "coda", "terra", "fields", "ggplot2", "RANN"
  ),
  .inorder = FALSE, .errorhandling = "pass") %dopar% {
    cat(sprintf("Worker PID %d running taxon: %s\n", Sys.getpid(), taxon))
    flush.console()
    
    ############################ Set input ############################
    taxon_clean <- gsub(" ", "_", taxon)
    out_dir <- paste0(r_output_dir, "/spatial_mat/",city,"/", taxon_col,"/")
    # if (!dir.exists(out_dir)) {
    #   dir.create(out_dir, recursive = TRUE)
    # }
    model_path <- paste0(out_dir, taxon_clean, "_", phe_type, ".rds")
    
    if (file.exists(model_path)) {
      next  # Skip the current iteration of the loop
    }
    model_input_sf <- model_input_all |>
      filter(!!sym(taxon_col) == !!taxon) |>
      mutate(
        lon_jitter = lon,
        lat_jitter = lat
      )
    if (phe_type %in% c("SOS","EOS")) {
      phe_metric <- "doy"
    } else {
      phe_metric <- "doy_diff"
    }
    tavg_phe <- model_input_sf |> st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F) |>
      st_transform(crs = st_crs(boundary)) #jitter
    temp_col <- names(tavg_phe)[grepl("^temp", names(tavg_phe))]
    prcp_col <- names(tavg_phe)[grepl("^prcp", names(tavg_phe))]
    srad_col <- names(tavg_phe)[grepl("^srad", names(tavg_phe))]
    coords_ori <- st_coordinates(model_input_sf |> st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F) |>
                                   st_transform(crs = st_crs(boundary)))
    
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
    coords <- st_coordinates(tavg_phe)
    # dup_idx <- duplicated(coords_ori)
    # coords <- coords_ori
    # n_dup <- sum(dup_idx, na.rm = TRUE)
    # if (n_dup > 0) {
    #   coords[dup_idx, ] <- coords[dup_idx, ] +
    #     matrix(rnorm(n_dup * 2, mean = 0, sd = 1), ncol = 2)
    # }
    
    dist_mat <- as.matrix(dist(coords))
    max.d <- max(iDist(coords))
    N_knots <- 50
    N_unique <- nrow(unique(coords_ori))
    if (N_unique < 30) {
      message(paste("Not enough unique coordinates for taxon:", taxon, "with", N_unique, "unique points. Skipping..."))
      next  # Skip the current iteration of the loop
    }
    message(paste("Running spatial model for taxon:", taxon, "with", N_unique, "unique points."))
    bbox_sfc <- tavg_phe |>
      st_union() |>
      st_convex_hull()
    grid_poly_all <- st_make_grid(bbox_sfc, cellsize = 1000, what = "polygons")
    grid_poly_in <- st_intersection(grid_poly_all, bbox_sfc)
    grid_pts <- st_centroid(grid_poly_in)
    
    pred.coords <- st_coordinates(grid_pts)
    pred.covar <- matrix(0, nrow = nrow(pred.coords), ncol = p)
    colnames(pred.covar) <- colnames(X)
    pred.covar[, "(Intercept)"] <- 1  
    
    
    if (50 < N_unique){
      use_knots <- TRUE
      bbox_sfc2 <- tavg_phe |>
        st_union() |>
        st_concave_hull(ratio = 0.4)
      grid_knots <- st_make_grid(bbox_sfc2, n = c(15, 15), what = "centers")
      grid_knots_in <- st_intersection(grid_knots, bbox_sfc2)
      knots <- st_coordinates(grid_knots_in)
      tmp <- RANN::nn2(coords, knots, searchtype = 'radius', k = 100, radius = 20000)
      np  <- tmp[[1]]
      np[np > 0] <- 1
      coverage <- rowSums(np)
      ntot <- 50
      idx_nonzero <- which(coverage > 0)
      wt <- order(coverage[idx_nonzero], decreasing = TRUE)[1:ntot]
      sel <- idx_nonzero[wt]
      knots <- knots[sel[!is.na(sel)], ]
      dist_mat_nm <- as.matrix(rdist(coords, knots))
      dist_mat_mm <- as.matrix(rdist(knots, knots))
    } else {
      use_knots <- FALSE
      knots <- NULL
      dist_mat_nm <- dist_mat
      dist_mat_mm <- dist_mat
    }
    
    tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1, "nu"=0.1)
    ig_from_mean_cv <- function(m, cv = 1) {
      stopifnot(m > 0, cv > 0)
      shape <- 1 + 1/(cv^2)
      scale <- m * (shape - 1)
      c(shape = shape, scale = scale)
    }
    build_spbayes_ig_priors <- function(ols_model, tau_ratio = 0.1, sigma_cv = 1, tau_cv = 1) {
      sigma2_resid <- (summary(ols_model)$sigma)^2
      sigma_sq_ig <- ig_from_mean_cv(m = sigma2_resid, cv = sigma_cv)
      tau_sq_ig   <- ig_from_mean_cv(m = tau_ratio * sigma2_resid, cv = tau_cv)
      list(
        sigma.sq.IG = unname(sigma_sq_ig),
        tau.sq.IG   = unname(tau_sq_ig),
        resid_var   = sigma2_resid
      )
    }
    ig_priors <- build_spbayes_ig_priors(
      ols_model,
      tau_ratio = 0.1,  
      sigma_cv  = 1.0, 
      tau_cv    = 1.0
    )
    
    priors <- list(
      "beta.Norm" = list(rep(0, p), diag(100, p)),
      "phi.Unif" = c(3/(0.75*max.d), 3/(0.001*max.d)),  
      "sigma.sq.IG" = ig_priors$sigma.sq.IG,
      "tau.sq.IG" = ig_priors$tau.sq.IG,
      "nu.Unif" = c(0.5, 2.5)
    )
    
    n.samples = 10000
    burn.in = 5000
    starting1 <- list("beta" = as.numeric(ols_model$coefficients), 
                      "sigma.sq" = summary(ols_model)$sigma^2, "tau.sq" = 0.1*summary(ols_model)$sigma^2, 
                      "phi" = 3/(0.1*max.d), "nu" = 1.5)
    starting2 <- list("beta" = rep(0, p), 
                      "sigma.sq" = summary(ols_model)$sigma^2, "tau.sq" = 0.1*summary(ols_model)$sigma^2, 
                      "phi" = 3/(0.1*max.d), "nu" = 1.5)
    
    starting_list <- list(chain1 = starting1, chain2 = starting2)
    
    
    
    ############################ Run model ############################
    
    # safe_spRecover_beta <- function(sp_model, burn.in = 1, n.samples = NULL) {
    #   beta_list <- list()
    #   theta_list <- list()
    #   w_list <- list()
    #   used_idx <- integer(0)
    #   for (i in seq(from = burn.in, to = n.samples)) {
    #     tryCatch({
    #       co <- spRecover(sp_model, start = i, end = i, 
    #                       get.w = T, verbose = F, n.report = 1)
    #       beta_list[[length(beta_list) + 1]] <- as.matrix(co$p.beta.recover.samples)
    #       w_list[[length(w_list) + 1]] <- as.numeric(co$p.w.recover.samples)
    #       theta_list[[length(theta_list) + 1]] <- as.matrix(co$p.theta.recover.samples)
    #       used_idx <- c(used_idx, i)
    #     }, error = function(e) {
    #     })
    #   }
    #   if (length(beta_list) == 0) {
    #     message("No samples recovered successfully.")
    #     return(NULL)
    #   }
    #   beta_mat <- do.call(rbind, beta_list)
    #   w_mat <- do.call(cbind, w_list)
    #   theta_mat <- do.call(rbind, theta_list)
    #   beta_mcmc <- mcmc(beta_mat, start = 1, end = nrow(beta_mat), thin = 1)
    #   theta_mcmc <- mcmc(theta_mat, start = 1, end = nrow(theta_mat), thin = 1)
    #   return(list(
    #     p.beta.recover.samples = beta_mcmc,
    #     p.theta.recover.samples = theta_mcmc,
    #     p.w.recover.samples = w_mat,
    #     used_samples = nrow(beta_mat)
    #   ))
    # }
    
    run_spLM_chains <- function(
    formula, data, coords, knots,
    starting_list, tuning, priors,
    cov.model = "matern", n.samples = n.samples,
    burn.in = burn.in,
    param_pattern = "temp_",
    verbose = TRUE, use_knots,
    pred.coords, pred.covar
    ) {
      n.chains <- length(starting_list)
      if (is.null(burn.in)) burn.in <- floor(0.5 * n.samples)
      chain_results <- vector("list", n.chains)
      chain_recover <- vector("list", n.chains)
      predict_w_list <- vector("list", n.chains)
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
          amcmc = list("n.batch" = 5,"batch.length" = 2000, accept.rate = 0.43),
          verbose = F
        )
        if (use_knots) args_list$knots <- knots
        if (verbose) message("Running spLM")
        sp_model <- do.call(spLM, args_list)
        # chain_object <- safe_spRecover_beta(sp_model, burn.in = burn.in, n.samples = n.samples)
        if (verbose) message("Recovering samples")
        chain_object <- spRecover(sp_model, start = burn.in, get.w = T, verbose = F)
        if (is.null(chain_object)) {
          message("No samples recovered successfully for chain ", i)
          next  # Skip the current iteration of the loop
        }
        if (verbose) message("Predicting spatial random effect")
        predict_w <- spPredict(sp_model, pred.coords, pred.covar, start = 9000, verbose = F)
        predict_w_list[[i]] <- predict_w
        chain_recover[[i]] <- chain_object
        # param_cols <- grep(param_pattern, colnames(chain_recover[[i]]$p.beta.recover.samples))
        chain <- chain_object$p.beta.recover.samples
        # varnames(chain) <- "beta_temp"
        chain_results[[i]] <- mcmc(chain)
      }
      if ((length(chain_results[[1]]) + length(chain_results[[2]])) == 0) {
        return(NULL)
      }
      # used_sets <- lapply(chain_recover, `[[`, "used_samples")
      # min_len <- min(sapply(chain_recover, function(x) x$used_samples))
      # chain_results_aligned <- lapply(chain_recover, function(ch) {
      #   beta_mat <- as.matrix(ch$p.beta.recover.samples)[1:min_len, , drop = FALSE]
      #   mcmc(beta_mat, start = 1, end = min_len, thin = 1)
      # })
      beta_mcmc_list <- mcmc.list(chain_results)
      
      return(list(
        beta_mcmc_list = beta_mcmc_list,
        chains_recover = chain_recover,
        predict_w = predict_w_list
      ))
    }
    
    res <- run_spLM_chains(
      formula = sp_formula,
      data = tavg_phe_df,
      coords = coords,
      knots = knots,
      starting_list = starting_list,
      tuning = tuning,
      priors = priors,
      cov.model = "matern",
      n.samples = n.samples,
      burn.in = burn.in,
      param_pattern = "temp_",
      use_knots = use_knots,
      pred.coords = pred.coords, 
      pred.covar = pred.covar
    )
    if (is.null(res)) {
      next  # Skip the current iteration of the loop
    }
    beta2_mcmc <- lapply(res$beta_mcmc_list, function(chain) {
      temp_cols <- grep("^temp_", colnames(chain), value = TRUE)  
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
    
    #   ############################ Calculate marginal log-likelihood ############################  
    # log_posterior_spLM <- function(params, data) {
    #   y <- data$y
    #   X <- data$X
    #   priors <- data$priors
    #   dist_mat_nm <- data$dist_mat_nm
    #   dist_mat_mm <- data$dist_mat_mm
    #   
    #   beta_names <- colnames(X)
    #   beta <- as.numeric(params[beta_names])
    #   sigma.sq <- params["sigma.sq"]
    #   tau.sq   <- params["tau.sq"]
    #   phi      <- params["phi"]
    #   nu       <- params["nu"]
    #   
    #   if (any(is.na(c(beta, sigma.sq, tau.sq, phi, nu))) || sigma.sq <= 0 || tau.sq <= 0) {
    #     return(-Inf)
    #   }
    #   n <- length(y)
    #   m <- nrow(dist_mat_mm)
    #   matern_cov <- function(d, phi, nu) {
    #     d <- ifelse(d == 0, 1e-10, d)
    #     factor <- (2^(1-nu))/gamma(nu)
    #     cov <- factor * ( (sqrt(2*nu)*d/phi)^nu ) * besselK((sqrt(2*nu)*d/phi), nu)
    #     return(cov)
    #   }
    #   Rnm <- matern_cov(dist_mat_nm, phi, nu)
    #   Rmm <- matern_cov(dist_mat_mm, phi, nu) + 1e-6 * diag(m)
    #   Lmm <- chol(Rmm)
    #   Q <- (sigma.sq^(0.5)) * forwardsolve(t(Lmm), t(Rnm))  # m x n
    #   mu <- X %*% beta
    #   
    #   logdet_Sigma <- n * log(tau.sq) +
    #     determinant(diag(m) + (1/tau.sq) * Q %*% t(Q), logarithm=TRUE)$modulus
    #   
    #   resid <- y - mu
    #   alpha <- resid / tau.sq -
    #     (1 / tau.sq^2) * t(Q) %*% solve(diag(m) + (1/tau.sq) * Q %*% t(Q), Q %*% resid)
    #   
    #   quadform <- sum(resid * alpha)
    #   loglik <- -0.5 * (quadform + logdet_Sigma + n * log(2 * pi))
    #   
    #   logprior <- mvtnorm::dmvnorm(beta,
    #                                mean = priors$beta.Norm[[1]],
    #                                sigma = priors$beta.Norm[[2]],
    #                                log = TRUE)
    #   
    #   # 先验：sigma.sq 和 tau.sq 逆伽马
    #   logprior <- logprior + invgamma::dinvgamma(sigma.sq,
    #                                              shape = priors$sigma.sq.IG[1],
    #                                              scale = priors$sigma.sq.IG[2],
    #                                              log = TRUE)
    #   logprior <- logprior + invgamma::dinvgamma(tau.sq,
    #                                              shape = priors$tau.sq.IG[1],
    #                                              scale = priors$tau.sq.IG[2],
    #                                              log = TRUE)
    #   
    #   # 先验：phi uniform
    #   if (phi < priors$phi.Unif[1] || phi > priors$phi.Unif[2]) {
    #     return(-Inf)
    #   } else {
    #     logprior <- logprior + log(1 / (priors$phi.Unif[2] - priors$phi.Unif[1]))
    #   }
    #   if (!is.null(priors$nu.Unif)) {
    #     if (nu < priors$nu.Unif[1] || nu > priors$nu.Unif[2]) return(-Inf)
    #     logprior <- logprior + log(1 / (priors$nu.Unif[2] - priors$nu.Unif[1]))
    #   }
    #   return(loglik + logprior)
    # }
    # 
    # calc_marginal_likelihood <- function(chains_recover, log_posterior_fun, y, X, priors, dist_mat_nm, dist_mat_mm) {
    #   message("Calculating marginal likelihood...")
    #   sample_list <- lapply(chains_recover, function(chain) {
    #     cbind(chain$p.beta.recover.samples, chain$p.theta.recover.samples)
    #   })
    #   samples_all <- do.call(rbind, sample_list)
    #   samples <- as.matrix(samples_all)
    #   param_names <- colnames(samples)
    #   
    #   lb <- setNames(rep(-Inf, length(param_names)), param_names)
    #   ub <- setNames(rep( Inf, length(param_names)), param_names)
    #   if ("sigma.sq" %in% param_names) lb["sigma.sq"] <- 0
    #   if ("tau.sq"   %in% param_names) lb["tau.sq"]   <- 0
    #   if ("phi"      %in% param_names) {
    #     lb["phi"] <- priors$phi.Unif[1]
    #     ub["phi"] <- priors$phi.Unif[2]
    #   }
    #   if ("nu" %in% param_names) {
    #     lb["nu"] <- priors$nu.Unif[1]
    #     ub["nu"] <- priors$nu.Unif[2]
    #   }
    #   
    #   bridge_result <- bridgesampling::bridge_sampler(
    #     samples = samples,
    #     log_posterior = log_posterior_fun,
    #     data = list(y = y, X = X, priors = priors, dist_mat_nm = dist_mat_nm, dist_mat_mm = dist_mat_mm),
    #     lb = lb,
    #     ub = ub,
    #     param_types = rep("real", ncol(samples)),
    #     repetitions = 5,
    #     method = "normal",
    #     cores = 1,
    #     use_neff = TRUE,
    #     silent = FALSE,
    #     verbose = TRUE,
    #     varlist = ls(),
    #     envir = environment()
    #   )
    #   return((bridge_result$logml))
    # }
    # 
    # res_logml <- calc_marginal_likelihood(chains_recover = res$chains_recover, log_posterior_fun = log_posterior_spLM, 
    #                                       y = y, X = X, priors = priors, dist_mat_nm = dist_mat_nm, dist_mat_mm = dist_mat_mm)
    
    ############################ Calculate WAIC ############################
    
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
    
    ############################ Plot covariance function ############################
    
    samps <- do.call(
      rbind,
      lapply(res$chains_recover, function(chain) chain$p.theta.recover.samples)
    ) |> as.data.frame()
    matern_cor <- function(h, phi, nu) {
      a <-  h * phi
      out <- (2^(1 - nu) / gamma(nu)) * (a^nu) * besselK(a, nu)
      out[a == 0] <- 1
      out
    }
    dseq <- seq(0, 10000, length.out = 200)
    
    ## Latent field correlation ρ(h) (no sigma.sq, no tau.sq)
    rho_mat <- matrix(NA_real_, nrow = nrow(samps), ncol = length(dseq))
    for (k in seq_len(nrow(samps))) {
      rho_mat[k, ] <- matern_cor(dseq, samps[k, "phi"], samps[k, "nu"])
    }
    
    ## Credible band for correlation
    ci <- apply(rho_mat, 2, quantile, c(0.5, 0.025, 0.975), na.rm = TRUE)
    CorDF <- data.frame(
      Distance = dseq,
      Median   = ci[1, ],
      Lower    = ci[2, ],
      Upper    = ci[3, ]
    )
    plot_cor <- ggplot(CorDF, aes(x = Distance, y = Median)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey", alpha = 0.4) +
      geom_line(color = "blue", size = 1) +
      labs(x = "Distance", y = "Correlation") +
      theme_minimal(base_size = 14)
    
    
    ############################ Plot spatial random effect - w ############################      
    
    Y_pred <- do.call(
      rbind,
      lapply(res$predict_w, function(chain) chain$p.y.predictive.samples |> t())
    ) |> t()
    B_draws <- do.call(rbind, lapply(res$chains_recover, function(ch) {
      as.matrix(ch$p.beta.recover.samples[4001:5001,])
    }))  # S_total x p
    Xbeta <- pred.covar %*% t(B_draws)
    W_mat <- Y_pred - Xbeta     # n_pred x S
    w_mean    <- rowMeans(W_mat)
    
    surf_sf <- st_as_sf(data.frame(w_mean = w_mean), geometry = grid_poly_in)
    plot_spatial_field <- ggplot(surf_sf) +
      geom_sf(aes(fill = w_mean), color = NA) +
      scale_fill_viridis_c() +
      coord_sf(expand = FALSE) +
      labs(fill = "E[w|data]") +
      theme_minimal(base_size = 12)
    
    
    ############################ Save results ############################   
    
    modelMCMCrun <- list(
      samples = res$beta_mcmc_list,
      WAIC = list(WAIC = res_waic),
      converge = converge$psrf[1],
      plot_cor = plot_cor,
      plot_spatial_field = plot_spatial_field
      # mean_logml = mean(res_logml)
    )
    write_rds(modelMCMCrun, file = model_path, compress = "gz")
    
    list(taxon = taxon, ok = TRUE)
  }
  
  res
}

