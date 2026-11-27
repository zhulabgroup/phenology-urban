# NewYork SOS / per-species INLA diagnostics
# Stage 1: non-spatial vs SPDE (residual Moran's I)
# Stage 2: linear vs natural-spline temperature effect, CONDITIONING on chosen spatial structure
# Reuses the saved model mesh/SPDE; diagnostic coordinates are projected to EPSG:5070.
# This is a diagnostic script: it does not overwrite the manuscript model files.

suppressPackageStartupMessages({
  library(INLA)
  library(sf)
  library(spdep)
  library(ggplot2)
  library(splines)
  library(tidyverse)
  library(here)
})

source(here::here("scripts/04_stats_model/fit_spatial/prepare_model_input.R"))

# =========== SETTINGS  =========================================
CITY           <- "NewYork"
SOURCE_NAME    <- "Grid"     # Set to the source used for the INLA manuscript fit
PHE_TYPE       <- "SOS"         # Confirm the actual SOS label
SPECIES        <- "Acer saccharinum"
FITTING_TYPE   <- "doy"            # "doy", "start", or "end"

# Use "batch_moran" for every saved city-species model.
RUN_STAGE      <- "batch_moran"    # "spatial", "linearity", "both", "batch_moran"
LINEARITY_WITH_SPDE <- TRUE

# Leave blank to autodetect only when EXACTLY ONE column has the prefix.
TEMP_COLUMN    <- ""               # e.g. "temp_winter"
PRCP_COLUMN    <- ""               # e.g. "prcp_winter"
SRAD_COLUMN    <- ""               # e.g. "srad_winter"
SPLINE_DF      <- 4L               # Natural spline with 4 fixed-effect basis columns
MORAN_K        <- 8L
MORAN_NSIM     <- 999L
MORAN_MAX_N    <- 15000L            # Deterministic subsample if > 15000 unique points/year
CPO_COMPUTE    <- FALSE            # Original model had TRUE; disable for faster diagnostics
MODEL_ROOT     <- "data/intermediate/spatial/inla/spatial_spde"
RESULT_DIR     <- "data/output/sensitive_analysis"
MODEL_RDS      <- file.path("data/intermediate/spatial/inla/spatial_spde", CITY,
                            "canonicalname",
                            paste0(gsub(" ", "_", SPECIES), "_", PHE_TYPE, ".rds"))

# =========== FUNCTIONS  =========================================
say <- function(...) cat(paste0(..., "\n"))

get_column <- function(d, given, prefix, season=NULL) {
  if (nzchar(given)) {
    if (!given %in% names(d)) stop("Missing requested ", prefix, " column: ", given)
    return(given)
  }
  if (!is.null(season)) {
    seasonal_column <- paste0(prefix, "_", season)
    if (seasonal_column %in% names(d)) return(seasonal_column)
  }
  hits <- grep(paste0("^", prefix), names(d), value=TRUE)
  if (length(hits) != 1L)
    stop("Expected exactly one '", prefix, "*' column; found: ",
         paste(hits, collapse=", "), ". Set the corresponding *_COLUMN explicitly.")
  hits
}

# Do the SAME response IQR exclusion, factor(year), and lon/lat projection
# as the uploaded alg_sp_inla.R, lines 50-75. Keep WGS84 coordinates separately.
prepare_input <- function(original, species=SPECIES, phe_type=PHE_TYPE,
                          validate_for_spline=TRUE) {
  if (inherits(original, "sf")) original <- sf::st_drop_geometry(original)
  d <- as.data.frame(original)
  if (!"canonicalname" %in% names(d)) stop("Missing canonicalname")
  if (!all(c("lon","lat","year") %in% names(d)))
    stop("Missing one of lon, lat, year")
  d <- d[!is.na(d$canonicalname) & as.character(d$canonicalname)==species, , drop=FALSE]
  if (!nrow(d)) stop("No rows found for ", species)

  response_col <- if (FITTING_TYPE == "start") "start" else if (FITTING_TYPE == "end") "end" else
    if (grepl("^[SE]OS", phe_type)) "doy" else "doy_diff"
  if (!response_col %in% names(d)) stop("Missing response column: ", response_col)
  season <- if (grepl("^SOS", phe_type) || phe_type == "GSL_winter") "winter" else
    if (grepl("^EOS", phe_type) || phe_type == "GSL_summer") "summer" else
      if (phe_type == "GSL_range_season") "growing" else
        if (phe_type == "GreenUp") "spring" else "fall"
  tmp <- get_column(d, TEMP_COLUMN, "temp", season)
  pc <- get_column(d, PRCP_COLUMN, "prcp", season)
  sr <- get_column(d, SRAD_COLUMN, "srad", season)
  say("Variables: response=", response_col, ", temp=", tmp, ", prcp=", pc, ", srad=", sr)

  yy <- suppressWarnings(as.numeric(as.character(d[[response_col]])))
  qs <- stats::quantile(yy, c(.25,.75), na.rm=TRUE)
  iqr <- diff(qs)
  lo <- as.numeric(qs[1] - 1.5*iqr)
  hi <- as.numeric(qs[2] + 1.5*iqr)
  keep_y <- is.finite(yy) & yy >= lo & yy <= hi
  d <- d[keep_y, ,drop=FALSE]
  yy <- yy[keep_y]
  say("After the original per-species 1.5*IQR filter: n=", nrow(d))

  vals <- data.frame(
    y=yy,
    temp=as.numeric(d[[tmp]]),
    prcp=as.numeric(d[[pc]]),
    srad=as.numeric(d[[sr]]),
    lon_wgs84=as.numeric(d$lon),
    lat_wgs84=as.numeric(d$lat),
    year=as.factor(d$year)
  )
  good <- stats::complete.cases(vals) &
    is.finite(vals$y) & is.finite(vals$temp) & is.finite(vals$prcp) &
    is.finite(vals$srad) & is.finite(vals$lon_wgs84) & is.finite(vals$lat_wgs84) &
    vals$lon_wgs84 >= -180 & vals$lon_wgs84 <= 180 &
    vals$lat_wgs84 >= -90 & vals$lat_wgs84 <= 90
  vals <- droplevels(vals[good, ,drop=FALSE])
  if (validate_for_spline && nrow(vals)<50L) stop("Too few complete records after preprocessing")
  if (validate_for_spline && length(unique(vals$temp))<SPLINE_DF+3L)
    stop("Too few distinct temperatures for spline")
  if (validate_for_spline &&
      (length(unique(vals$lon_wgs84))<3L || length(unique(vals$lat_wgs84))<3L))
    stop("Insufficient spatial extent for SPDE mesh")

  # Exact spatial logic from alg_sp_inla.R:
  # st_as_sf(coords=c('lon','lat'), crs=4326) |> st_transform(crs=5070)
  xy_sf <- sf::st_as_sf(vals, coords=c("lon_wgs84","lat_wgs84"), crs=4326, remove=FALSE)
  xy_sf <- sf::st_transform(xy_sf, crs=5070)
  xy <- sf::st_coordinates(xy_sf)
  vals$x <- as.numeric(xy[,1])
  vals$z <- as.numeric(xy[,2])  # Northing; call z to avoid conflict with response y
  vals$Intercept <- 1
  vals$obs_id <- seq_len(nrow(vals))
  say("Complete matched data: n=", nrow(vals), ", years=", length(unique(vals$year)),
      ", unique locations=", nrow(unique(vals[, c("x","z")])))
  vals
}

# Load the saved spatial fit and verify that the prepared response rows match
# the original model's prediction rows in both count and order.
load_original_fit <- function(dat) {
  if (!file.exists(MODEL_RDS)) {
    stop("Original INLA model file not found: ", MODEL_RDS)
  }
  original_fit <- readRDS(MODEL_RDS)
  if (!identical(as.character(original_fit$species), SPECIES)) {
    stop("Species in original model file does not match SPECIES: ", MODEL_RDS)
  }
  if (is.null(original_fit$mesh) || is.null(original_fit$spde)) {
    stop("Original model file must contain both mesh and spde: ", MODEL_RDS)
  }
  if (is.null(original_fit$compare) ||
      !all(c("value", "true_value") %in% rownames(original_fit$compare))) {
    stop("Original model file must contain compare rows 'value' and 'true_value'")
  }
  original_y <- suppressWarnings(as.numeric(original_fit$compare["true_value", ]))
  if (length(original_y) != nrow(dat) ||
      !isTRUE(all.equal(dat$y, original_y, tolerance = 1e-10))) {
    stop("Prepared observations do not match the saved model's response rows/order. ",
         "Do not reuse saved predictions until the input filtering/order is aligned.")
  }
  say("Loaded saved spatial fit from ", MODEL_RDS,
      " (matched observations: ", nrow(dat), ")")
  original_fit
}

saved_spatial_fit <- function(dat, original_fit) {
  pred <- suppressWarnings(as.numeric(original_fit$compare["value", ]))
  if (length(pred) != nrow(dat) || any(!is.finite(pred))) {
    stop("Saved spatial predictions do not match the prepared observations")
  }
  res <- original_fit$res
  list(res=res, fitted=pred, residual=dat$y-pred,
       waic=as.numeric(res$waic$waic), dic=as.numeric(res$dic$dic),
       beta_temp=as.numeric(res$summary.fixed["temp", "mean"]),
       spatial=TRUE, nonlinear=FALSE)
}

# Reuse the exact mesh/SPDE from the saved model; create A for current rows.
build_spde <- function(dat, original_fit) {

  coord <- as.matrix(dat[,c("x","z")])
  mesh <- original_fit$mesh
  spde <- original_fit$spde
  A <- INLA::inla.spde.make.A(mesh=mesh, loc=coord)
  if (nrow(A) != nrow(dat)) {
    stop("Projection matrix rows do not match diagnostic observations")
  }
  if (any(rowSums(A) == 0)) {
    stop("Some diagnostic coordinates fall outside the saved model mesh")
  }
  s_index <- INLA::inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde)
  say("Reusing original mesh/SPDE from ", MODEL_RDS,
      " (mesh nodes: ", mesh$n, ")")
  list(mesh=mesh,spde=spde,A=A,s_index=s_index)
}

# Fit all models using the same Gaussian likelihood, covariates, year RE,
# mesh, SPDE hyperpriors and INLA strategy. No new independent smooth in stage 2.
fit_inla <- function(dat, geo=NULL, spatial=FALSE, nonlinear=FALSE, basis=NULL) {
  cols <- if (nonlinear) paste0("ns",seq_len(ncol(basis))) else "temp"
  fx <- data.frame(Intercept=1, temp=dat$temp, prcp=dat$prcp,
                   srad=dat$srad,year=dat$year)
  if (nonlinear) {
    if (is.null(basis) || nrow(basis)!=nrow(dat)) stop("Bad spline basis")
    fx$temp <- NULL
    for (j in seq_len(ncol(basis))) fx[[cols[j]]] <- as.numeric(basis[,j])
  }
  # INLA 25.10.19 fails on the original shape=.01 iid precision prior here.
  formula_terms <- c("Intercept",cols,"prcp","srad",
    "f(year, model='iid', hyper=list(prec=list(prior='loggamma', param=c(1,5e-5), initial=4, fixed=FALSE)))")
  if (spatial) formula_terms <- c(formula_terms,"f(spatial.field, model=spde)")
  fml <- stats::as.formula(paste("y ~ 0 +", paste(formula_terms,collapse=" + ")),
                           env=environment())

  all_names <- c("Intercept",cols,"prcp","srad")
  prior_mean <- stats::setNames(rep(0,length(all_names)),all_names)
  prior_mean["Intercept"] <- 180
  prior_prec <- stats::setNames(rep(1/(5^2),length(all_names)),all_names)
  prior_prec["Intercept"] <- 1/(40^2)
  prior_prec[c("prcp","srad")] <- 1/(1^2)
  # Spline basis columns are dimensionless; use broader priors than temp slope.
  if (nonlinear) prior_prec[cols] <- 1/(20^2)
  priors <- list(mean=as.list(prior_mean),prec=as.list(prior_prec))

  if (spatial) {
    if (is.null(geo)) stop("geo required for spatial fit")
    spde <- geo$spde   # available in formula environment
    stk <- INLA::inla.stack(tag="est", data=list(y=dat$y),
      A=list(1, geo$A),effects=list(fx, geo$s_index))
    res <- INLA::inla(fml,family="gaussian",data=INLA::inla.stack.data(stk),
      control.inla=list(strategy="laplace",interpolator="gaussian",fast=FALSE,restart=5),
      control.fixed=priors,
      control.compute=list(waic=TRUE,dic=TRUE,cpo=CPO_COMPUTE),
      control.predictor=list(compute=TRUE,A=INLA::inla.stack.A(stk)))
    ind <- INLA::inla.stack.index(stack=stk,tag="est")$data
  } else {
    res <- INLA::inla(fml,family="gaussian",data=cbind(data.frame(y=dat$y),fx),
      control.inla=list(strategy="laplace",interpolator="gaussian",fast=FALSE,restart=5),
      control.fixed=priors,
      control.compute=list(waic=TRUE,dic=TRUE,cpo=CPO_COMPUTE),
      control.predictor=list(compute=TRUE))
    ind <- seq_len(nrow(dat))
  }
  pred <- as.numeric(res$summary.fitted.values[ind,"mean"])
  if (length(pred)!=nrow(dat) || any(!is.finite(pred)))
    stop("Fitted values do not align with input rows; inspect INLA stack/index")
  residual <- dat$y - pred
  list(res=res, fitted=pred, residual=residual,
       waic=as.numeric(res$waic$waic),dic=as.numeric(res$dic$dic),
       beta_temp=if (!nonlinear) as.numeric(res$summary.fixed["temp","mean"]) else NA_real_,
       spatial=spatial,nonlinear=nonlinear)
}

# Aggregate coincident coordinates WITHIN year so each site/year is not counted
# multiple times in nearest-neighbour graph (no pseudo-replication by location).
moran_by_year <- function(dat, residual, model_name) {
  dd <- data.frame(year=as.character(dat$year), x=dat$x, z=dat$z,
                   residual=as.numeric(residual))
  years <- sort(unique(dd$year))
  set.seed(3107)
  ans <- lapply(years, function(yr) {
    s <- dd[dd$year==yr,c("x","z","residual"),drop=FALSE]
    s <- stats::aggregate(residual ~ x+z, data=s, FUN=mean)
    n_all <- nrow(s)
    if (n_all > MORAN_MAX_N) s <- s[sample.int(n_all,MORAN_MAX_N),,drop=FALSE]
    n_use <- nrow(s)
    k <- min(MORAN_K,n_use-1L)
    if (k < 2L || stats::sd(s$residual)==0) {
      return(data.frame(model=model_name,year=yr,n_unique_all=n_all,n_used=n_use,
                        k=k,moran_I=NA_real_,p_value=NA_real_))
    }
    xy <- as.matrix(s[,c("x","z")])
    tryCatch({
      knn <- spdep::knearneigh(xy,k=k,longlat=FALSE)
      nb <- spdep::knn2nb(knn,sym=TRUE)
      lw <- spdep::nb2listw(nb,style="W",zero.policy=TRUE)
      mc <- spdep::moran.mc(s$residual,lw,nsim=MORAN_NSIM,
                            alternative="greater",zero.policy=TRUE)
      data.frame(model=model_name,year=yr,n_unique_all=n_all,n_used=n_use,
                 k=k,moran_I=unname(as.numeric(mc$statistic)),p_value=mc$p.value)
    },error=function(e) {
      warning("Moran's I failed for ", yr, ": ", conditionMessage(e))
      data.frame(model=model_name,year=yr,n_unique_all=n_all,n_used=n_use,
                 k=k,moran_I=NA_real_,p_value=NA_real_)
    })
  })
  do.call(rbind,ans)
}

moran_pair_by_year <- function(dat, residual_spatial, residual_nonspatial) {
  dd <- data.frame(year=as.character(dat$year), x=dat$x, z=dat$z,
                   residual_spatial=as.numeric(residual_spatial),
                   residual_nonspatial=as.numeric(residual_nonspatial))
  years <- sort(unique(dd$year))
  set.seed(3107)
  ans <- lapply(years, function(yr) {
    site_residuals <- stats::aggregate(
      cbind(residual_spatial, residual_nonspatial) ~ x + z,
      data=dd[dd$year == yr, ], FUN=mean
    )
    n_all <- nrow(site_residuals)
    if (n_all > MORAN_MAX_N) {
      site_residuals <- site_residuals[sample.int(n_all, MORAN_MAX_N), , drop=FALSE]
    }
    n_use <- nrow(site_residuals)
    k <- min(MORAN_K, n_use - 1L)
    result <- data.frame(
      year=yr, n_unique_sites=n_all, n_used=n_use, k=k,
      moran_spatial=NA_real_, moran_nonspatial=NA_real_,
      p_analytic_spatial=NA_real_, p_analytic_nonspatial=NA_real_,
      p_two_sided_spatial=NA_real_, p_two_sided_nonspatial=NA_real_,
      delta_nonspatial_minus_spatial=NA_real_
    )
    if (k < 2L) return(result)

    xy <- as.matrix(site_residuals[, c("x", "z")])
    tryCatch({
      nb <- spdep::knn2nb(spdep::knearneigh(xy, k=k), sym=TRUE)
      lw <- spdep::nb2listw(nb, style="W", zero.policy=TRUE)
      moran_stat <- function(residual) {
        centered <- residual - mean(residual)
        denominator <- sum(centered^2)
        if (denominator == 0) return(NA_real_)
        s0 <- sum(vapply(lw$weights, sum, numeric(1)))
        if (s0 == 0) return(NA_real_)
        lagged <- spdep::lag.listw(lw, centered, zero.policy=TRUE)
        n_use / s0 * sum(centered * lagged) / denominator
      }
      result$moran_spatial <- moran_stat(site_residuals$residual_spatial)
      result$moran_nonspatial <- moran_stat(site_residuals$residual_nonspatial)
      result$p_analytic_spatial <- spdep::moran.test(
        site_residuals$residual_spatial, lw,
        alternative="greater", zero.policy=TRUE
      )$p.value
      result$p_analytic_nonspatial <- spdep::moran.test(
        site_residuals$residual_nonspatial, lw,
        alternative="greater", zero.policy=TRUE
      )$p.value
      result$p_two_sided_spatial <- spdep::moran.test(
        site_residuals$residual_spatial, lw,
        alternative="two.sided", zero.policy=TRUE
      )$p.value
      result$p_two_sided_nonspatial <- spdep::moran.test(
        site_residuals$residual_nonspatial, lw,
        alternative="two.sided", zero.policy=TRUE
      )$p.value
      result$delta_nonspatial_minus_spatial <-
        result$moran_nonspatial - result$moran_spatial
      result
    }, error=function(e) {
      warning("Paired Moran's I failed for year ", yr, ": ", conditionMessage(e))
      result
    })
  })
  do.call(rbind, ans)
}

save_stage1_tables <- function(moran, model_comparison, result_dir=RESULT_DIR) {
  dir.create(result_dir, recursive=TRUE, showWarnings=FALSE)
  utils::write.csv(moran, file.path(result_dir, "stage1_moran_by_year.csv"),
                   row.names=FALSE)
  utils::write.csv(model_comparison,
                   file.path(result_dir, "stage1_model_fit_comparison.csv"),
                   row.names=FALSE)
  invisible(result_dir)
}

run_batch_moran <- function(model_root=MODEL_ROOT, source=SOURCE_NAME,
                            model_files=NULL) {
  if (!exists("prepare_data", mode="function", inherits=TRUE)) {
    stop("Source prepare_model_input.R first so prepare_data() is available")
  }
  if (is.null(model_files)) {
    model_files <- list.files(
      model_root, pattern="_(SOS|EOS)(_[0-9]+)?[.]rds$",
      recursive=TRUE, full.names=TRUE
    )
  }
  if (!length(model_files)) stop("No city-species model RDS files found in ", model_root)

  input_cache <- new.env(parent=emptyenv())
  result_list <- vector("list", length(model_files))
  error_list <- list()
  error_count <- 0L

  for (i in seq_along(model_files)) {
    model_path <- model_files[[i]]
    city <- basename(dirname(dirname(model_path)))
    phe_type <- sub(".*_((SOS|EOS)(_[0-9]+)?)[.]rds$", "\\1", basename(model_path))
    species <- NA_character_

    result_list[[i]] <- tryCatch({
      original_fit <- readRDS(model_path)
      species <- as.character(original_fit$species)
      if (length(species) != 1L || is.na(species) || !nzchar(species)) {
        stop("Model RDS does not contain one valid species name")
      }
      cache_key <- paste(city, phe_type, sep="__")
      if (!exists(cache_key, envir=input_cache, inherits=FALSE)) {
        city_input <- prepare_data(city, source=source,
                                   taxa_col="canonicalname", phe_type=phe_type)
        assign(cache_key, city_input, envir=input_cache)
      }
      dat <- prepare_input(get(cache_key, envir=input_cache),
                           species=species, phe_type=phe_type,
                           validate_for_spline=FALSE)
      if (is.null(original_fit$compare) ||
          !all(c("value", "true_value") %in% rownames(original_fit$compare))) {
        stop("Model RDS lacks compare rows 'value' and 'true_value'")
      }
      original_y <- as.numeric(original_fit$compare["true_value", ])
      spatial_pred <- as.numeric(original_fit$compare["value", ])
      if (length(original_y) != nrow(dat) ||
          !isTRUE(all.equal(dat$y, original_y, tolerance=1e-10))) {
        stop("Prepared response rows/order do not match saved model predictions")
      }
      if (all(c("Lon", "Lat") %in% rownames(original_fit$compare))) {
        original_x <- as.numeric(original_fit$compare["Lon", ])
        original_z <- as.numeric(original_fit$compare["Lat", ])
        if (length(original_x) != nrow(dat) || length(original_z) != nrow(dat) ||
            !isTRUE(all.equal(dat$x, original_x, tolerance=1e-6)) ||
            !isTRUE(all.equal(dat$z, original_z, tolerance=1e-6))) {
          stop("Prepared coordinates/order do not match the saved model")
        }
      }
      if (length(spatial_pred) != nrow(dat) || any(!is.finite(spatial_pred))) {
        stop("Saved spatial predictions do not match prepared observations")
      }

      stage1(dat, original_fit, city, species, phe_type, model_path)
    }, error=function(e) {
      error_count <<- error_count + 1L
      error_list[[error_count]] <<- data.frame(
        city=city, species=species, phe_type=phe_type,
        model_file=model_path, error=conditionMessage(e),
        stringsAsFactors=FALSE
      )
      NULL
    })

    if (i %% 100L == 0L) {
      say("Processed ", i, " / ", length(model_files), " model files; failures: ",
          error_count)
      invisible(gc())
    }
  }

  valid_results <- Filter(Negate(is.null), result_list)
  moran_results <- if (length(valid_results)) {
    do.call(rbind, lapply(valid_results, `[[`, "moran"))
  } else data.frame()
  model_comparison <- if (length(valid_results)) {
    do.call(rbind, lapply(valid_results, `[[`, "model_comparison"))
  } else data.frame()
  errors <- if (length(error_list)) do.call(rbind, error_list) else data.frame()
  save_stage1_tables(moran_results, model_comparison)
  say("Stage 1 tables saved: ", RESULT_DIR)
  say("Successful combinations: ", nrow(model_comparison),
      "; failed combinations: ", nrow(errors))
  list(moran=moran_results, model_comparison=model_comparison, errors=errors)
}

stage1 <- function(dat, original_fit, city=CITY, species=SPECIES,
                   phe_type=PHE_TYPE, model_path=MODEL_RDS) {
  nonspatial_fit <- fit_inla(dat, spatial=FALSE)
  spatial_fit <- saved_spatial_fit(dat, original_fit)
  spatial_residual <- spatial_fit$residual
  nonspatial_residual <- nonspatial_fit$residual

  moran <- moran_pair_by_year(dat, spatial_residual, nonspatial_residual)
  moran$city <- city
  moran$species <- species
  moran$phe_type <- phe_type
  moran$n_observations <- nrow(dat)
  moran$nonspatial_year_re_prior <- "loggamma(1, 5e-5)"
  moran$model_file <- model_path
  moran <- moran[, c("city", "species", "phe_type", "year", "n_observations",
                     "n_unique_sites", "n_used", "k", "moran_spatial",
                     "p_analytic_spatial", "moran_nonspatial",
                     "p_analytic_nonspatial", "p_two_sided_spatial",
                     "p_two_sided_nonspatial", "delta_nonspatial_minus_spatial",
                     "nonspatial_year_re_prior", "model_file")]

  spatial_rmse <- sqrt(mean(spatial_residual^2))
  nonspatial_rmse <- sqrt(mean(nonspatial_residual^2))
  model_comparison <- data.frame(
    city=city, species=species, phe_type=phe_type,
    n_observations=nrow(dat),
    WAIC_spatial=spatial_fit$waic,
    WAIC_nonspatial=nonspatial_fit$waic,
    delta_WAIC_nonspatial_minus_spatial=
      nonspatial_fit$waic - spatial_fit$waic,
    DIC_spatial=spatial_fit$dic,
    DIC_nonspatial=nonspatial_fit$dic,
    delta_DIC_nonspatial_minus_spatial=nonspatial_fit$dic - spatial_fit$dic,
    RMSE_spatial=spatial_rmse,
    RMSE_nonspatial=nonspatial_rmse,
    delta_RMSE_nonspatial_minus_spatial=nonspatial_rmse - spatial_rmse,
    nonspatial_year_re_prior="loggamma(1, 5e-5)",
    model_file=model_path, stringsAsFactors=FALSE
  )

  list(moran=moran, model_comparison=model_comparison)
}

stage2 <- function(dat, original_fit) {
  say("--- STAGE 2: temperature linearity under CHOSEN spatial structure ---")
  linear <- if (LINEARITY_WITH_SPDE) {
    saved_spatial_fit(dat, original_fit)
  } else {
    fit_inla(dat, spatial=FALSE)
  }
  spatial_geo <- if (LINEARITY_WITH_SPDE) build_spde(dat, original_fit) else NULL
  basis <- splines::ns(dat$temp,df=SPLINE_DF)
  spline_fit <- fit_inla(dat,geo=spatial_geo,spatial=LINEARITY_WITH_SPDE,
                          nonlinear=TRUE,basis=basis)
  beta <- linear$beta_temp
  coef_name <- paste0("ns",seq_len(ncol(basis)))
  coef <- as.numeric(spline_fit$res$summary.fixed[coef_name,"mean"])
  if (any(!is.finite(coef))) stop("Nonfinite spline coefficients")
  T0 <- as.numeric(stats::median(dat$temp))
  temps <- sort(unique(c(seq(as.numeric(stats::quantile(dat$temp,.05)),
                             as.numeric(stats::quantile(dat$temp,.95)),length.out=201L),T0)))
  make_basis <- function(t) as.matrix(stats::predict(basis,newx=as.numeric(t)))
  # All curves are TEMPERATURE-ONLY partial effects, centered at T0;
  # spatial and year terms are held constant and cancel in contrasts.
  E0 <- as.numeric(make_basis(T0) %*% coef)
  E <- as.numeric(make_basis(temps) %*% coef) - E0
  h <- max(1e-4, diff(range(dat$temp))*1e-4)
  deriv <- as.numeric((make_basis(temps+h)-make_basis(temps-h))%*%coef)/(2*h)
  deriv_T0 <- as.numeric((make_basis(T0+h)-make_basis(T0-h))%*%coef)/(2*h)
  dat_effect <- as.numeric(make_basis(dat$temp)%*%coef)-E0
  linear_effect <- beta*(dat$temp-T0)
  metrics <- data.frame(model=c("Linear","Natural_spline"),
    spatial=LINEARITY_WITH_SPDE,n=nrow(dat),WAIC=c(linear$waic,spline_fit$waic),
    DIC=c(linear$dic,spline_fit$dic),
    in_sample_RMSE=c(sqrt(mean(linear$residual^2)),sqrt(mean(spline_fit$residual^2))))
  # write_csv(metrics,"02_linearity_model_comparison.csv")
  contrasts <- data.frame(
    reference_temperature=T0, OLS_slope=beta,
    spline_derivative_at_median=deriv_T0,
    absolute_slope_difference_at_median=abs(deriv_T0-beta),
    RMSE_partial_effect_difference_days=sqrt(mean((dat_effect-linear_effect)^2)),
    RMSE_derivative_difference_days_per_C=sqrt(mean((
      as.numeric((make_basis(dat$temp+h)-make_basis(dat$temp-h))%*%coef)/(2*h)-beta)^2)),
    delta_WAIC_linear_minus_spline=linear$waic-spline_fit$waic)
  # write_csv(contrasts,"02_linearity_key_metrics.csv")
  curve <- data.frame(temp_C=temps,linear_effect_days=beta*(temps-T0),
                      spline_effect_days=E, spline_Taylor_at_T0_days=deriv_T0*(temps-T0),
                      spline_derivative_days_per_C=deriv,
                      linear_slope_days_per_C=beta)
  # write_csv(curve,"02_temperature_effect_and_derivative.csv")

  curves_long <- rbind(
    data.frame(temp_C=temps,response=curve$linear_effect_days,model="Linear"),
    data.frame(temp_C=temps,response=curve$spline_effect_days,model="Natural spline"),
    data.frame(temp_C=temps,response=curve$spline_Taylor_at_T0_days,model="Spline tangent at median T"))
  p1 <- ggplot(curves_long,aes(temp_C,response,color=model))+
    geom_hline(yintercept=0,color="grey75")+geom_vline(xintercept=T0,linetype="dotted")+
    geom_line(linewidth=1)+theme_bw()+
    labs(title=paste(CITY,SPECIES,"temperature partial response"),
         subtitle=paste("Same",if(LINEARITY_WITH_SPDE)"INLA-SPDE" else "non-spatial INLA",
                        "structure; curves centered at median temperature"),
         x="Observed winter temperature (°C)",y="Effect on SOS relative to median T (days)")
  # save_plot(p1,"02_partial_temperature_response.png")
  p2 <- ggplot(curve,aes(x=temp_C,y=spline_derivative_days_per_C))+
    geom_hline(yintercept=beta,linetype="dashed",color="#3569b0",linewidth=.9)+
    geom_line(linewidth=1,color="black")+
    geom_vline(xintercept=T0,linetype="dotted")+
    theme_bw()+labs(title="Conditional temperature derivative vs. linear slope",
                    subtitle="Posterior-mean functions; bands NOT estimated in this script",
                    x="Winter temperature (°C)",y="Sensitivity (days / °C)")
  # save_plot(p2,"02_spline_derivative_vs_linear_slope.png")

  dd <- rbind(data.frame(temp=dat$temp, fitted=linear$fitted,residual=linear$residual,model="Linear"),
              data.frame(temp=dat$temp, fitted=spline_fit$fitted,residual=spline_fit$residual,model="Natural spline"))
  p3 <- ggplot(dd,aes(temp,residual))+
    geom_point(alpha=.13,size=.35)+
    geom_smooth(method="loess",se=FALSE,color="#3569b0",linewidth=.9)+
    geom_hline(yintercept=0,linetype="dashed")+
    facet_wrap(~model)+theme_bw()+
    labs(title="Residuals along the temperature gradient",x="Winter temperature (°C)",y="Residual (days)")
  # save_plot(p3,"02_residuals_vs_temperature.png")
  p4 <- ggplot(dd,aes(fitted,residual))+
    geom_point(alpha=.12,size=.35)+
    geom_smooth(method="loess",se=FALSE,color="#3569b0",linewidth=.9)+
    geom_hline(yintercept=0,linetype="dashed")+
    facet_wrap(~model)+theme_bw()+
    labs(title="Residuals versus fitted SOS",x="Fitted SOS (DOY)",y="Residual (days)")
  # save_plot(p4,"02_residuals_vs_fitted.png")
  qq <- do.call(rbind,lapply(split(dd,dd$model),function(d) {
    d <- d[order(d$residual),,drop=FALSE]
    data.frame(theoretical=stats::qnorm(stats::ppoints(nrow(d))),
               residual=d$residual,model=d$model[1])
  }))
  p5 <- ggplot(qq,aes(theoretical,residual))+geom_point(size=.5,alpha=.3)+
    geom_qq_line(aes(sample=residual),color="grey35")+
    facet_wrap(~model)+theme_bw()+
    labs(title="Residual Q–Q",x="Theoretical quantiles",y="Residual quantiles")
  # save_plot(p5,"02_residual_QQ.png")
  moran <- rbind(moran_by_year(dat,linear$residual,"Linear"),
                 moran_by_year(dat,spline_fit$residual,"Natural_spline"))
  # write_csv(moran,"02_moran_linear_vs_spline.csv")
  # write_csv(data.frame(obs_id=dat$obs_id,year=as.character(dat$year),
  #                      x=dat$x,z=dat$z,temp=dat$temp,y=dat$y,
  #                      fitted_linear=linear$fitted,residual_linear=linear$residual,
  #                      fitted_spline=spline_fit$fitted,residual_spline=spline_fit$residual),
  #           "02_linearity_residuals_matched.csv")
  # saveRDS(list(linear=linear,spline=spline_fit,basis=basis,
  #              spatial=LINEARITY_WITH_SPDE),
  #         file.path(OUT_DIR,"stage2_linearity_models.rds"),compress="xz")
  say("Temperature linearity: WAIC(linear) = ",round(linear$waic,2),
      "; WAIC(spline) = ",round(spline_fit$waic,2))
  say("beta_linear = ",signif(beta,4),
      "; derivative_spline(T0) = ",signif(deriv_T0,4)," days/°C")
  say("Linearity diagnostics complete; output files were not written.")
  invisible(TRUE)
}

# =========== RUN  =========================================

# Batch mode processes all saved city-species model files and keeps results in memory.
if (RUN_STAGE == "batch_moran") {
  batch_moran <- run_batch_moran()
  moran_all_city_species <- batch_moran$moran
  model_fit_comparison_all_city_species <- batch_moran$model_comparison
  moran_batch_errors <- batch_moran$errors
} else if (RUN_STAGE %in% c("spatial","linearity","both")) {
  # Prepare the shared observations once; both stages use the saved spatial fit.
  if (!exists("model_input_all",inherits=TRUE)) {
    if (!exists("prepare_data",mode="function",inherits=TRUE))
      stop("Provide model_input_all or load prepare_data() before sourcing script")
    model_input_all <- prepare_data(CITY,source=SOURCE_NAME,
                                    taxa_col="canonicalname",phe_type=PHE_TYPE)
  }
  dat <- prepare_input(model_input_all)
  original_fit <- load_original_fit(dat)
  if (RUN_STAGE %in% c("spatial","both")) {
    stage1_result <- stage1(dat, original_fit)
    moran_all_city_species <- stage1_result$moran
    model_fit_comparison_all_city_species <- stage1_result$model_comparison
    save_stage1_tables(moran_all_city_species,
                       model_fit_comparison_all_city_species)
  }
  if (RUN_STAGE %in% c("linearity","both")) stage2(dat, original_fit)
}
if (!RUN_STAGE %in% c("spatial","linearity","both","batch_moran")) {
  stop("Unknown RUN_STAGE value")
}
