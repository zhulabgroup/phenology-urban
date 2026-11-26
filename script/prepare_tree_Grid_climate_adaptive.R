# Adaptive SOS/EOS preseason climate at tree locations.
# Source this file from the phenology-urban project. No processing on source().
# Dependencies: here, readr, data.table, dplyr, sf, terra
suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(data.table)
  library(dplyr)
  library(sf)
  library(terra)
})

.adapt_col <- function(x, explicit, candidates, label) {
  n <- names(x)
  if (!is.null(explicit)) {
    if (length(explicit) != 1L || !explicit %in% n)
      stop(label, " column not found: ", explicit, "; available: ", paste(n, collapse = ", "))
    return(explicit)
  }
  norm <- function(z) gsub("[^a-z0-9]", "", tolower(z))
  hits <- unique(unlist(lapply(candidates, function(z) n[norm(n) == norm(z)])))
  if (!length(hits)) stop("Cannot identify ", label, "; available columns: ", paste(n, collapse = ", "),
                          ". Set the corresponding *_col argument.")
  if (length(hits) > 1L) stop("Ambiguous ", label, ": ", paste(hits, collapse = ", "),
                              ". Set the corresponding *_col argument.")
  hits[[1L]]
}

.adapt_require <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) stop(label, " missing columns: ", paste(miss, collapse = ", "))
}
.adapt_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
.adapt_dates <- function(x, label) {
  # ISO dates or ISO datetimes; preserve the date printed in the input (no TZ shift).
  d <- as.IDate(substr(as.character(x), 1L, 10L), format = "%Y-%m-%d")
  if (anyNA(d)) stop(label, " contains missing or non-ISO dates; check the input.")
  d
}
.adapt_calendar <- function(yrs, anchor_doy, length_days) {
  jan1 <- as.IDate(sprintf("%04d-01-01", yrs))
  next_jan1 <- as.IDate(sprintf("%04d-01-01", yrs + 1L))
  # Day 366 in a non-leap year is clamped to December 31.
  anchor <- jan1 + pmin(anchor_doy, as.integer(next_jan1 - jan1)) - 1L
  data.table(year = yrs, anchor_date = anchor,
             preseason_start = anchor - length_days, preseason_end = anchor - 1L)
}

.adapt_uhat_grid <- function(raw, template, alignment) {
  grid <- unique(raw[is.finite(lon) & is.finite(lat), .(lon, lat)])
  if (!nrow(grid)) stop("UHAT contains no valid coordinates")
  p <- terra::vect(as.data.frame(grid), geom = c("lon", "lat"), crs = "EPSG:4326")
  xy <- terra::crds(terra::project(p, terra::crs(template)))
  global_cell <- terra::cellFromXY(template, xy)
  keep <- !is.na(global_cell) & is.finite(xy[,1L]) & is.finite(xy[,2L])
  if (!any(keep)) stop("All UHAT coordinates fall outside the MODIS template")
  if (any(!keep)) warning(sum(!keep), " UHAT grid coordinates outside template; excluded")
  grid <- grid[keep]
  xy <- xy[keep, , drop = FALSE]
  rc <- terra::rowColFromCell(template, global_cell[keep])
  r0 <- min(rc[,1L]); r1 <- max(rc[,1L]); c0 <- min(rc[,2L]); c1 <- max(rc[,2L])
  nr <- r1 - r0 + 1L; nc <- c1 - c0 + 1L
  rx <- terra::res(template)[1L]; ry <- terra::res(template)[2L]
  xx <- terra::xmin(template) + (c0 - 1L) * rx
  yy <- terra::ymax(template) - r1 * ry
  # City-only raster: no full MODIS-sized array, and no intermediate file.
  base <- terra::rast(nrows = nr, ncols = nc, xmin = xx, xmax = xx + nc*rx,
                      ymin = yy, ymax = yy + nr*ry, crs = terra::crs(template))
  grid[, cell := as.numeric((rc[,1L] - r0) * nc + rc[,2L] - c0 + 1L)]
  if (alignment == "legacy") {
    # Same shift principle as build_UHAT_tif(); compute it ONCE using the
    # complete coordinate grid, not the non-NA extent of each annual window.
    dx <- min(xy[,1L]) - rx/2 - xx
    dy <- min(xy[,2L]) - ry/2 - yy
    base <- terra::shift(base, dx = dx, dy = dy)
    message(sprintf("UHAT grid shift: dx=%.3f, dy=%.3f (template CRS units)", dx, dy))
  }
  list(grid = grid, raster = base)
}

.adapt_daymet_window <- function(dm, ids, calendar, n_required) {
  # dm keyed by id; only access trees with this anchor DOY.
  z <- dm[.(ids), nomatch = 0L]
  if (!nrow(z)) return(data.table(id_key = character(), year = integer(),
                                  Daymet_srad_preseason = numeric(), Daymet_prcp_preseason = numeric(),
                                  Daymet_n_days_srad = integer(), Daymet_n_days_prcp = integer()))
  lookup <- rbindlist(lapply(seq_len(nrow(calendar)), function(i) {
    data.table(date = seq(calendar$preseason_start[i], calendar$preseason_end[i], by = "day"),
               year = calendar$year[i])
  }))
  z[, year := lookup$year[match(date, lookup$date)]]
  z <- z[!is.na(year)]
  if (anyDuplicated(z, by = c("id_key", "date")))
    stop("Duplicate Daymet id-date records within an adaptive window; resolve before summing precipitation")
  z[, {
    ns <- sum(is.finite(srad)); np <- sum(is.finite(prcp))
    list(Daymet_srad_preseason = if (ns >= n_required) .adapt_mean(srad) else NA_real_,
         Daymet_prcp_preseason = if (np >= n_required) sum(prcp[is.finite(prcp)]) else NA_real_,
         Daymet_n_days_srad = ns, Daymet_n_days_prcp = np)
  }, by = .(id_key, year)]
}

.adapt_uhat_window <- function(raw, base, pts, start, end, n_required) {
  dates <- as.IDate(seq(start, end, by = "day"))
  z <- raw[.(dates), nomatch = 0L] # keyed by date
  # Calculate separate means of daytime/nighttime predictions, consistent
  # with summarise_UHAT_seasons(), then convert Kelvin to Celsius.
  if (nrow(z) && anyDuplicated(z, by = c("lon", "lat", "date"))) {
    z <- z[, .(t2m_pred_day = .adapt_mean(t2m_pred_day),
               t2m_pred_night = .adapt_mean(t2m_pred_night)),
           by = .(lon, lat, cell, date)]
  }
  out <- data.table(UHAT_temp_preseason = rep(NA_real_, nrow(pts)),
                    UHAT_n_days_day = rep(NA_integer_, nrow(pts)),
                    UHAT_n_days_night = rep(NA_integer_, nrow(pts)))
  if (!nrow(z)) return(out)
  per_location <- z[, {
    nd <- sum(is.finite(t2m_pred_day)); nn <- sum(is.finite(t2m_pred_night))
    list(temp = if (nd >= n_required && nn >= n_required)
      (.adapt_mean(t2m_pred_day) + .adapt_mean(t2m_pred_night))/2 - 273.15 else NA_real_,
      n_day = nd, n_night = nn)
  }, by = .(lon, lat, cell)]
  # Duplicate projected locations in a cell are averaged as in build_UHAT_tif.
  # Coverage counts are the minimum among locations represented in the cell.
  per_cell <- per_location[, .(temp = .adapt_mean(temp),
                               n_day = min(n_day), n_night = min(n_night)), by = cell]
  r <- terra::rast(base)
  v <- rep(NA_real_, terra::ncell(r)); v[per_cell$cell] <- per_cell$temp
  terra::values(r) <- v
  names(r) <- "UHAT_temp_preseason"
  out[, UHAT_temp_preseason := terra::extract(r, pts, ID = FALSE)[[1L]]]
  # Extract coverage on the same aligned geometry, without extra rasters.
  cell <- terra::cellFromXY(r, terra::crds(pts))
  ii <- match(cell, per_cell$cell)
  out[, `:=`(UHAT_n_days_day = per_cell$n_day[ii], UHAT_n_days_night = per_cell$n_night[ii])]
  out
}

prepare_tree_Grid_climate_adaptive <- function(
    city, raw_uhat, daymet_raw, preseason_length = 30L, years = NULL,
    summary_path = NULL,
    summary_city_col = NULL, summary_species_col = NULL, summary_doy_col = NULL,
    tree_species_col = "taxa", min_coverage = 0.25, overwrite = FALSE,
    alignment = c("legacy", "template"),
    phenophase = c("SOS", "EOS")) {
  phenophase <- match.arg(phenophase)
  alignment <- match.arg(alignment)
  if (is.null(summary_path)) {
    summary_path <- paste0("~/lab-data/datasets/vegetation/PS/CONUS_urban/no_outlier/", "adapt_", phenophase, "_summary.rds")
  }
  if (!file.exists(summary_path)) stop("Missing ", phenophase, " summary: ", summary_path)
  mean_doy_name <- paste0("mean_", phenophase, "_DOY")
  if (length(preseason_length) != 1L || !is.numeric(preseason_length) ||
      !is.finite(preseason_length) || preseason_length != floor(preseason_length) ||
      preseason_length < 1 || preseason_length > 365) stop("preseason_length must be an integer in 1:365")
  preseason_length <- as.integer(preseason_length)
  if (length(min_coverage) != 1L || !is.numeric(min_coverage) || !is.finite(min_coverage) ||
      min_coverage <= 0 || min_coverage > 1) stop("min_coverage must be in (0,1]")
  n_required <- ceiling(preseason_length * min_coverage)
  out_dir <- here("data/proc/tree_climate", paste0("Grid_", phenophase, "_pre", preseason_length))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("Processing adaptive ", phenophase, " climate: ", city)
  
  # 1. One fixed reference phenophase DOY for each city-species pair.
  s <- readRDS(summary_path)
  if (inherits(s, "sf")) s <- sf::st_drop_geometry(s)
  s <- as.data.frame(s)
  cc <- .adapt_col(s, summary_city_col, c("city", "city_code", "city_name"), "summary city")
  sc <- .adapt_col(s, summary_species_col, c("species", "canonicalname", "scientific_name", "sci_name", "Factor"), "summary species")
  dc <- .adapt_col(s, summary_doy_col,
                   c("mean_doy", "doy_mean", mean_doy_name,
                     paste0("mean_", phenophase), paste0(phenophase, "_mean"),
                     paste0(phenophase, "_mean_doy"), "mean"),
                   paste("summary mean", phenophase, "DOY"))
  s <- s[!is.na(s[[cc]]) & trimws(as.character(s[[cc]])) == city, , drop = FALSE]
  ref <- unique(data.table(species_key = trimws(as.character(s[[sc]])),
                           mean_DOY = suppressWarnings(as.numeric(as.character(s[[dc]])))))
  ref <- ref[!is.na(species_key) & nzchar(species_key)]
  if (any(!is.finite(ref$mean_DOY) | ref$mean_DOY < 1 | ref$mean_DOY > 366))
    stop("Invalid mean ", phenophase, " DOY for ", city, "; check summary")
  if (anyDuplicated(ref$species_key)) stop("Multiple distinct mean ", phenophase, " DOYs for a city-species pair: ", city)
  if (!nrow(ref)) stop("No ", phenophase, " summary rows for ", city)
  ref[, anchor_DOY := as.integer(floor(mean_DOY + 0.5))]
  
  trees <- readRDS(paste0("~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/",city, ".rds"))
  if (inherits(trees, "sf")) trees <- sf::st_drop_geometry(trees)
  trees <- as.data.frame(trees)
  .adapt_require(trees, c("id", "lon", "lat"), "tree_list")
  # Resolve tree species field by matching against summary species names.
  if (is.null(tree_species_col)) {
    candidates <- names(trees)[gsub("[^a-z0-9]", "", tolower(names(trees))) %in%
                                 c("species", "canonicalname", "scientificname", "sciname", "speciesname", "factor", "taxon", "binomial")]
    scores <- vapply(candidates, function(n) sum(unique(trimws(as.character(trees[[n]]))) %in% ref$species_key), integer(1L))
    if (!length(scores) || max(scores) == 0L) stop("No species column matches ", phenophase, " summary; set tree_species_col. Columns: ", paste(names(trees), collapse = ", "))
    best <- candidates[scores == max(scores)]
    if (length(best) > 1L) {
      same <- vapply(best[-1L], function(n) identical(trimws(as.character(trees[[n]])), trimws(as.character(trees[[best[1L]]]))), logical(1L))
      if (!all(same)) stop("Ambiguous tree species columns: ", paste(best, collapse = ", "), "; set tree_species_col")
    }
    tree_species_col <- best[1L]
  }
  .adapt_require(trees, tree_species_col, "tree_list")
  trees$lon <- suppressWarnings(as.numeric(as.character(trees$lon)))
  trees$lat <- suppressWarnings(as.numeric(as.character(trees$lat)))
  valid <- is.finite(trees$lon) & is.finite(trees$lat) & abs(trees$lon) <= 180 & abs(trees$lat) <= 90 & !is.na(trees$id)
  trees <- trees[valid, , drop = FALSE]
  if (anyDuplicated(as.character(trees$id))) stop("Duplicate tree IDs in ", city)
  tree_species <- trimws(as.character(trees[[tree_species_col]]))
  match_ref <- match(tree_species, ref$species_key)
  if (anyNA(match_ref)) message(sum(is.na(match_ref)), " trees without matching species-city ", phenophase, " summary; excluded")
  trees <- trees[!is.na(match_ref), , drop = FALSE]
  tree_species <- tree_species[!is.na(match_ref)]
  match_ref <- match_ref[!is.na(match_ref)]
  if (!nrow(trees)) stop("No valid trees matched to ", phenophase, " summary for ", city)
  reserved <- c("year", "mean_SOS_DOY", "mean_EOS_DOY", "phenophase", "anchor_DOY", "anchor_date", "preseason_start", "preseason_end", "preseason_length",
                "UHAT_temp_preseason", "Daymet_srad_preseason", "Daymet_prcp_preseason", "Daymet_n_days_srad", "Daymet_n_days_prcp", "UHAT_n_days_day", "UHAT_n_days_night")
  if (length(intersect(names(trees), reserved))) stop("tree_list already contains output column(s): ", paste(intersect(names(trees), reserved), collapse = ", "))
  tree_sf <- sf::st_as_sf(trees, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  ids_all <- as.character(trees$id)
  anchors <- ref$anchor_DOY[match_ref]
  active_species <- unique(tree_species)
  safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
  paths <- setNames(file.path(out_dir, paste0(safe(city), "__", safe(active_species), ".rds")), active_species)
  if (anyDuplicated(tolower(unname(paths)))) stop("Species filenames collide after sanitization; revise species labels")
  pending <- if (overwrite) active_species else active_species[!file.exists(paths[active_species])]
  if (!length(pending)) { message("All species outputs already exist; skipped ", city); return(invisible(paths)) }
  message("Species to process: ", length(pending), "; summary columns: ", cc, ", ", sc, ", ", dc,
          "; tree species column: ", tree_species_col)
  
  # 2. Daily Daymet by tree ID. Read city data only once.
  dm <- as.data.table(daymet_raw)
  setnames(dm, "id", "id_key")
  dm[, date := .adapt_dates(date, "Daymet")]
  if (is.null(years)) years <- sort(unique(as.integer(format(unique(dm$date), "%Y"))))
  if (!is.numeric(years) || !length(years) || any(!is.finite(years)) || any(years != floor(years)) ||
      any(years < 1 | years > 9998)) stop("years must be a vector of valid integer years")
  years <- sort(unique(as.integer(years)))
  dm <- dm[id_key %in% ids_all[tree_species %in% pending]]
  setkey(dm, id_key)
  
  # 3. Read daily UHAT predictions and align the city template once.
  raw <- as.data.table(raw_uhat)
  raw[, date := .adapt_dates(time, "UHAT")]
  raw[, time := NULL]
  available_dates <- sort(unique(raw$date))
  if (length(available_dates) < 2L || median(as.integer(diff(available_dates))) > 1)
    stop("UHAT raw timestamps do not appear daily; daily adaptive windows require daily UHAT data")
  tfiles <- list.files("~/lab-data/datasets/climate/UHAT/grid",
                       pattern = "^MYD11A1\\..*LST_Day_1km\\.tif$", full.names = TRUE)
  if (!length(tfiles)) stop("No MODIS UHAT grid template found")
  grid_info <- .adapt_uhat_grid(raw, terra::rast(tfiles[1L]), alignment)
  raw[grid_info$grid, on = .(lon, lat), cell := i.cell]
  raw <- raw[!is.na(cell)]
  setkey(raw, date)
  pts_all <- terra::project(terra::vect(trees[, c("lon", "lat"), drop = FALSE],
                                        geom = c("lon", "lat"), crs = "EPSG:4326"),
                            terra::crs(grid_info$raster))
  
  # 4. Reuse each year-window raster for every species sharing an anchor DOY.
  for (doy in sort(unique(anchors[tree_species %in% pending]))) {
    idx <- which(anchors == doy & tree_species %in% pending)
    ids <- ids_all[idx]
    spp <- tree_species[idx]
    cal <- .adapt_calendar(years, doy, preseason_length)
    message("Anchor DOY ", doy, ": ", length(unique(spp)), " species, ", length(idx), " trees")
    dm_stats <- .adapt_daymet_window(dm, ids, cal, n_required)
    ans <- rbindlist(lapply(seq_len(nrow(cal)), function(i) {
      message("  ", cal$year[i], ": ", cal$preseason_start[i], " to ", cal$preseason_end[i])
      u <- .adapt_uhat_window(raw, grid_info$raster, pts_all[idx, ],
                              cal$preseason_start[i], cal$preseason_end[i], n_required)
      cbind(data.table(id_key = ids, year = cal$year[i], anchor_DOY = doy,
                       anchor_date = cal$anchor_date[i], preseason_start = cal$preseason_start[i],
                       preseason_end = cal$preseason_end[i], preseason_length = preseason_length), u)
    }))
    ans <- merge(ans, dm_stats, by = c("id_key", "year"), all.x = TRUE, sort = FALSE)
    for (nm in c("Daymet_n_days_srad", "Daymet_n_days_prcp")) set(ans, which(is.na(ans[[nm]])), nm, 0L)
    # Attach original tree attributes and geometry; split into city-species RDS files.
    for (sp in unique(spp)) {
      j <- idx[spp == sp]
      tree_sp <- tree_sf[j, ]
      tab <- ans[id_key %in% ids_all[j]]
      ii <- match(tab$id_key, ids_all[j])
      out <- tree_sp[ii, ]
      out[[mean_doy_name]] <- ref$mean_DOY[match(sp, ref$species_key)]
      out$phenophase <- phenophase
      for (nm in setdiff(names(tab), "id_key")) out[[nm]] <- tab[[nm]]
      out <- out[order(match(as.character(out$id), ids_all[j]), out$year), ]
      target <- unname(paths[sp])
      # Write to a temporary RDS in the same directory, then rename.
      tmp <- tempfile(pattern = ".adaptive_", tmpdir = out_dir, fileext = ".rds")
      saveRDS(out, tmp, compress = "gzip")
      if (!file.rename(tmp, target)) { unlink(tmp); stop("Cannot finalize output: ", target) }
      message("  Saved ", basename(target), " (", nrow(out), " rows; missing UHAT=", sum(is.na(out$UHAT_temp_preseason)), ")")
    }
  }
  invisible(paths)
}
