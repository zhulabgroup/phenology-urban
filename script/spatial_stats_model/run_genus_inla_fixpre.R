

prepare_genus_data <- function(phe_type, points_in_buffer, all_sites_temp, srad_daily, source) {
  if (source == "WU") {
    ### WU sites
    points_in_buffer_joined <- points_in_buffer %>%
      left_join((srad_daily %>%
                   select(id, year, season, Avg_srad) %>%
                   pivot_wider(names_from = season, values_from = Avg_srad, names_prefix = "srad_")),
                by = c("id", "year")) %>%
      left_join(all_sites_temp %>%
                  select(name, year, season, Avg_AvgTemp) %>%
                  pivot_wider(names_from = season, values_from = Avg_AvgTemp, names_prefix = "temp_"),
                by = c("Site" = "name", "year")) %>%
      left_join(all_sites_temp %>%
                  select(name, year, season, Sum_Sum_mm) %>%
                  pivot_wider(names_from = season, values_from = Sum_Sum_mm, names_prefix = "prcp_"),
                by = c("Site" = "name", "year"))
  } else if (source == "Daymet") {
    ## All Daymet
    points_in_buffer_joined <- points_in_buffer %>%
      left_join((srad_daily %>%
                   select(id, year, season, Avg_srad) %>%
                   pivot_wider(names_from = season, values_from = Avg_srad, names_prefix = "srad_")),
                by = c("id", "year")) %>%
      left_join(all_sites_temp %>%
                  mutate(Avg_AvgTemp = temp) %>%
                  select(id, year, season, Avg_AvgTemp) %>%
                  pivot_wider(names_from = season, values_from = Avg_AvgTemp, names_prefix = "temp_"),
                by = c("id", "year")) %>%
      left_join(all_sites_temp %>%
                  select(id, year, season, prcp) %>%
                  pivot_wider(names_from = season, values_from = prcp, names_prefix = "prcp_"),
                by = c("id", "year"))
  } else {
    stop("Invalid source. Must be 'WU' or 'Daymet'.")
  }
  
  if (phe_type == "SOS") {
    phe_metric <- "doy"
    direction <- "up"
    thres <- 0.5
    model_input <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, taxa,genus, lon, lat, year, doy, hex_id, ends_with("winter")) %>%
      drop_na()
  } else if (phe_type == "EOS") {
    phe_metric <- "doy"
    direction <- "down"
    thres <- 0.5
    model_input <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres == !!thres) %>%
      select(id, taxa, genus, lon, lat, year, doy, hex_id, ends_with("summer")) %>%
      drop_na()
  } else if (phe_type == "GreenUp") {
    phe_metric <- "doy_diff"
    direction <- "up"
    diff <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = points_in_buffer_joined %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, taxa, genus, lon, lat, year, doy_diff, hex_id, ends_with("spring")) %>%
      drop_na()
    
  } else {
    direction <- "down"
    diff <- points_in_buffer_joined %>%
      filter(direction == !!direction & thres %in% c(0.2, 0.8)) %>%
      group_by(id, year) %>%
      summarise(
        doy_diff = abs(diff(doy[order(thres)])),
        .groups = "drop"
      )
    model_input = points_in_buffer_joined %>%
      filter(direction == !!direction & thres == 0.2) %>%
      left_join(diff, by = c("year", "id")) %>%
      select(id, taxa, genus, lon, lat, year, doy_diff, hex_id, ends_with("fall")) %>%
      drop_na()
  }
  
  return(model_input)
}



run_inla_genus = function(model_input_all, r_output_dir, phe_type, city){
  # options(repos = c(
  #   INLA = "https://inla.r-inla-download.org/R/testing",
  #   getOption("repos")
  # ))
  # install.packages("inlabru")
  # install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  
  genus_list <- unique(model_input_all$genus)
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
  
  for (genus in genus_list) {
    
    model_path <- paste0(r_output_dir, genus, "_", phe_type, ".rds")
    
    if (file.exists(model_path)) {
      cat(sprintf("\n File already exists for genus: %s, phe_type: %s, city: %s. Skipping INLA run.\n", 
                  genus, phe_type, city))
      next  # Skip the current iteration of the loop
    } # if the file exists, skip to the next iteration  
    
    cat(sprintf("\n Running INLA for genus: %s, phe_type: %s, city: %s\n", 
                genus, phe_type, city))
    if (phe_type %in% c("SOS","EOS")) {
      phe_metric <- "doy"
    } else {
      phe_metric <- "doy_diff"
    }
    
    
    
    model_input_sf <- model_input_all %>%
      filter(genus == !!genus) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      st_transform(crs = st_crs(boundary)) %>%
      mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2],
             year = as.factor(year)) 
    
    model_input <- model_input_sf %>% 
      st_drop_geometry()
    
    
    ### Prediction (In-sample)
    prediction <- model_input
    coorp <- cbind(prediction$lon, prediction$lat) 
    
    # Load required libraries
    ######################################## Prepare the mesh and spde ######################################## 
    
    ### Build the mesh 
    coord <- cbind(model_input$lon, model_input$lat) # coordinates of observations
    domain <- inla.nonconvex.hull(as.matrix(coord),
                                  concave = -.07,convex = -0.05, resolution=c(100,100))
    
    mesh <- inla.mesh.2d(boundary = domain, loc = coord,
                         max.edge = c(10000, 20000), 
                         cutoff = 1000)  # prevent close points from being connected
    
    ### Build the mesh 
    coord <- cbind(model_input$lon, model_input$lat) # coordinates of observations
    domain <- inla.nonconvex.hull(as.matrix(coord),
                                  concave = -.07,convex = -0.05, resolution=c(100,100))
    
    mesh <- inla.mesh.2d(boundary = domain, loc = coord,
                         max.edge = c(10000, 20000), 
                         cutoff = 1000)  # prevent close points from being connected
    
    ### Build the spde 
    
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
    A.phe <- inla.spde.make.A(mesh = mesh, loc = coord)
    s.index <- inla.spde.make.index(name = "spatial.field",
                                    n.spde = spde$n.spde)
    
    A.phe.p <- inla.spde.make.A(mesh = mesh, loc = coorp)
    
    ######################################## Stack the data ########################################
    
    stk.e <- inla.stack(
      tag = "est",
      data = list(y = model_input[[phe_metric]]),
      A = list(1, A.phe),
      effects = list(data.frame(Intercept = rep(1, nrow(coord)),
                                temp = (model_input %>% select(starts_with("temp")) %>% pull()),
                                prcp = (model_input %>% select(starts_with("prcp")) %>% pull()),
                                srad = (model_input %>% select(starts_with("srad")) %>% pull()),
                                year = model_input$year), 
                     s = s.index)
    )
    
    
    stk.p <- inla.stack(
      tag = "pred",
      data = list(y = NA),
      A = list(1, A.phe.p),
      effects = list(data.frame(Intercept = rep(1, nrow(coorp)),
                                temp = (prediction %>% select(starts_with("temp")) %>% pull()),
                                prcp = (prediction %>% select(starts_with("prcp")) %>% pull()),
                                srad = (prediction %>% select(starts_with("srad")) %>% pull()),
                                year = prediction$year), 
                     s = s.index)
    )
    
    # stk.full has stk.e and stk.p
    stk.full <- inla.stack(stk.e, stk.p)
    
    ######################################## Run/Save the model ########################################
    
    ### Run the model
    formula <- y ~ 0 + Intercept + temp + prcp + srad + f(year, model = 'iid') + f(spatial.field, model = spde)
    
    res <- inla(formula,
                family = "gaussian",
                data = inla.stack.data(stk.full),
                control.compute = list(return.marginals = TRUE, waic = TRUE),
                control.predictor = list(
                  compute = TRUE,
                  A = inla.stack.A(stk.full)))
    # res$waic$waic
    # summary(res)
    
    fix_effect <- as.data.frame(summary(res)$fixed)%>%
      rename(Estimate = mean,Lower = `0.025quant`,Upper = `0.975quant`) %>%
      mutate(Sig = ifelse(Lower*Upper>0, "*", "")) 
    fix_effect <- data.frame(
      Factor = rownames(fix_effect),  # Add Factor column
      fix_effect, 
      row.names = NULL  # Remove row names to avoid duplication
    ) %>%
      filter(Factor == "temp") %>%
      mutate(Factor = ifelse(Factor == "temp", genus, Factor))
    
    index <- inla.stack.index(stack = stk.full, tag = "pred")$data
    pred_mean <- res$summary.fitted.values$mean[index]
    
    compare<- rbind(
      Lon = coorp[, 1], Lat = coorp[, 2],
      value = pred_mean, variable = "pred_mean", true_value = model_input[[phe_metric]]
    )
    
    ######################################## Check autocorrelation ###########
    
    spatial_field_mean <- res$summary.random$spatial.field$mean
    
    fitted_spatial <- as.vector(A.phe %*% spatial_field_mean)
    temp_vec <- model_input %>% select(starts_with("temp")) %>% pull()
    
    cor_temp_spatial <- cor(temp_vec, fitted_spatial, use = "complete.obs")
    print(paste("Correlation between temp and spatial random effect:", round(cor_temp_spatial, 3)))
    ggplot(data.frame(temp = temp_vec, spatial = fitted_spatial), aes(x = temp, y = spatial)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      labs(title = paste("Spatial Confounding: cor =", round(cor_temp_spatial, 2)),
           x = "Temperature", y = "Spatial Random Effect") +
      theme_minimal()
    
    ######################################## Compare with non-spatial model ########################################
    formula_nosp <- y ~ 0 + Intercept + temp + prcp + srad + f(year, model = 'iid')
    
    res_nosp <- inla(formula_nosp,
                     family = "gaussian",
                     data = inla.stack.data(stk.full), 
                     control.predictor = list(
                       compute = TRUE,
                       A = inla.stack.A(stk.full)
                     ),
                     control.compute = list(dic = TRUE, waic = TRUE))
    # summary(res_nosp)
    # Extract variance components
    var_fixed_with_spatial <- res$summary.fixed$mean^2 %>% sum()
    var_fixed_nospatial <- res_nosp$summary.fixed$mean^2 %>% sum()
    
    index_est <- inla.stack.index(stack = stk.full, tag = "est")$data
    fitted_total <- res$summary.fitted.values$mean[index_est]
    y_obs <- model_input[[phe_metric]]
    var_y <- var(y_obs, na.rm = TRUE)
    residuals <- y_obs - fitted_total
    var_resid <- var(residuals, na.rm = TRUE)
    R2_conditional <- 1 - var_resid / var_y
    
    
    fitted_fixed_only <- res_nosp$summary.fitted.values$mean[index_est]
    resid_fixed <- y_obs - fitted_fixed_only
    var_resid_fixed <- var(resid_fixed, na.rm = TRUE)
    R2_marginal <- 1 - var_resid_fixed / var_y
    
    cat("Marginal R² (fixed only):", round(R2_marginal, 3), "\n")
    cat("Conditional R² (fixed + random):", round(R2_conditional, 3), "\n")
    
    resid_df <- data.frame(
      lon = model_input$lon,
      lat = model_input$lat,
      resid = resid_fixed
    )
    
    library(ggplot2)
    ggplot(resid_df, aes(x = lon, y = lat, color = resid)) +
      geom_point() +
      scale_color_gradient2() +
      coord_equal() +
      labs(title = "Residuals from non-spatial model",
           fill = "Residual") +
      theme_minimal()
    
    ######################################## Save the model results ########################################
    
    list(mesh = mesh,
         res = res,
         spde = spde,
         stk.full = stk.full,
         fix_effect = fix_effect,
         compare = compare,
         genus = genus) %>%
      saveRDS(model_path)
  }
}


plot_coef_genus_city <- function(r_output_dir, phe_type){
  pattern = paste0(phe_type,"\\.rds$")
  files <- list.files(r_output_dir, pattern = pattern, full.names = TRUE)
  plots_tavg <- list()
  coef_summary <- data.frame()
  plots_mean <- list()
  plots_sd <- list()
  plots_semi <- list()
  
  # for one genus
  for (file in files) {
    model_data <- readRDS(file)
    genus <- model_data$genus
    # beta_names <- c("beta[0]", "beta[AvgTemp_mean]", "beta[Sum_mm_sum]", "beta[srad_mean]")
    # df_list <- lapply(seq_along(model_data$res$marginals.fixed), function(i) {
    #   data.frame(x = model_data$res$marginals.fixed[[i]][,1], 
    #              density = model_data$res$marginals.fixed[[i]][,2], 
    #              beta = beta_names[i])
    # })
    # df <- bind_rows(df_list)
    beta_avg_temp <- model_data$res$marginals.fixed[[3]]
    df <- data.frame(
      x = beta_avg_temp[, 1],
      density = beta_avg_temp[, 2],
      beta = "beta[temp]"
    )
    plot_tavg <- ggplot(df, aes(x = x, y = density)) +
      geom_line(linewidth = 1, color = "blue") +
      # facet_wrap(~beta, scales = "free", labeller = label_parsed) + 
      theme_minimal() +
      labs(x = "Coefficient Value", y = "Density") +
      ggtitle(paste(genus)) +
      theme(legend.position = "none",
            plot.title = element_text(size = 8))
    
    plots_tavg[[genus]] <- plot_tavg
    
    coef_summary <- coef_summary %>%
      rbind(model_data$fix_effect)
    
    
    proj <- inla.mesh.projector(model_data$mesh,dims = c(300, 300))
    mean_s <- inla.mesh.project(proj, model_data$res$summary.random$s$mean)
    sd_s <- inla.mesh.project(proj, model_data$res$summary.random$s$sd)
    df <- expand.grid(x = proj$x, y = proj$y)
    df$mean_s <- as.vector(mean_s)
    df$sd_s <- as.vector(sd_s)
    
    p_mean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) +
      geom_raster() +
      scale_fill_viridis(na.value = "transparent") +
      ggtitle(paste(genus)) +
      coord_fixed(ratio = 0.8) + theme_bw() +
      theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
            plot.title = element_text(size = 8))
    
    p_sd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) +
      geom_raster() +
      scale_fill_viridis(na.value = "transparent") +
      ggtitle(paste(genus)) +
      coord_fixed(ratio = 1) + theme_bw()+
      theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
            plot.title = element_text(size = 8))
    
    #Plot the decay of spatial autocorrelation across a user-defined range
    
    spde.est <- inla.spde2.result(inla = model_data$res, name = "spatial.field",
                                  spde = model_data$spde, do.transf = TRUE)
    Kappa = inla.zmarginal(spde.est$marginals.kappa[[1]], silent = TRUE)[1]$mean
    MaxRange <- 50000
    Resolution <- 100
    d.vec <- seq(0, MaxRange, length = Resolution)
    Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1)
    Cor.M[1] <- 1
    CorData <- data.frame(Distance = d.vec, Correlation = Cor.M)
    
    p_semi <- ggplot(CorData, aes(x = Distance, y = Correlation)) +
      geom_line(linewidth = 1, alpha = 0.8, color = "blue") +
      coord_fixed(ratio = MaxRange) +
      labs(x = "Distance (m)", y = "Correlation", title = "Matern Correlation Function") +
      ggtitle(paste(genus)) +
      theme_minimal()+
      theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
            plot.title = element_text(size = 8))
    
    plots_mean[[genus]] <- p_mean
    plots_sd[[genus]] <- p_sd
    plots_semi[[genus]] <- p_semi
  }
  
  combined_plot_tavg <- wrap_plots(plots_tavg, ncol = 5)
  combined_plot_semi <- wrap_plots(plots_semi, ncol = 5)
  combined_plot_mean <- wrap_plots(plots_mean, ncol = 5)
  combined_plot_sd <- wrap_plots(plots_sd, ncol = 5)
  
  coef_summary <- coef_summary %>%
    mutate(SigAlpha = ifelse(Sig == "*", 1, 0.2))%>%
    filter(Factor != "Alnus")
  
  p_coef <-  ggplot(coef_summary, aes(x = factor(Factor, levels = rev(Factor)), 
                                      y = Estimate)) +
    geom_point(aes(alpha = SigAlpha), position = position_dodge(w = 0.5), size = 3) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper, alpha = SigAlpha), 
                  position = position_dodge(w = 0.5), 
                  width = 0.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = NULL, y = "Estimate") +
    coord_flip() +
    scale_alpha_identity() +
    theme(legend.position = "right") +
    geom_text(aes(label = Sig, y = Upper + 0.05  ), hjust = -0.5, show.legend = FALSE) +
    geom_text(aes(label = Estimate, y = Estimate  ), vjust = -0.7, show.legend = FALSE) +
    guides(alpha = "none") +
    ggtitle(bquote(beta[Temp]~"for"~.(phe_type))) +
    theme_minimal()
  
  return(list(
    combined_plot_tavg = combined_plot_tavg,
    combined_plot_semi = combined_plot_semi,
    combined_plot_mean = combined_plot_mean,
    combined_plot_sd = combined_plot_sd,
    coef_summary_plot = p_coef
  ))
}

plot_coef_genus_cities <- function(r_output_dir, phe_type){
  pattern = paste0(phe_type,"\\.rds$")
  files <- list.files(r_output_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  coef_summary <- data.frame()
  
  # for one 
  for (file in files) {
    model_data <- readRDS(file)
    genus <- model_data$genus
    rel_path <- str_remove(file, paste0("^", r_output_dir))
    city <- str_split(rel_path, "/")[[1]][12]
    
    coef_summary <- coef_summary %>%
      rbind(model_data$fix_effect %>% mutate(city = city)) 
  }
  
  coef_summary <- coef_summary %>%
    mutate(SigAlpha = ifelse(Sig == "*", 1, 0.2)) %>%
    filter(Factor != "Alnus")
  
  
  p_coef <-  ggplot(coef_summary, aes(x = factor(Factor, levels = rev(unique(Factor))), 
                                      y = Estimate,
                                      color = city)) +
    geom_point(aes(alpha = SigAlpha), position = position_dodge(w = 0.5), size = 1.2) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper, alpha = SigAlpha), 
                  position = position_dodge(w = 0.5), 
                  linewidth = 0.5,
                  width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha =0.5) +
    # coord_flip() +
    labs(x = NULL, y = "Estimate") +
    scale_alpha_identity() +
    scale_color_manual(
      values = c("NY" = "darkolivegreen4", "HT" = "navyblue", "DV" = "chocolate", "ST" = "purple")
    ) +
    theme(legend.position = "right") +
    # geom_text(aes(label = Sig, y = Upper + 0.05), hjust = -0.5, show.legend = FALSE, size = 3) +
    # geom_text(aes(label = Estimate, y = Estimate), vjust = -0.9, show.legend = FALSE, size = 3) +
    guides(alpha = "none") +
    theme_minimal()+
    ggtitle(bquote(beta[Temp]~"for"~.(phe_type))) +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      strip.text = element_text(size = 12, face = "bold"))
  
  return(p_coef)
}

# 
# 
# 
# mesh$n
# ggplot() +
#   geom_sf(data=boundary,color='turquoise',fill='transparent') +
#   gg(mesh)
#   # geom_sf(data=tavg_allyear_var_sf,col='purple',size=1.7,alpha=0.5)
# 
# 
# # dim(A.phe)
# 
# ######################################## Inference from model ########################################
# 
# ########################################
# # Spatial field
# ########################################
# # projection matrix for 300 x 300 locations on a grid that covers the mesh region.
# 
# model_result <- readRDS(model_path)
# 
# proj <- inla.mesh.projector(model_result$mesh,dims = c(300, 300))
# mean_s <- inla.mesh.project(proj, model_result$res$summary.random$s$mean)
# sd_s <- inla.mesh.project(proj, model_result$res$summary.random$s$sd)
# df <- expand.grid(x = proj$x, y = proj$y)
# df$mean_s <- as.vector(mean_s)
# df$sd_s <- as.vector(sd_s)
# 
# gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) +
#   geom_raster() +
#   scale_fill_viridis(na.value = "transparent") +
#   coord_fixed(ratio = 1) + theme_bw()
# 
# gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) +
#   geom_raster() +
#   scale_fill_viridis(na.value = "transparent") +
#   coord_fixed(ratio = 1) + theme_bw()
# 
# gmean + gsd
# 
# #Plot the decay of spatial autocorrelation across a user-defined range
# 
# spde.est <- inla.spde2.result(inla = model_result$res, name = "spatial.field",
#                               spde = model_result$spde, do.transf = TRUE)
# Kappa = inla.zmarginal(spde.est$marginals.kappa[[1]])[1]$mean
# MaxRange <- 50000
# Resolution <- 100
# d.vec <- seq(0, MaxRange, length = Resolution)
# Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1)
# Cor.M[1] <- 1
# CorData <- data.frame(Distance = d.vec, Correlation = Cor.M)
# ggplot(CorData, aes(x = Distance, y = Correlation)) +
#   geom_line(linewidth = 1, alpha = 0.8, color = "blue") +
#   coord_fixed(ratio = MaxRange) +
#   labs(x = "Distance (m)", y = "Correlation", title = "Matern Correlation Function") +
#   theme_minimal()
# 
# ########################################
# # Random year effect
# ########################################
# 
# model_result$res$summary.random$year
# 
# ########################################
# # Fixed effect
# ########################################
# 
# # for one genus
# 
# beta_names <- c("beta[0]", "beta[temp]", "beta[prcp]", "beta[srad]")
# df_list <- lapply(seq_along(model_result$res$marginals.fixed), function(i) {
#   data.frame(x = model_result$res$marginals.fixed[[i]][,1], 
#              density = model_result$res$marginals.fixed[[i]][,2], 
#              beta = beta_names[i])
# })
# df <- bind_rows(df_list)
# ggplot(df, aes(x = x, y = density)) +
#   geom_line(linewidth = 1, color = "blue") +
#   facet_wrap(~beta, scales = "free", labeller = label_parsed) + 
#   theme_minimal() +
#   labs(x = "Coefficient Value", y = "Density") +
#   theme(legend.position = "none")  # Remove legend since we use facets
# 
# # for comparison among genus
# 
# ggplot(model_result$fix_effect, aes(x = factor(Factor, levels = rev(Factor)), 
#                                         y = Estimate, 
#                                         colour = Factor)) +
#   geom_point(alpha = SigAlpha, position = position_dodge(w = 0.5), size = 3) +
#   geom_errorbar(aes(ymin = Lower, ymax = Upper, alpha = SigAlpha), 
#                 position = position_dodge(w = 0.5), 
#                 width = 0.2) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   labs(x = NULL, y = "Estimate") +
#   coord_flip() +
#   scale_alpha_identity() +
#   theme(legend.position = "right") +
#   geom_text(aes(label = Sig, y = Upper + 0.05  ), hjust = -0.5, show.legend = FALSE) +
#   guides(alpha = "none") +
#   theme_minimal()
# 
# ######################################## Prediction from model ########################################
# 
### In sample - Prediction
# 
# compare <- data.frame(t(compare)) %>%
#   mutate(value = as.numeric(value), true_value = as.numeric(true_value))
# compare <- compare %>% mutate(diff = value - true_value)
# lm_fit <- lm(value ~ true_value, data = compare)
# 
# slope <- round(coef(lm_fit)[2], 3)
# 
# ggplot(compare, aes(x = true_value, y = value)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
#   geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue") +
#   annotate("text", x = -Inf, y = Inf,
#            label = paste("Slope =", round(slope, 2)),
#            hjust = -0.1, vjust = 1.5,
#            color = "blue", size = 5) +
#   labs(x = "True value", y = "Predicted value") +
#   theme_classic()
# # spatial
# 
# ggplot(compare, 
#        aes(x = Lat, y = Lon, color=diff)) +
#   geom_point(alpha = 0.6) 
#   theme_minimal()
