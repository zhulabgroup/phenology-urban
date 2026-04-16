fit_uni <- function(phe_type, input_data, city_feature, species_name, clim = "mean_temp", bootstrap){
  combine_feature <- input_data |>
    filter(Factor == species_name) |>
    mutate(w = 1/(sd^2)) |>
    left_join(
      city_feature,
      by = "city"
    )
  
  df_phe <- combine_feature |>
    filter(phe_type == !!phe_type)
  
  if (nrow(df_phe) < 2) return(NULL)
  
  if (bootstrap == F){
    # lm
    fmla <- as.formula(paste0("Estimate ~ ", clim))
    fit_w <- lm(fmla, data = df_phe, weights = w)
    
    summary_fit <- summary(fit_w)
    R2   <- summary_fit$r.squared
    pval <- summary_fit$coefficients[2, 4]
    slope <- summary_fit$coefficients[2, 1]
    label_text <- sprintf("italic(slope) == %.3f*','~italic(p) == %.3f", slope, pval)
    
    # interval
    newdat <- data.frame(
      clim_var = seq(min(df_phe[[clim]], na.rm = TRUE),
                     max(df_phe[[clim]], na.rm = TRUE),
                     length.out = 200)
    )
    names(newdat) <- clim
    
    pred <- predict(fit_w, newdata = newdat, se.fit = TRUE)
    newdat$y_pred <- pred$fit
    newdat$y_lwr  <- pred$fit - 1.96 * pred$se.fit
    newdat$y_upr  <- pred$fit + 1.96 * pred$se.fit
    newdat$species  <- species_name
    
    # plot
    p <- ggplot(df_phe, aes(x = .data[[clim]], y = Estimate, alpha = SigAlpha)) +
      geom_pointrange(aes(ymin = Lower, ymax = Upper, colour = Color)) +
      geom_ribbon(data = newdat,
                  aes(x = .data[[clim]], ymin = y_lwr, ymax = y_upr),
                  inherit.aes = FALSE,
                  fill = "grey80", alpha = 0.4) +
      geom_line(data = newdat,
                aes(x = .data[[clim]], y = y_pred),
                inherit.aes = FALSE,
                colour = "grey40", linewidth = 0.8) +
      annotate("text", x = Inf, y = Inf,
               label = label_text, parse = TRUE,
               hjust = 1.2, vjust = 1.5, size = 4) +
      scale_alpha_identity() +
      scale_colour_identity() +
      guides(alpha = "none") +
      theme_minimal() +
      labs(
        title = species_name,
        x = case_when(
          (rlang::as_name(clim) %in% c("mean_prcp", "bio12")) ~ "Mean annual precipitation (mm)",
          (rlang::as_name(clim) %in% c("mean_temp", "bio1")) ~ "Mean annual temperature (°C)",
          (rlang::as_name(clim) %in% c("hot_edge")) ~ "Safety margin for temperature (°C)",
          (rlang::as_name(clim) %in% c("dry_edge")) ~ "Safety margin for precipitation (mm)",
          (rlang::as_name(clim) %in% c("d_temp")) ~ "Difference between city's temperature and niche centroid (°C)",
          (rlang::as_name(clim) %in% c("d_prcp")) ~ "Difference between city's precipitation and niche centroid (mm)"
        ),
        y = bquote(beta[Temp]~" estimate for "~.(phe_type))
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "italic"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
        # axis.title = element_blank()
      )
    plot <- p
    model_summary <- summary_fit
    pred_table <- newdat
    full_data <- NA
    
  } else {
    set.seed(1)
    B <- 500
    boot_fits <- map_dfr(seq_len(B), function(b) {
      dat_r <- df_phe |>
        mutate(beta_draw = rnorm(n(), mean = Estimate, sd = sd))
      fmla <- as.formula(paste0("beta_draw ~ ", clim))
      fit <- lm(fmla, data = dat_r)
      coef_tab <- summary(fit)$coefficients
      tibble(
        boot = b,
        term = rownames(coef_tab),
        estimate = coef_tab[, "Estimate"],
        se       = coef_tab[, "Std. Error"],
        p        = coef_tab[, "Pr(>|t|)"]
      )
    })
    boot_summary <- boot_fits %>%
      filter(!is.na(term)) %>%
      group_by(term) %>%
      summarise(
        mean = mean(estimate, na.rm = TRUE),
        sd = sd(estimate, na.rm = TRUE),
        q025 = quantile(estimate, 0.025, na.rm = TRUE),
        q50  = quantile(estimate, 0.50,  na.rm = TRUE),
        q975 = quantile(estimate, 0.975, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(species_name = species_name,
             phe_type = phe_type,
             size = nrow(df_phe))
    model_summary <- boot_summary
    full_data <- boot_fits
    
    # plot
    coef_df <- boot_fits %>%
      filter(term %in% c("(Intercept)", clim)) %>%
      select(boot, term, estimate) %>%
      tidyr::pivot_wider(names_from = term, values_from = estimate) %>%
      rename(b0 = `(Intercept)`, b1 = !!clim)
    xgrid <- tibble(x = seq(min(df_phe[[clim]], na.rm = TRUE),
                            max(df_phe[[clim]], na.rm = TRUE),
                            length.out = 200))
    pred_band <- coef_df %>%
      tidyr::crossing(xgrid) %>%
      mutate(yhat = b0 + b1 * x) %>%
      group_by(x) %>%
      summarise(
        y_pred = median(yhat, na.rm = TRUE),
        y_lwr = quantile(yhat, 0.025, na.rm = TRUE),
        y_upr = quantile(yhat, 0.975, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(species = species_name) |>
      rename(!!clim := x)
    
    label_text <- sprintf(
      "Slope = %.3f [%.3f, %.3f]",
      boot_summary$mean[boot_summary$term == clim],
      boot_summary$q025[boot_summary$term == clim],
      boot_summary$q975[boot_summary$term == clim]
    )
    
    p <- ggplot(df_phe, aes(x = .data[[clim]], y = Estimate, alpha = SigAlpha)) +
      geom_pointrange(aes(ymin = Lower, ymax = Upper, colour = Color)) +
      
      geom_ribbon(data = pred_band, aes(x = .data[[clim]], ymin = y_lwr, ymax = y_upr),
                  inherit.aes = FALSE, 
                  fill = "grey80", alpha = 0.4) +
      geom_line(data = pred_band, aes(x = .data[[clim]], y = y_pred),
                inherit.aes = FALSE, 
                colour = "grey40", linewidth = 0.8) +
      annotate("text", x = Inf, y = Inf,
               label = label_text, 
               hjust = 1.2, vjust = 1.5, size = 4) +
      scale_alpha_identity() +
      scale_colour_identity() +
      guides(alpha = "none") +
      theme_minimal() +
      labs(
        title = species_name,
        x = case_when(
          (rlang::as_name(clim) %in% c("mean_prcp", "bio12")) ~ "Mean annual precipitation (mm)",
          (rlang::as_name(clim) %in% c("mean_temp", "bio1")) ~ "Mean annual temperature (°C)",
          (rlang::as_name(clim) %in% c("hot_edge")) ~ "Safety margin for temperature (°C)",
          (rlang::as_name(clim) %in% c("dry_edge")) ~ "Safety margin for precipitation (mm)",
          (rlang::as_name(clim) %in% c("d_temp")) ~ "Difference between city's temperature and niche centroid (°C)",
          (rlang::as_name(clim) %in% c("d_prcp")) ~ "Difference between city's precipitation and niche centroid (mm)"
        ),
        y = bquote(beta[Temp]~" estimate for "~.(phe_type))
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "italic"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
        # axis.title = element_blank()
      )
    plot <- p
    pred_table <- pred_band
  }
  
  return(list(model_summary = model_summary, 
              plot = plot,
              pred_table = pred_table,
              full_data = full_data))
}


fit_multi_pglmm <- function(phe_type, input_data, city_feature, package, bootstrap, prcp, temp) {
  D <- input_data %>%
    rename(species = Factor) %>%
    mutate(species = gsub("_", " ", species)) %>%
    filter(phe_type == !!phe_type,
           species %in% phylo_tree$tip.label) %>%
    left_join(city_feature, by = "city")
  
  sp_data <- sort(unique(na.omit(D$species)))
  sp_tree <- sort(phylo_tree$tip.label)
  keep_sp <- intersect(sp_data, sp_tree)
  
  D2 <- D %>%
    filter(species %in% keep_sp) %>%
    droplevels() %>%
    mutate(
      species = factor(species, levels = keep_sp),
      city = factor(city)
    )
  
  tr2 <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, keep_sp))
  A2 <- ape::vcv(tr2, corr = TRUE)
  
  D2 <- D2 %>%
    filter(species %in% rownames(A2)) %>%
    mutate(
      species  = factor(species),
      prcp_s  = as.numeric(scale(.data[[prcp]])), # standardize the climate variables prcp
      temp_s   = as.numeric(scale(.data[[temp]])), # standardize the climate variables temp
      w        = 1/(sd^2)
    ) %>%
    droplevels()
  
  A2 <- A2[levels(D2$species), levels(D2$species)]
  
  if(package == "phyr"){
    
    if(bootstrap == TRUE){
      set.seed(1)
      B <- 500
      boot_fits <- map_dfr(seq_len(B), function(b) {
        dat_r <- D2 %>%
          mutate(beta_draw = rnorm(n(), mean = Estimate, sd = sd))
        fit <- pglmm(
          beta_draw ~ 1 + temp_s + prcp_s + (1 | species__),
          data = dat_r,
          family = "gaussian",
          cov_ranef = list(species = A2),
          bayes = FALSE,
          REML = FALSE
        )
        fit_nophy <- pglmm(
          beta_draw ~ 1 + temp_s + prcp_s + (1 | species),
          data = dat_r,
          family = "gaussian",
          bayes = FALSE,
          REML = FALSE
        )
        tibble(
          boot = b,
          term = rownames(fit[["B"]]),
          estimate = as.numeric(fit[["B"]][,1]),
          se = as.numeric(fit[["B.se"]][,1]),
          p = as.numeric(fit[["B.pvalue"]][,1]),
          
          delta_AIC = fit_nophy$AIC - fit$AIC,
          delta_logLik = fit$logLik - fit_nophy$logLik,
          delta_BIC = fit_nophy$BIC - fit$BIC
        )
        })
      boot_summary <- boot_fits %>%
        filter(!is.na(term)) %>%
        group_by(term) %>%
        summarise(
          mean = mean(estimate, na.rm = TRUE),
          sd = sd(estimate, na.rm = TRUE),
          q025 = quantile(estimate, 0.025, na.rm = TRUE),
          q50  = quantile(estimate, 0.50,  na.rm = TRUE),
          q975 = quantile(estimate, 0.975, na.rm = TRUE),
          .groups = "drop"
        )
      phy_structure <- boot_fits %>%
        filter(term == "(Intercept)") %>%
        summarise(
          mean_delta_AIC = mean(delta_AIC, na.rm = TRUE),
          mean_delta_logLik = mean(delta_logLik, na.rm = TRUE),
          mean_delta_BIC = mean(delta_BIC, na.rm = TRUE),
          ttest_mean_delta_AIC = t.test(delta_AIC, mu = 0)$p.value,
          ttest_mean_delta_logLik = t.test(delta_logLik, mu = 0)$p.value,
          ttest_mean_delta_BIC = t.test(delta_BIC, mu = 0)$p.value
        )
      model <- boot_summary
      full_data <- boot_fits
      p_SOS_hist <- boot_fits %>%
        mutate(term = dplyr::case_when(
          term == "temp_s"    ~ temp,   # temp = "d_temp"
          term == "prcp_s"    ~ prcp,   # prcp = "d_prcp"
          TRUE ~ term
               )) %>%
        group_by(term) %>%
        mutate(
          lower_95 = quantile(estimate, 0.025),
          upper_95 = quantile(estimate, 0.975),
          Significance = if_else(
            (lower_95 > 0 | upper_95 < 0),
            "Significant",
            "Non-Significant"
          )
        ) %>%
        ungroup() %>%
        ggplot(aes(x = estimate, y = term)) +
        # 2. Re-add the density as a filled geom_density to control the alpha/color better
        geom_density_ridges(
          aes(fill = Significance),
          alpha = 0.5,
          scale = 0.9,
          rel_min_height = 0.01,
          color = "black" # Optional: add a thin black border
        ) +
        
        # 3. Add the Median and 95% CI on top of the density
        stat_pointinterval(
          .width = c(0.95),
          color = "black",
          size = 2,
          point_size = 1.5
        ) +
        scale_fill_manual(values = c(
          "Significant" = "steelblue",
          "Non-Significant" = "lightcoral"
        )) +
        
        # Add the dashed line at zero for reference
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        
        # Set the title and theme
        labs(
          title = paste0(phe_type, ", Fixed effects"),
          x = "Posterior Estimate (95% CI)",
          y = NULL,
          fill = NULL
        ) +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12)
        )
    } else if(bootstrap == FALSE){
      model <- pglmm(
        Estimate ~ 1+ (1|species__) + temp_s + prcp_s,
        data = D2,
        family = "gaussian",
        cov_ranef = list(species = A2),
        bayes = T
      )
      p_SOS_hist <- plot_bayes(model, sort = F)
      full_data <- NA
    }
    
  } else if(package == "brms"){
  # Fit the brms model
    model <- brm(
      # Formula: The syntax (1 | gr(factor, cov = matrix)) specifies a
      # random effect structured by a custom covariance matrix (A2).
      Estimate ~ temp_s + prcp_s + (1 | gr(species, cov = A2)) ,
      data = D2,
      family = gaussian(),
      prior = set_prior("normal(0, 5)", class = "b"),
      chains = 4,
      cores = 4,
      iter = 4000,
      warmup = 2000,
      data2 = list(A2 = A2),
      control = list(adapt_delta = 0.95, max_treedepth = 15)
    )
  
    post_draws <- as_tibble(model)
    p_SOS_hist <- post_draws %>%
      select(b_Intercept, b_temp_s, b_prcp_s) %>%
      pivot_longer(
        cols = starts_with("b_"),
        names_to = "Parameter",
        values_to = "Estimate"
      ) %>%
      mutate(Parameter = gsub("b_", "", Parameter),
             Parameter = dplyr::case_when(
               Parameter == "Intercept" ~ "Intercept",
               Parameter == "temp_s"    ~ temp,   # temp = "d_temp"
               Parameter == "prcp_s"    ~ prcp,   # prcp = "d_prcp"
               TRUE ~ Parameter
             )) %>%
      group_by(Parameter) %>%
      mutate(
        lower_95 = quantile(Estimate, 0.025),
        upper_95 = quantile(Estimate, 0.975),
        Significance = if_else(
          (lower_95 > 0 | upper_95 < 0),
          "Significant",
          "Non-Significant"
        )
      ) %>%
      ungroup() %>%
      ggplot(aes(x = Estimate, y = Parameter)) +
      # 2. Re-add the density as a filled geom_density to control the alpha/color better
      geom_density_ridges(
        aes(fill = Significance),
        alpha = 0.5,
        scale = 0.9,
        rel_min_height = 0.01,
        color = "black" # Optional: add a thin black border
      ) +
  
      # 3. Add the Median and 95% CI on top of the density
      stat_pointinterval(
        .width = c(0.95),
        color = "black",
        size = 2,
        point_size = 1.5
      ) +
      scale_fill_manual(values = c(
        "Significant" = "steelblue",
        "Non-Significant" = "lightcoral"
      )) +
  
      # Add the dashed line at zero for reference
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
      # Set the title and theme
      labs(
        title = paste0(phe_type, ", Fixed effects"),
        x = "Posterior Estimate (95% CI)",
        y = NULL,
        fill = NULL
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)
      )
    full_data = post_draws
  }
  return(list(model = model, plot = p_SOS_hist, full_data = full_data, phy_structure = phy_structure) )
}
