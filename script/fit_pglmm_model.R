fit_pglmm <- function(phe_type, coef_summary_all) {
  D <- coef_summary_all %>%
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
      bio12_s  = as.numeric(scale(bio12)),
      bio1_s   = as.numeric(scale(bio1)),
      w        = 1/(sd^2)
    ) %>%
    droplevels()
  
  A2 <- A2[levels(D2$species), levels(D2$species)]
  
  # 2. Fit the brms model
  m_brms <- brm(
    # Formula: The syntax (1 | gr(factor, cov = matrix)) specifies a
    # random effect structured by a custom covariance matrix (A2).
    Estimate ~ bio1_s + bio12_s + (1 | gr(species, cov = A2)),
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
  
  post_draws <- as_tibble(m_brms)
  p_SOS_hist <- post_draws %>%
    select(b_Intercept, b_bio1_s, b_bio12_s) %>%
    pivot_longer(
      cols = starts_with("b_"),
      names_to = "Parameter",
      values_to = "Estimate"
    ) %>%
    mutate(Parameter = gsub("b_", "", Parameter)) %>%
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
  return(list(model = m_brms, plot = p_SOS_hist))
}