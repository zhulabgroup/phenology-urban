library(scales)
plot_city_coef <- function(
    coef_summary_all,
    level = "species",
    city_order_bio = "bio1" 
) {
  
  # ------------------------------------------------------------
  # 1. Select the corresponding phylogenetic tree
  # ------------------------------------------------------------
  class_tree <- switch(
    level,
    species = phylo_tree,
    genus   = tree_genus,
    family  = tree_family,
    stop("level must be 'species', 'genus', or 'family'")
  )
  
  # ------------------------------------------------------------
  # 2. Prepare coefficients
  # ------------------------------------------------------------
  df <- coef_summary_all |>
    group_by(
      city_label,
      .data[[level]],
      phe_type
    ) |>
    summarise(
      Estimate = mean(Estimate, na.rm = TRUE),
      
      # Keep the existing significance aesthetics
      Sig = first(Sig),
      SigAlpha = first(SigAlpha),
      Color = first(Color),
      
      .groups = "drop"
    ) |>
    left_join(
      city_feature,
      by = "city_label"
    ) |>
    mutate(
      city = case_when(
        city == "DT" ~ "Detroit",
        city == "TP" ~ "Tampa",
        TRUE ~ city
      )
    )
  
  # ------------------------------------------------------------
  # 3. Keep taxa included in the corresponding tree
  # ------------------------------------------------------------
  df_filtered <- df
  
  # ------------------------------------------------------------
  # 4. Order cities by the selected bioclimatic variable
  # ------------------------------------------------------------
  city_order <- df_filtered |>
    distinct(
      city_label,
      .data[[city_order_bio]]
    ) |>
    filter(
      !is.na(.data[[city_order_bio]])
    ) |>
    arrange(.data[[city_order_bio]]) |>
    pull(city_label)
  
  coef_df_plot <- df_filtered |>
    mutate(
      city_label = factor(city_label, levels = city_order)
    )
  
  level_title <- case_when(
    level == "species" ~ "species",
    level == "genus"   ~ "genus",
    level == "family"  ~ "family"
  )
  
  facet_limits <- coef_df_plot |>
    summarise(
      max_abs = max(abs(Estimate), na.rm = TRUE),
      .by = phe_type
    ) |>
    tidyr::expand_grid(sign = c(-1, 1)) |>
    mutate(y_limit = sign * max_abs)
  
  # ------------------------------------------------------------
  # 5. Plot
  # ------------------------------------------------------------
  ggplot(
    coef_df_plot,
    aes(
      x = city_label,
      y = Estimate,
      group = .data[[level]]
    )
  ) +
    geom_blank(
      data = facet_limits,
      aes(y = y_limit),
      inherit.aes = FALSE
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_point(
      aes(
        color = Color,
        alpha = SigAlpha
      ),
      position = position_jitter(
        width = 0.16,
        height = 0,
        seed = 1
      ),
      size = 2
    ) +
    scale_color_identity() +
    scale_alpha_identity() +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 1),
      breaks = scales::breaks_pretty(n = 6),
      limits = c(-20, 20)
    ) +
    labs(
      x = NULL,
      y = bquote(beta[T]~" estimate")
    ) +
    facet_wrap(
      vars(phe_type),
      ncol = 1,
      scales = "free_y",
      labeller = as_labeller(function(x) {
        dplyr::case_when(
          stringr::str_detect(x, "SOS") ~ "a SOS",
          stringr::str_detect(x, "EOS") ~ "b EOS",
          TRUE ~ x
        )
      })
    ) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(
        face = "bold",
        hjust = 0,
        size = 16
      ),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 1
      ),
      plot.title = element_text(
        face = "bold"
      )
    )
}


sig_to_logical <- function(x) {
  if (is.logical(x)) {
    return(replace_na(x, FALSE))
  }
  tolower(trimws(as.character(x))) %in%
    c("*", "significant", "true", "t", "yes", "1")
}

######## summarize

plot_hist_effect_type <- function(
    overlap_city_species_range_model
) {
  effect_type_df <- overlap_city_species_range_model |>
    mutate(
      is_sig = sig_to_logical(Sig),
      effect_type = case_when(
        is_sig  & Estimate > 0 ~ "Significant positive",
        is_sig  & Estimate < 0 ~ "Significant negative",
        !is_sig & Estimate > 0 ~ "Non-significant positive",
        !is_sig & Estimate < 0 ~ "Non-significant negative",
        TRUE                      ~ NA_character_
      ),
      effect_type = factor(
        effect_type,
        levels = c(
          "Non-significant positive",
          "Significant positive",
          "Non-significant negative",
          "Significant negative"
        )
      ),
      SigAlpha = as.numeric(SigAlpha)
    ) |>
    filter(!is.na(effect_type)) |>
    select(
      Estimate,
      phe_type,
      effect_type,
      Color,
      SigAlpha
    )
  
  effect_aes <- effect_type_df |>
    group_by(effect_type) |>
    summarise(
      Color = first(Color),
      SigAlpha = first(SigAlpha),
      .groups = "drop"
    )
  fill_values <- effect_aes |>
    select(effect_type, Color) |>
    deframe()
  alpha_values <- effect_aes |>
    select(effect_type, SigAlpha) |>
    deframe()
  facet_limits <- effect_type_df |>
    summarise(
      max_abs = max(abs(Estimate), na.rm = TRUE),
      .by = phe_type
    ) |>
    tidyr::expand_grid(sign = c(-1, 1)) |>
    mutate(y_limit = sign * max_abs)
  
  p_effect_hist <- ggplot(effect_type_df,
    aes(
      y = Estimate,
      fill = effect_type,
      group = effect_type
    )
  ) +
    geom_histogram(
      aes(alpha = effect_type),
      bins = 150,
      linewidth = 0.4,
      position = "stack"
    ) +
    geom_blank(
      data = facet_limits,
      aes(y = y_limit),
      inherit.aes = FALSE
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    scale_fill_manual(
      values = fill_values,
      drop = FALSE,
      name = NULL
    ) +
    scale_alpha_manual(
      values = alpha_values,
      drop = FALSE,
      guide = "none"
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 1),
      breaks = scales::breaks_pretty(n = 6),
      limits = c(-20, 20)
    ) +
    scale_x_continuous(
      breaks = c(0, 10, 100, 200),
      trans = "sqrt"
    ) +
    labs(
      x = NULL, 
      y = NULL
    ) +
    facet_wrap(
      vars(phe_type),
      ncol = 1,
      scales = "free_y"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.x = element_text(
        color = "white",
        face = "bold",
        hjust = 0,
        size = 16
      ),
      plot.title = element_text(
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 1
      ),
      axis.title.x = element_blank(),
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank(),
      # axis.line.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  return(p_effect_hist)
}
